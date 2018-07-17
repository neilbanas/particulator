classdef modelRun_sinmod2d < modelRun

	% for working with output from SINMOD (Slagstad et al. 2015).
	%
	% assumes vertical averaging over some depth range: doesn't implement
	% vertical advection or diffusion.	
	%
	% similar to modelRun_biomas2d and modelRun_nemoGlobal2d, which are also
	% curvilinear, z-level models.

	properties
		filename
		
		F0, F1					% two frames of in-memory storage

		vars					% lookup table for associating filenames with
								% standard variable names
		
		avg						% setup for depth-averaging

		si						% scatteredInterpolant objects for fields
								% that don't change (H,mask) 
	end
	
	
	methods
	
		function run = modelRun_sinmod2d(filename,depthRange);
			%internal name, name inside nc file
			tab = {...
				'u',	'v_north'; ...	% note that u,v are switched!
				'v',	'u_east'; ...
				'Ks',	'????'; ...
				'temp',	'temperature'; ...
				'salt',	'salinity'; ...
				'ice',	'ice_compactness'; ...	% fractional ice cover
				'iceh',	'ice_thickness'};	 	% ice thickness
			run.vars.local = tab(:,1);
			run.vars.ncname = tab(:,2);

			run.filename = filename;			
			nc = netcdf.open(filename,'NOWRITE');
			
			% read grid
			grid.x = netcdf.getVar(nc,netcdf.inqVarID(nc,'gridLons'),'double');
			grid.y = netcdf.getVar(nc,netcdf.inqVarID(nc,'gridLats'),'double');
			grid.zw = - netcdf.getVar(nc,netcdf.inqVarID(nc,'zc'),'double');
			grid.zw = [0; grid.zw];
			grid.dz = - diff(grid.zw);
			grid.z = 0.5.*(grid.zw(1:end-1) + grid.zw(2:end));
			grid.H = netcdf.getVar(nc,netcdf.inqVarID(nc,'depth'),'double');
			grid.mask = double(grid.H < 1e10);
			grid.H(grid.mask==0) = 0;
			[I,J] = size(grid.x); % I, J match the dimensions called xc, yc
								  % inside the netcdf file
			K = length(grid.z);
			dz3 = repmat(reshape(grid.dz,[1 1 K]),[I J 1]);
			z3 = repmat(reshape(grid.z,[1 1 K]),[I J 1]);
			H3 = repmat(grid.H,[1 1 K]);
			mask3 = z3 > -H3;
			
			% read timebase
			t = netcdf.getVar(nc,netcdf.inqVarID(nc,'time'),'double');
			units = netcdf.getAtt(nc,netcdf.inqVarID(nc,'time'),'units');
%			timeref = strrep(units, 'seconds since ', '');
%			timeref = datenum(timeref,'yyyy-mm-dd HH:MM:SS');
%			run.t = linspace(t0,t1,run.numFrames) ./ 86400 + timeref;
			run.t = t - datenum(2006,4,1);
			netcdf.close(nc);

			run.grid = grid;
			
			% make scatteredInterpolants for variables that don't change
			warning off
			run.si.H = scatteredInterpolant(grid.x(:),grid.y(:),grid.H(:));
			run.si.mask = ...
				scatteredInterpolant(grid.x(:),grid.y(:),grid.mask(:));
			[jj,ii] = meshgrid(1:J,1:I);
			run.si.ii = scatteredInterpolant(grid.x(:),grid.y(:),ii(:));
			run.si.jj = scatteredInterpolant(grid.x(:),grid.y(:),jj(:));
			warning on		
			
			% setup for depth averaging
			if length(depthRange)==1
				depthRange = depthRange.*[1 1];
			end
			k = find(abs(grid.zw-depthRange(1)) == ...
					 min(abs(grid.zw-depthRange(1))));
			k2 = find(abs(grid.zw-depthRange(2)) == ...
					 min(abs(grid.zw-depthRange(2))));
			run.avg.kRange = sort([k(1) k2(1)]);
			kk = run.avg.kRange(1) : run.avg.kRange(2) - 1;
				% the - 1 is for matching a zw cell bottom to a z cell center
			run.avg.mask = zeros(I,J,K);
			run.avg.mask(:,:,kk) = mask3(:,:,kk);
			run.avg.dz = dz3;
			run.avg.dz(run.avg.mask==0) = 0;
			run.avg.h = sum(run.avg.dz,3);
		end % constructor
		
		
		% reading from model files ---------------------------------------------
		
		
		function c = read(run,localVarname,n);
			[I,J] = size(run.grid.x);
			vi = strmatch(localVarname,run.vars.local,'exact');
			if ~isempty(vi) 
				nc = netcdf.open(run.filename,'NOWRITE');
				varid = netcdf.inqVarID(nc,run.vars.ncname{vi});
				[~,~,dimids,~] = netcdf.inqVar(nc,varid);
				if length(dimids)==4 % 3D x time
					K = length(run.grid.z);
					c = netcdf.getVar(nc,varid,[0 0 0 n-1],[I J K 1],'double');
				else % 2D x time
					c = netcdf.getVar(nc,varid,[0 0 n-1],[I J 1],'double');
				end
			else
				warning(['don''t know what file contains ' localVarname '.']);
				c = [];
			end
			netcdf.close(nc);
			c(c>1e16) = nan;
		end
		
		
		function run = loadFrame(run,n,tracers);
			run.loadedN(2) = n;
			% load U, V
			run.F1.u = run.read('u',n);
			run.F1.v = run.read('v',n);
			% load Ks
%			run.F1.Ks = run.read('Ks',n);
			run.F1.Ks = zeros(size(run.F1.u));
			% load tracers
			for m=1:length(tracers)
				run.F1.(tracers{m}) = run.read(tracers{m},n);
			end
			
			% Create scatteredInterpolant objects for depth averages of
			% all fields. This makes repeated calls to interp() much faster.
			% An even faster method would be to include ii, jj in
			% par_integrate.interpEverything() and then work entirely in
			% index units, not x, y. 
			fields = fieldnames(run.F1);
			for i=1:length(fields)
				if isnumeric(run.F1.(fields{i}))
					% depth-average if necessary
					if ndims(run.F1.(fields{i})) == 2
						C = run.F1.(fields{i});
					else
						C = sum(run.F1.(fields{i}) .* run.avg.dz, 3) ...
							./ run.avg.h;
						% this may not be correct for Ks, if it's defined on
						% the w grid. But depth averages of Ks don't
						% mean much anyway (interpProfile is where this
						% issue needs to be dealt with better)
					end
					% now make the scatteredInterpolant
					run.F1.si.(fields{i}) = scatteredInterpolant( ...
						run.grid.x(:), run.grid.y(:), C(:));
				end
			end	
		end
	
		
		function run = advanceTo(run,n,tracers);
			run.F0 = run.F1;
			run.loadedN(1) = run.loadedN(2);
			run.loadFrame(n,tracers);
		end
		

		% interpolating model variables ----------------------------------------
		
			
		function c = interp(run,name,x,y,sigma,t);
			if strcmpi(name,'H')
				c = run.si.H(x,y);
			elseif strcmpi(name,'zeta')
				c = 0;
			elseif strcmpi(name,'mask')
				c = run.si.mask(x,y);
			elseif strcmpi(name,'w')
				c = 0;
			else	
				% note that sigma is ignored in all of these--
				% no 3D interpolations at all!		
				c0 = run.F0.si.(name)(x,y);
				c1 = run.F1.si.(name)(x,y);
				c = run.tinterp(t, c0, c1);
			end
		end
		
		
		function c = interpDepthAverage(run,name,x,y,zMinMax,t);
			c = run.interp(name,x,y,[],t);
		end
		
		
		function v_axis = verticalAxisForProfiles(run);
			v_axis = run.grid.z; % maybe not correct for vertical diffusivity
		end
		
		function c = interpProfile(run,name,x,y,t);
			% make sure the variable is the right shape and size for
			% extracting a profile from
			K = size(run.F1.u,3);
			if size(run.F1.(name),3) ~= K
%				warning([name ' is the wrong size for profiles.']);
				c = [];
				return
			end	
			% find fractional indices of all the (x,y) points
			ii = run.si.ii(x,y);
			jj = run.si.jj(x,y);			
			% find the indices of the 2x2 columns of points that surround
			% each (x,y) 
			[I,J] = size(run.grid.x);
			ii0 = min(I-1, max(1, floor(ii)));
			a = repmat(ii - ii0, [1 K]);
			ii0 = repmat(ii0(:), [1 K]);
			jj0 = min(J-1, max(1, floor(jj)));
			b = repmat(jj - jj0, [1 K]);
			jj0 = repmat(jj0(:), [1 K]);
			kk = repmat((1:K), [length(x) 1]);
			ind00 = sub2ind([I J K], ii0,   jj0,   kk);
			ind01 = sub2ind([I J K], ii0,   jj0+1, kk);
			ind10 = sub2ind([I J K], ii0+1, jj0,   kk);
			ind11 = sub2ind([I J K], ii0+1, jj0+1, kk);
			c_n0 = (1-a) .* (1-b).* run.F0.(name)(ind00) ...
				 + (1-a) .*    b .* run.F0.(name)(ind10) ...
				 + a     .* (1-b).* run.F0.(name)(ind01) ...
				 + a     .*    b .* run.F0.(name)(ind11);
			c_n1 = (1-a) .* (1-b).* run.F1.(name)(ind00) ...
				 + (1-a) .*    b .* run.F1.(name)(ind10) ...
				 + a     .* (1-b).* run.F1.(name)(ind01) ...
				 + a     .*    b .* run.F1.(name)(ind11);
			c = run.tinterp(t,c_n0,c_n1);
			% blank out the points below the seabed, since for some tracers
			% the fill value of 0 is a plausible value
			% use the shallowest depth at any of the 2x2 columns, since
			% otherwise we're allowing interpolation between a real value
			% and the land-mask value
			ind00 = sub2ind([I J], ii0(:,1),   jj0(:,1)  );
			ind01 = sub2ind([I J], ii0(:,1),   jj0(:,1)+1);
			ind10 = sub2ind([I J], ii0(:,1)+1, jj0(:,1)  );
			ind11 = sub2ind([I J], ii0(:,1)+1, jj0(:,1)+1);
			ind = [ind00(:) ind01(:) ind10(:) ind11(:)];
			Hmin = min(run.grid.H(ind),[],2);
			Hmin = repmat(Hmin,[1 K]);
			zz = run.verticalAxisForProfiles;
			zz = repmat(zz(:)',[size(Hmin,1) 1]);
			c(zz < -Hmin) = nan;
		end
		
		
		function us = scaleU(run,u,x,y); % m/s -> deg lon per day
			us = u .* 86400 ./ 111325 ./ cos(y./180.*pi);
		end
		function vs = scaleV(run,v,x,y); % m/s -> deg lat per day
			vs = v .* 86400 ./ 111325;
		end
		
		function [x1,y1,active] = filterCoordinates(run,x,y);
			active = y > 30; % model fields end at 30N
			y1 = y;
			x1 = x;
			% deal with the case where a point goes over the pole (y>90).
			f = find(y1>90);
			y1(f) = 180 - y1(f);
			x1(f) = x1(f) + 180;
			% confine x coordinates to 0...360
			x1 = mod(x1,360);
		end		


		function ci = tinterp(run,ti,c0,c1);
			% given values c0, c1 at the times of frame0 and frame1, returns
			% an interpolated value at a time ti in between
			t0 = run.t(run.loadedN(1));
			t1 = run.t(run.loadedN(2));
			f = (ti-t0) ./ (t1-t0); % fraction of the way from frame0 to frame1
			f = max(min(f,1),0); % note that running over bounds in time is
								 % handled differently from bounds in space
			ci = c0 + (c1-c0) .* f;
		end
		
		
	end % methods

end % classdef




