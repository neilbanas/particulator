classdef modelRun_nemoGlobal2d < modelRun

	% for working with output from NOC Southampton NEMO/MEDUSA global model.
	% so far works on test files from Andy Yool, Nov 2017, from a 1/12 deg
	% hindcast
	%
	% assumes vertical averaging over some depth range: doesn't implement
	% vertical advection or diffusion.	

	properties
		maskFile, gridFile, filenames, filedir;
		
		region					% row and column bounds of the subset of the
								% grid to use. These aren't indices into the
								% full 4322 x 3059 global grid, but rather the
								% 3522 x 1189 subset that Andy Yool passed on
								% to Neil in Jan 2018: north of 30N only, and
								% with a chunk of Eurasian land omitted
								
		oneFrameOverAndOver		% flag making it possible to use the
								% annual average Andy sent as a test file
								% as if it were a full model run
								
		box						% start, count vectors to make it easy to
								% work with only a subset of the global grid
		
		F0, F1					% two frames of in-memory storage

		vars					% lookup table for associating filenames with
								% standard variable names
		
		avg						% setup for depth-averaging

		si						% scatteredInterpolant objects for fields
								% that don't change (H,mask) 
	end
	
	
	methods
	
		function run = modelRun_nemoGlobal2d(filenameTemplate,depthRange, ...
											 region, griddir);
			if nargin < 4 || isempty(griddir)
				griddir = '/Users/neil/Dropbox/particulator/data/nemo-global/';
			end
			run.maskFile = [griddir 'mask3D_30N.mat'];
			run.gridFile = [griddir 'mesh_hgr_withholes_AXY.nc'];
			% filenameTemplate is something like
			% .../nemo-global/2009/banas_ORCA0083-N06_[date]d05[var].nc

			run.oneFrameOverAndOver = 0;
			
			if nargin < 3 || isempty(region)
				region = [1 3522 1 1189]; % matlab indices (1-based),
										  % not netcdf indices!
			end
			% 'region' can be a 4-element vector containing the row and column
			% bounds of a subset of the grid, or else one of these named
			% presets
			if strcmpi(region,'pan-arctic')
				run.region = [300 3350 300 1189];
			elseif strcmpi(region,'bering-chukchi')
				run.region = [600 1200 300 1000];
			elseif strcmpi(region,'atl-arctic')
				run.region = [1900 3350 300 1189]; 
			else
				run.region = region;
			end
			
				%internal name, nc file, name inside nc file
			tab = {...
				'u',	'U',	'uo'; ...		% m/s
				'v',	'V',	'vo'; ...
				'Ks',	'W',	'difvho'; ...	% m2/s
				'temp',	'T',	'potemp'; ...
				'salt',	'T',	'salin'; ...
				'ice',	'I',	'ice_pres'; ... % fractional ice cover
				'iceh',	'I',	'sit'; ...		% ice thickness
				'snow',	'I',	'snd'; ... 		% snow thickness
				'chl',	'P',	'CHL'; ... 		% chlorophyll, mg chl/m3
				'P',	'P',	'PHY'; ... 		% phytop carbon, mmolC/m3
				'NO3',	'P',	'DIN'; ...		% mmol N/m3
				'MZ',	'P',	'ZMC'};			% microzooplankton carbon
			run.vars.local = tab(:,1);
			run.vars.ncfile = tab(:,2);
			run.vars.ncname = tab(:,3);
			
			% assemble lists of filenames for the *U.nc, *V.nc, etc files.
			% in the 2009 dataset we're starting with, the dates on the *P.nc
			% files are off by one day from the others--the difference in
			% timebase is ignored. But this issue is why we're assembling
			% full lists of filenames for each variable rather than assuming
			% they follow a template.
			i_dateBit = strfind(filenameTemplate,'[date]');
			i_middleBit = i_dateBit + 6;
			i_varBit = strfind(filenameTemplate,'[var]');
			i_endBit = i_varBit+5;
			i_dirs = strfind(filenameTemplate,'/');
			run.filedir = filenameTemplate(1:i_dirs(end));
			firstBit = filenameTemplate(1:i_dateBit-1);
			middleBit = filenameTemplate(i_middleBit:i_varBit-1);
			endBit = filenameTemplate(i_endBit:end);
			theVars = unique(run.vars.ncfile);
			for i=1:length(theVars)
				theFiles = dir([firstBit '*' middleBit theVars{i} endBit]);
				run.filenames.(theVars{i}) = {theFiles(:).name};
			end
			
			% read tracer grid + timebase
			ncname = [run.filedir run.filenames.T{1}];
			nc = netcdf.open(ncname,'NOWRITE'); % first tracer file
			run.box.start = [run.region(1) run.region(3)] - [1 1];
			run.box.count = [run.region(2) run.region(4)] - run.box.start;
			grid.x = netcdf.getVar(nc,netcdf.inqVarID(nc,'nav_lon'),...
				run.box.start, run.box.count, 'double');
			grid.x = mod(grid.x, 360); % switch from -180..180 to 0..360
			grid.y = netcdf.getVar(nc,netcdf.inqVarID(nc,'nav_lat'),...
				run.box.start, run.box.count, 'double');
			grid.z = - netcdf.getVar(nc,netcdf.inqVarID(nc,'deptht'),'double');			
			t = netcdf.getVar(nc,netcdf.inqVarID(nc,'time_counter'),'double');
			run.t = t./86400 + datenum('1/1/1950');
			netcdf.close(nc);	
			if run.oneFrameOverAndOver
				run.t = run.t + (0 : 5 : 365); % just pretend
				run.numFrames = length(run.t);
			else
				dt = str2num(middleBit(2:end)); % interval between files
				run.t = run.t + dt .* (0:length(run.filenames.T)-1);
				run.numFrames = length(run.t);
			end
			% read w grid		
			ncname = [run.filedir run.filenames.W{1}];
			nc = netcdf.open(ncname,'NOWRITE'); % first w file
			grid.zw = - netcdf.getVar(nc,netcdf.inqVarID(nc,'depthw'),'double');			
			netcdf.close(nc);
			[I,J] = size(grid.x); % these match run.region
			K = length(grid.z);
			% read 3D mask and convert into 2D mask and H
			load(run.maskFile,'mask3');
			mask3 = mask3(run.region(1):run.region(2), ...
						  run.region(3):run.region(4), :);
			mask3 = double(mask3);
			dz3 = - diff(grid.zw);
			dz3 = dz3([1:end end]);
			dz3 = repmat(reshape(dz3,[1 1 K]),[I J 1]);
			grid.H = sum(dz3.*mask3,3); % convert the 3D mask into a
										% 2D water depth variable
			grid.mask = mask3(:,:,1); % surface layer of the 3D mask = 
									  % the 2D land mask	
			% dz3, mask3 will be used again in the run.avg setup, but not 
			% otherwise saved		
					
			run.grid = grid;
			
			% make scatteredInterpolants for variables that don't change
			warning off
			run.si.H = scatteredInterpolant(grid.x(:),grid.y(:),grid.H(:));
			run.si.mask = ...
				scatteredInterpolant(grid.x(:),grid.y(:),grid.mask(:));
			[jj,ii] = meshgrid(1:J,1:I);
			run.si.ii = scatteredInterpolant(grid.x(:),grid.y(:),ii(:));
			run.si.jj = scatteredInterpolant(grid.x(:),grid.y(:),jj(:));
				% i is a row index and longitude-like
				% j is a column index and latitude-like 
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
			kk = run.avg.kRange(1) : run.avg.kRange(2);
			run.avg.mask = zeros(I,J,K);
			run.avg.mask(:,:,kk) = mask3(:,:,kk);
			run.avg.dz = dz3;
			run.avg.dz(run.avg.mask==0) = 0;
			run.avg.h = sum(run.avg.dz,3);

		end % constructor
		
		
		% reading from model files ---------------------------------------------
		
		
		function c = read(run,localVarname,n);
			vi = strmatch(localVarname,run.vars.local,'exact');
			if ~isempty(vi) 
				ncfile = run.vars.ncfile{vi};
				ncname = [run.filedir run.filenames.(ncfile){n}];
				nc = netcdf.open(ncname,'NOWRITE');
				disp(['reading ' run.vars.ncname{vi} ' from ' ncname]);
				varid = netcdf.inqVarID(nc,run.vars.ncname{vi});
				[~,~,dimids,~] = netcdf.inqVar(nc,varid);
				if length(dimids)==4 % 3D x time
					K = length(run.grid.z);
					c = netcdf.getVar(nc,varid, [run.box.start 0 0], ...
									[run.box.count K 1], 'double');
				else % 2D x time
					c = netcdf.getVar(nc,varid, [run.box.start 0], ...
									[run.box.count 1], 'double');
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
			if run.oneFrameOverAndOver, n=1; end
			% load U, V
			run.F1.u = run.read('u',n);
			run.F1.v = run.read('v',n);
			% load Ks
			run.F1.Ks = run.read('Ks',n);
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
						% this isn't correct for Ks, which is defined on
						% the w grid. But depth averages of Ks don't
						% mean much anyway (interpProfile is where this
						% issue needs to be dealt with better)
					end
					% now make the scatteredInterpolant
					warning off
					run.F1.si.(fields{i}) = scatteredInterpolant( ...
						run.grid.x(:), run.grid.y(:), C(:));
					warning on
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
			v_axis = run.grid.z; % not correct for vertical diffusivity,
								 % which is on the w grid
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




