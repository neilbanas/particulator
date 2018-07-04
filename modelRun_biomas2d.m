classdef modelRun_biomas2d < modelRun

	% for working with output from the Biomas model.
	% pre-averages over a specified depth range, so that all interpolations
	% are 2D.

	properties
		F0, F1					% two frames of in-memory storage
		year
		griddir = 'data/biomas/';
		dirname, basename
		
		localVars, fileVars		% lookup table for associating filenames with
								% standard variable names
		tracerDims				% are the named variables 2D or 3D
		
		avg, avgu				% setup for depth-averaging
		pad						% setup for padding grid with strips at x=0,360
								% for wraparound interpolation
		si						% scatteredInterpolant objects for fields
								% that don't change (H,mask) 
	end
	
	methods
	
		function run = modelRun_biomas2d(dirname,year,griddir,depthRange);
			run.year = year;
			run.numFrames = 365;
			run.t = datenum(year,0,0) + (1:365);
			
			run.dirname = dirname;
			run.basename = ['_600_300.H' num2str(run.year)];
			run.localVars = ...
				{'uv', 'w',    'Ks',    'temp', 'ice',  'iceh', 'swrad'};
			run.fileVars = ...
				{'uo', 'woday','vdcday','to',  'aiday', 'hiday','osswday'};
			run.tracerDims = ...
				[ 3     3       3        3      2        2       2]; % 2D or 3D
				% physical variables have been renamed according to a personal,
				% ROMS-ish convention.
				%	ice = fractional ice cover
				%	iceh = ice thickness
				% 	swrad = shortwave radiation, of which PAR is 0.43
				% any var not in this list (like all the bio vars) is presumed
				% to be in a file with the same prefix as the var name--and also
				% presumed to be 3D.
			
			% load grid.
			% for now the grid is stored with the particulator code, not with 
			% each run itself: hardwired to the 600x300x40 configuration
			if nargin > 2
				if ~isempty(griddir)
					run.griddir = griddir;
				end
			end
			NN = [600 300];
			% scalar grid x,y
			% a column is a bit like a longitude line and a row is a bit like
			% a latitude line, but the grid is stretched and the pole is over
			% land in Alaska
			a = load('-ascii',[run.griddir 'grid.dat.rot']);
			a = a';
			grid.x = reshape(a(1:end/2),NN);
			grid.y = reshape(a(end/2+1:end),NN);
			% velocity grid xu,yu and cell edge lengths
			a = load('-ascii',[run.griddir 'grid.dat.pop']);
			a = a';
			a = a(:);
			N = prod(NN);
			grid.yu = reshape(a(1:N),NN);
			grid.xu = reshape(a(N+1:2*N),NN);
			% lengths of grid cell edges in km
			grid.hun = reshape(a(2*N+1:3*N),NN); 
			grid.hue = reshape(a(3*N+1:4*N),NN);
			grid.hus = reshape(a(4*N+1:5*N),NN);
			grid.huw = reshape(a(5*N+1:6*N),NN);
			% layer depths
			grid.dz = load('-ascii',[run.griddir 'dz.dta40']);
			grid.dz = grid.dz./100; % cm -> m
			grid.zw = -cumsum(grid.dz(:));
			grid.z = grid.zw + grid.dz./2;
			% depths and land mask on the tracer grid
			fid = fopen([run.griddir 'levels_40_t_aBering1']);
			a = fread(fid);
			fclose(fid);
			a = a(a~=10);
			a = reshape(a,[2 length(a)/2])';
			a(a==' ') = '0';
			a = a - '0';
			k = 10.*a(:,1) + a(:,2);
			grid.mask_rho = reshape(k > 0,NN);
			grid.H = zeros(NN);
			grid.H(grid.mask_rho) = -grid.zw(k(grid.mask_rho));
			% depths and velocity-based land mask on the u grid.
			% this is more conservative (classifies more coastal points as land)
			load([griddir 'velocityBasedMask.mat'],'masku');
			grid.masku = masku;
			grid.Hu = griddata(grid.x,grid.y,grid.H,grid.xu,grid.yu);
			% for consistency's sake, we need a final mask variable called
			% "mask" that behaves well for both tracers and velocity--this is
			% what gets used by rel.avoidLand, for example
			% the naming would get confusing if this were on the u grid,
			% so it's on the tracer grid, even though it is controlled almost
			% entirely by the velocity-based mask
			grid.mask = double(grid.mask_rho==1 ...
				& griddata(grid.xu,grid.yu,grid.masku,grid.x,grid.y)>0.5);
			
			% the only bound where it's possible to be out of bounds
			grid.ymin = min(grid.y(:));
			
			run.grid = grid;
			
			% setup for depth averaging
			if length(depthRange)==1
				depthRange = depthRange.*[1 1];
			end
			k = find(abs(grid.zw-depthRange(1)) == ...
					 min(abs(grid.zw-depthRange(1))));
			k2 = find(abs(grid.zw-depthRange(2)) == ...
					 min(abs(grid.zw-depthRange(2))));
			run.avg.kRange = sort([k(1) k2(1)]);
			run.avg.zRange = run.grid.zw(run.avg.kRange);
			[I,J] = size(run.grid.x);
			K = length(run.grid.zw);
			zw3 = repmat(reshape(run.grid.zw,[1 1 K]),[I J 1]);
			dz3 = repmat(reshape(run.grid.dz,[1 1 K]),[I J 1]);
			% on the rho/tracer grid
			H3 = repmat(run.grid.H,[1 1 K]);
			run.avg.mask = zeros(I,J,K);
			run.avg.mask(:,:,run.avg.kRange(1):run.avg.kRange(2)) = 1;
			run.avg.mask(zw3 < -H3) = 0;
			run.avg.dz = dz3;
			run.avg.dz(run.avg.mask==0) = 0;
			run.avg.h = sum(run.avg.dz,3);
			% and again on the u/velocity grid
			Hu3 = repmat(run.grid.Hu,[1 1 K]);
			run.avgu.kRange = run.avg.kRange;
			run.avgu.zRange = run.avg.zRange;
			run.avgu.mask = zeros(I,J,K);
			run.avgu.mask(:,:,run.avg.kRange(1):run.avg.kRange(2)) = 1;
			run.avgu.mask(zw3 < -Hu3) = 0;
			run.avgu.dz = dz3;
			run.avgu.dz(run.avgu.mask==0) = 0;
			run.avgu.h = sum(run.avgu.dz,3);

			% setup for padding fields at x~0 and x~360
			padWidth = 1; % deg longitude
			fnear0 = find(grid.x < padWidth);
			fnear360 = find(grid.x > 360-padWidth);
			ffull = 1:prod(size(grid.x));
			run.pad.ind = [ffull(:); fnear0(:); fnear360(:)];
			run.pad.x = [run.grid.x(:); run.grid.x(fnear0) + 360; ...
						 run.grid.x(fnear360) - 360];
			run.pad.y = run.grid.y(run.pad.ind);
			fnear0 = find(grid.x < padWidth);
			fnear360 = find(grid.x > 360-padWidth);
			run.pad.indu = [ffull(:); fnear0(:); fnear360(:)];
			run.pad.xu = [run.grid.xu(:); run.grid.xu(fnear0) + 360; ...
						 run.grid.xu(fnear360) - 360];
			run.pad.yu = run.grid.y(run.pad.indu);
			
			% scatteredInterpolants for fields that don't change in time
			run.si.H = scatteredInterpolant(...
				run.pad.x(:),run.pad.y(:),run.grid.H(run.pad.ind));
			run.si.mask = scatteredInterpolant(...
				run.pad.x(:),run.pad.y(:),double(run.grid.mask(run.pad.ind)));
				
			% scatteredInterpolants for (i,j) indices in the horizontal. This
			% is setup for interpProfile()
			[ii,jj] = meshgrid(1:NN(2),1:NN(1));
			run.si.ii = scatteredInterpolant(...
				run.pad.x(:),run.pad.y(:),ii(run.pad.ind));
			run.si.jj = scatteredInterpolant(...
				run.pad.x(:),run.pad.y(:),jj(run.pad.ind));
		end % constructor
		
		
		% reading from model files ---------------------------------------------
		
		
		function fid = openFile(run,localVar);
			[I,J] = size(run.grid.x);
			K = length(run.grid.dz);
			framelength = I * J * K * 4;
			
			j = find(strcmp(localVar,run.localVars));
			prefix = [run.fileVars{j}];
			if isempty(prefix)
				prefix = localVar;
			end
			
			fid = fopen([run.dirname prefix run.basename]);
		end
		
		
		function C = read3D(run,localVar,n);
			[I,J] = size(run.grid.x);
			K = length(run.grid.dz);
			framelength = I * J * K * 4;
				% bytes per 600x300x40 frame of one variable
			fid = run.openFile(localVar);			
			fseek(fid,framelength*(n-1),-1);
			C = reshape(fread(fid,I*J*K,'real*4'),[I J K]);
			C(~isfinite(C)) = 0; % this isn't really right for all variables
			fclose(fid);
		end
		
		
		function C = read2D(run,localVar,n);
			[I,J] = size(run.grid.x);
			framelength = I * J * 4;
				% bytes per 600x300 frame of one variable
			fid = run.openFile(localVar);			
			fseek(fid,framelength*(n-1),-1);
			C = reshape(fread(fid,I*J,'real*4'),[I J]);
			C(~isfinite(C)) = 0;
			fclose(fid);
		end
				
				
		% consider wrapping each variable around the first dimension,
		% ([1:end 1],:), in order to make things interpolate correctly at
		% x ~ 210 in the N Pacific.

		function run = loadFrame(run,n,tracers);
			% u and v
			run.F1.u = run.read3D('uv',2*n-1);
			run.F1.v = run.read3D('uv',2*n); % because the two variables
										     % are striped in one file
			% Ks
			run.F1.Ks = run.read3D('Ks',n);	% assumes Ks is on the tracer grid,
									        % not the w grid
			% everything else
			for i=1:length(tracers)
				j = find(strcmp(tracers{i},run.localVars));
				if ~isempty(j) && run.tracerDims(j)==2
					run.F1.(tracers{i}) = run.read2D(tracers{i},n);
				else
					run.F1.(tracers{i}) = run.read3D(tracers{i},n);
				end
			end

			% create scatteredInterpolant objects for depth averages of
			% all fields. This makes repeated calls to interp() much faster.
			fields = fieldnames(run.F1);
			for i=1:length(fields)
				if isnumeric(run.F1.(fields{i}))
					% depth-average if necessary
					if ndims(run.F1.(fields{i})) == 2
						C = run.F1.(fields{i});
					elseif strcmpi(fields{i},'u') || strcmpi(fields{i},'v')
						C = sum(run.F1.(fields{i}) .* run.avgu.dz, 3) ...
							./ run.avgu.h;
					else
						C = sum(run.F1.(fields{i}) .* run.avg.dz, 3) ...
							./ run.avg.h;
					end
					% now make the scatteredInterpolant
					if strcmpi(fields{i},'u') || strcmpi(fields{i},'v')
						run.F1.si.(fields{i}) = scatteredInterpolant(...
							run.pad.xu, run.pad.yu, C(run.pad.indu));
					else
						run.F1.si.(fields{i}) = scatteredInterpolant(...
							run.pad.x, run.pad.y, C(run.pad.ind));
					end
				end
			end	

			% declare this frame loaded
			run.loadedN(2) = n;
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
			v_axis = run.grid.z; % ignoring everything on the w grid
		end
		
		function c = interpProfile(run,name,x,y,t);
			% find fractional indices of all the (x,y) points
			ii = run.si.ii(x,y);
			jj = run.si.jj(x,y);
			
			% find the indices of the 2x2 columns of points that surround
			% each (x,y) 
			[J,I] = size(run.grid.x);
			K = size(run.F1.u,3);
			ii0 = min(I-1, max(1, floor(ii)));
			a = repmat(ii - ii0, [1 K]);
			ii0 = repmat(ii0(:), [1 K]);
			jj0 = min(J-1, max(1, floor(jj)));
			b = repmat(jj - jj0, [1 K]);
			jj0 = repmat(jj0(:), [1 K]);
			kk = repmat((1:K), [length(x) 1]);
			ind00 = sub2ind([J I K], jj0,   ii0,   kk);
			ind01 = sub2ind([J I K], jj0,   ii0+1, kk);
			ind10 = sub2ind([J I K], jj0+1, ii0,   kk);
			ind11 = sub2ind([J I K], jj0+1, ii0+1, kk);
			c_n0 = (1-a) .* (1-b).* run.F0.(name)(ind00) ...
				 + (1-a) .*    b .* run.F0.(name)(ind10) ...
				 + a     .* (1-b).* run.F0.(name)(ind01) ...
				 + a     .*    b .* run.F0.(name)(ind11);
			c_n1 = (1-a) .* (1-b).* run.F1.(name)(ind00) ...
				 + (1-a) .*    b .* run.F1.(name)(ind10) ...
				 + a     .* (1-b).* run.F1.(name)(ind01) ...
				 + a     .*    b .* run.F1.(name)(ind11);
			c = run.tinterp(t,c_n0,c_n1);
		end
		
		
		function us = scaleU(run,u,x,y); % cm/s -> deg lon per day
			us = u ./ 100 .* 86400 ./ 111325 ./ cos(y./180.*pi);
		end
		function vs = scaleV(run,v,x,y); % cm/s -> deg lat per day
			vs = v ./ 100 .* 86400 ./ 111325;
		end
		function ws = scaleW(run,w); % cm/s -> m/day
			ws = w ./ 100 ./ 86400;
		end
		
		function [x1,y1,active] = filterCoordinates(run,x,y);
			active = y > run.grid.ymin;
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




