classdef modelRun_biomas < modelRun

	% for working with output from the Biomas model.

	properties
		F0, F1	 % two frames of in-memory storage
		outOfBoundsValue
		year
		griddir = 'data/biomas/';
		dirname, basename
		localVars, fileVars		% lookup table for associating files with
								% standard variables
		si	% scatteredInterpolant objects for H, mask
	end
	
	methods
	
		function run = modelRun_biomas(dirname,year,griddir);			
			run.year = year;
			run.numFrames = 365;
			run.t = datenum(2009,0,0) + (1:365);
			
			run.dirname = dirname;
			run.basename = ['_600_300.H' num2str(run.year)];
			run.localVars = {'uv','w','Ks','temp','ice','iceh','swrad'};
				% physical variables have been renamed according to a personal,
				% ROMS-ish convention (ice = fractional ice cover, iceh = ice
				% thickness, swrad = shortwave radiation, of which PAR is 0.43).
				% any var not in this list (like all the bio vars) is presumed
				% to be in a file with the same prefix as the var name.
			run.fileVars = ...
				{'uo','woday','vdcday','to','aiday','hiday','osswday'};
			
			run.nativeSigma = 0;
			run.outOfBoundsValue = 0;
			run.wScaleFactor = 86400;
			
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
			% depths and land mask
			fid = fopen([run.griddir 'levels_40_t_aBering1']);
			a = fread(fid);
			fclose(fid);
			a = a(a~=10);
			a = reshape(a,[2 length(a)/2])';
			a(a==' ') = '0';
			a = a - '0';
			k = 10.*a(:,1) + a(:,2);
			grid.mask = reshape(k > 0,NN);
			grid.H = zeros(NN);
			grid.H(grid.mask) = grid.zw(k(grid.mask));

			%{
			% 3D grids for convenience
			[I,J] = size(grid.H);
			K = length(grid.dz);
			% ----- tracer grid
			grid.x3 = repmat(grid.x,[1 1 K+2]);
			grid.y3 = repmat(grid.y,[1 1 K+2]);
			grid.z3 = repmat(....
				reshape([0; grid.z; grid.zw(end)],[1 1 K+2]),[I J 1]);
			grid.H3 = repmat(grid.H,[1 1 K+2]);
			grid.sigma3 = grid.z3 ./ grid.H3;
			grid.sigma3(~isfinite(grid.sigma3)) = 0;
			% ----- u,v grid
			grid.xu3 = repmat(grid.xu,[1 1 K+2]);
			grid.yu3 = repmat(grid.yu,[1 1 K+2]);
			grid.zu3 = repmat(....
				reshape([0; grid.z; grid.zw(end)],[1 1 K+2]),[I J 1]);
			grid.Hu3 = grid.H3; % not sure this is right
			grid.sigmau3 = grid.zu3 ./ grid.Hu3;
			grid.sigmau3(~isfinite(grid.sigmau3)) = 0;
			% ----- w grid
			grid.xw3 = repmat(grid.x,[1 1 K+1]);
			grid.yw3 = repmat(grid.y,[1 1 K+1]);
			grid.zw3 = repmat(reshape([0; grid.zw],[1 1 K+1]),[I J 1]);
			grid.Hw3 = repmat(grid.H,[1 1 K+1]);
			grid.sigmaw3 = grid.zw3 ./ grid.Hw3;
			grid.sigmaw3(~isfinite(grid.sigmaw3)) = 0;
			%}
			
			grid.ze = [0; grid.z; grid.zw(end)];
			grid.zwe = [0; grid.zw];
				% z,zw extrapolated to surface and bottom
			
			% the only bound where it's possible to be out of bounds
			grid.ymin = min(grid.y(:));
			
			% make scatteredInterpolant objects that don't change
			run.si.H = scatteredInterpolant(grid.x(:),grid.y(:),grid.H(:));
			run.si.mask = ...
				scatteredInterpolant(grid.x(:),grid.y(:),double(grid.mask(:)));
			
			run.grid = grid;
		end % constructor
		
		
		% reading from model files ---------------------------------------------
		
		
		function run = loadFrame(run,n,tracers);
			% read one day of data from each file
			[I,J] = size(run.grid.x);
			K = length(run.grid.dz);
			framelength = I * J * K * 4;
				% bytes per 600x300x40 frame of one variable
			% u and v
			fid = fopen([run.dirname ...
					run.fileVars{strmatch('uv',run.localVars,'exact')} ...
					run.basename]);
			fseek(fid,framelength*(n-1),-1);
			run.F1.u = reshape(fread(fid,I*J*K,'real*4'),[I J K]);
			fseek(fid,framelength*(365+n-1),-1);			
			run.F1.v = reshape(fread(fid,I*J*K,'real*4'),[I J K]);
			fclose(fid);
			run.F1.u = run.F1.u(:,:,[1 1:end end]);
			run.F1.v = run.F1.v(:,:,[1 1:end end]);
				% perhaps could wrap these around the first dimension...
				% ([1:end 1],:,[1 1:end end])
				% ... and likewise for all other fields, in order to make
				% things interpolate correctly at x ~ 210 in the N Pacific.
			% w
			fid = fopen([run.dirname ...
					run.fileVars{strmatch('w',run.localVars,'exact')} ...
					run.basename]);
			fseek(fid,framelength*(n-1),-1);
			run.F1.w = reshape(fread(fid,I*J*K,'real*4'),[I J K]);
			fclose(fid);
			run.F1.w = cat(3,zeros(I,J,1),run.F1.w);
			% Ks
			fid = fopen([run.dirname ...
					run.fileVars{strmatch('Ks',run.localVars,'exact')} ...
					run.basename]);
			fseek(fid,framelength*(n-1),-1);
			run.F1.Ks = reshape(fread(fid,I*J*K,'real*4'),[I J K]);
			fclose(fid);
			run.F1.Ks = run.F1.Ks(:,:,[1 1:end end]);
				% is Ks on the tracer grid or the w grid?
				% this version assumes tracer grid
			% everything else
			for i=1:length(tracers)
				prefix = [run.fileVars{...
							strmatch(tracers{i},run.localVars,'exact')}];
				if isempty(prefix), prefix = tracers{i}; end
				fid = fopen([run.dirname prefix run.basename]);
				fseek(fid,framelength*(n-1),-1);
				run.F1.(tracers{i}) = ...
					reshape(fread(fid,I*J*K,'real*4'),[I J K]);
				fclose(fid);
				run.F1.(tracers{i}) = run.F1.(tracers{i})(:,:,[1 1:end end]);
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

		% note: there's a slice down the N Pacific, through the Gulf of Alaska,
		% between grid.x(end,:) and grid.x(1,:) where things may not interpolate
		% correctly.
		% it also seems unlikely that things will interpolate correctly across
		% the dateline.
		% also at the pole.
		% so perhaps a pragmatic solution would be to take all points where the 
		% griddata() calls return NaNs and move them to the nearest non-masked
		% point.
	
		function H = interpH(run,x,y);
			H = run.si.H(x,y);
		end
		
		function zeta = interpZeta(run,x,y,t);
			zeta = 0;
		end
		
		function mask = interpMask(run,x,y,t);
			mask = run.si.mask(x,y);
			% note that t is not used here
		end
		
		function u = interpU(run,x,y,sigma,t);
			error('no support yet for interpolating 3D BIOMAS fields in sigma coordinates.');
		end
			
		function u = interpU_in_z(run,x,y,z,t);
			% no support yet for full 3D z variables. Assumes the zTrapLevel
			% case. Averages z, finds the model level closest to it, and
			% does a 2D interpolation there.
			z = mean(z);
			k = find(abs(run.grid.ze - z)==min(abs(run.grid.ze - z)));
			k = k(1);
			u0 = griddata(run.grid.xu,run.grid.yu,run.F0.u(:,:,k),x,y);
			u1 = griddata(run.grid.xu,run.grid.yu,run.F1.u(:,:,k),x,y);
			u = run.tinterp(t, u0, u1);
		end

		function v = interpV_in_z(run,x,y,z,t);
			z = mean(z);
			k = find(abs(run.grid.ze - z)==min(abs(run.grid.ze - z)));
			k = k(1);
			v0 = griddata(run.grid.xu,run.grid.yu,run.F0.v(:,:,k),x,y);
			v1 = griddata(run.grid.xu,run.grid.yu,run.F1.v(:,:,k),x,y);
			v = run.tinterp(t, v0, v1);
		end

		function w = interpW_in_z(run,x,y,z,t);
			z = mean(z);
			k = find(abs(run.grid.zwe - z)==min(abs(run.grid.zwe - z)));
			k = k(1);
			w0 = griddata(run.grid.x,run.grid.y,run.F0.w(:,:,k),x,y);
			w1 = griddata(run.grid.x,run.grid.y,run.F1.w(:,:,k),x,y);
			w = run.tinterp(t, w0, w1);
		end
		
		function Ks = interpKs_in_z(run,x,y,z,t);
			Ks = interpTracer_in_z(run,'Ks',x,y,z,t);
		end
		
		function c = interpTracer_in_z(run,name,x,y,z,t);
			if ndims(run.F0.(name))==2 % 2d
				k = 1;
			else
				z = mean(z);
				k = find(abs(run.grid.ze - z)==min(abs(run.grid.ze - z)));
				k = k(1);
			end
			c0 = griddata(run.grid.x,run.grid.y,run.F0.(name)(:,:,k),x,y);
			c1 = griddata(run.grid.x,run.grid.y,run.F1.(name)(:,:,k),x,y);
			c = run.tinterp(t, c0, c1);
		end
		
		
		function us = scaleU(run,u,x,y); % cm/s -> deg lon per day
			us = u ./ 100 .* 86400 ./ 111325 ./ cos(y./180.*pi);
		end
		function vs = scaleV(run,v,x,y); % cm/s -> deg lat per day
			vs = v ./ 100 .* 86400 ./ 111325;
		end
		
		
		function [x1,y1,active] = filterCoordinates(run,x,y);
			active = y > run.grid.ymin;
			% confine x coordinates to 0...360
			x1 = mod(x,360);
			% deal with the case where a point goes over the pole (y>90).
			% here's a bad way of dealing with it:
			y1 = min(y,90);
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




