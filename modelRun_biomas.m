classdef modelRun_biomas < modelRun

	% for working with output from the Biomas model.

	properties
		F0, F1	 % two frames of in-memory storage
		outOfBoundsValue;
		year;
		files;	% lookup table for associating files with variables
	end
	
	methods
	
		function run = modelRun_biomas(dirname,year);
			run.year = year;
			run.numFrames = 365;
			run.t = datenum(2009,0,0) + (1:365);
			
			files.vars = {'uv','w','Ks','temp','ice','iceh','swrad'};
				% physical variables have been renamed according to a personal,
				% ROMS-ish convention (ice = fractional ice cover, iceh = ice
				% thickness, swrad = shortwave radiation, of which PAR is 0.43).
				% any var not in this list (like all the bio vars) is presumed
				% to be in a file with the same prefix as the var name.
			files.prefixes = {'uo','woday','vdcday','to','ai','hi','osswday'};
			files.basename = ['_600_300.H' num2str(run.year)];
			
			% load grid.
			% for now the grid is stored with the particulator code, not with 
			% each run itself: hardwired to the 600x300x40 configuration
			griddir = 'data/biomas/';
			NN = [600 300];
			% scalar grid x,y
			% a column is a bit like a longitude line and a row is a bit like
			% a latitude line, but the grid is stretched and the pole is over
			% land in Alaska
			a = load('-ascii',[griddir 'grid.dat.rot']);
			a = a';
			grid.x = reshape(a(1:end/2),NN);
			grid.y = reshape(a(end/2+1:end),NN);
			% velocity grid xu,yu and cell edge lengths
			a = load('-ascii',[griddir 'grid.dat.pop']);
			a = a';
			a = a(:);
			N = prod(NN);
			grid.yu = reshape(a(1:N),NN);
			grid.xu = reshape(a(N+1:2*N),NN);
			grid.xu(grid.xu>180) = grid.xu(grid.xu>180) - 360;
			% lengths of grid cell edges in km
			grid.hun = reshape(a(2*N+1:3*N),NN); 
			grid.hue = reshape(a(3*N+1:4*N),NN);
			grid.hus = reshape(a(4*N+1:5*N),NN);
			grid.huw = reshape(a(5*N+1:6*N),NN);
			% layer depths
			grid.dz = load('-ascii',[griddir 'dz.dta40']);
			grid.dz = grid.dz./100; % cm -> m
			grid.zw = -cumsum(dz(:));
			grid.z = grid.zw + grid.dz./2;
			% depths and land mask
			fid = fopen([griddir 'levels_40_t_aBering1']);
			a = fread(fid);
			fclose(fid);
			a = a(a~=10);
			a = reshape(a,[2 length(a)/2])';
			a(a==' ') = '0';
			a = a - '0';
			k = 10.*a(:,1) + a(:,2);
			grid.mask = reshape(k > 0,NN);
			grid.H = nan.*ones(NN);
			grid.H(grid.mask) = grid.zw(k(grid.mask));

			% 3D grids for convenience
			[I,J] = size(grid.H);
			K = length(grid.dz);
			grid.H3 = repmat(grid.H,[1 1 K]);
			% ----- tracer grid
			grid.x3 = repmat(grid.x,[1 1 K+2]);
			grid.y3 = repmat(grid.y,[1 1 K+2]);
			grid.z3 = repmat(....
				reshape([0 grid.z grid.zw(end)],[1 1 K+2]),[I J 1]);
			grid.sigma3 = grid.z3 ./ grid.H3;
			% ----- u,v grid
			grid.xu3 = repmat(grid.xu,[1 1 K+2]);
			grid.yu3 = repmat(grid.yu,[1 1 K+2]);
			grid.zu3 = repmat(....
				reshape([0 grid.z grid.zw(end)],[1 1 K+2]),[I J 1]);
			grid.sigmau3 = grid.zu3 ./ grid.H3; % not sure this is right
			% ----- w grid
			grid.xw3 = repmat(grid.x,[1 1 K+1]);
			grid.yw3 = repmat(grid.y,[1 1 K+1]);
			grid.zw3 = repmat(reshape([0 grid.zw],[1 1 K+1]),[I J 1]);
			grid.sigmaw3 = grid.zw3 ./ grid.H3;

			% the only bound where it's possible to be out of bounds
			grid.ymin = min(grid.y(:));
			
			run.grid = grid;
		end % constructor
		
		
		% reading from model files ---------------------------------------------
		
		
		function run = loadFrame(run,n,tracers);
			% read one day of data from each file
			[I,J,K] = size(grid.mask3);
			framelength = I * J * K * 4;
				% bytes per 600x300x40 frame of one variable
			% u and v
			fid = fopen([dirname files.prefixes{strmatch('uv',files.vars)} ...
						 files.basename]);
			fseek(fid,framelength*(n-1),-1);
			run.F0.u = reshape(fread(fid,framelength,'real*4'),[I J K]);
			fseek(fid,framelength*(365+n-1),-1);			
			run.F0.v = reshape(fread(fid,framelength,'real*4'),[I J K]);
			fclose(fid);
			run.F0.u = run.F0.u(:,:,[1 1:end end]);
			run.F0.v = run.F0.v(:,:,[1 1:end end]);
				% perhaps could wrap these around the first dimension...
				% ([1:end 1],:,[1 1:end end])
				% ... and likewise for all other fields, in order to make
				% things interpolate correctly at x ~ 210 in the N Pacific.
			% w
			fid = fopen([dirname files.prefixes{strmatch('w',files.vars)} ...
						 files.basename]);
			fseek(fid,framelength*(n-1),-1);
			run.F0.w = reshape(fread(fid,framelength,'real*4'),[I J K]);
			fclose(fid);
			run.F0.w = cat(3,zeros(I,J 1),run.F0.w);
			% Ks
			fid = fopen([dirname files.prefixes{strmatch('Ks',files.vars)} ...
						 files.basename]);
			fseek(fid,framelength*(n-1),-1);
			run.F0.Ks = reshape(fread(fid,framelength,'real*4'),[I J K]);
			fclose(fid);
			run.F0.Ks = run.F0.Ks(:,:,[1 1:end end]);
				% is Ks on the tracer grid or the w grid?
				% this version assumes tracer grid
			% everything else
			for i=1:length(tracers)
				prefix = [files.prefixes{strmatch(tracers{i},files.vars)}];
				if isempty(prefix), prefix = tracers{i}; end
				fid = fopen([dirname prefix files.basename]);
				fseek(fid,framelength*(n-1),-1);
				run.F0.(tracers{i}) = ...
					reshape(fread(fid,framelength,'real*4'),[I J K]);
				fclose(fid);
				run.F0.(tracers{i}) = run.F0.(tracers{i})(:,:,[1 1:end end]);
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
			H(isin) = griddata(grid.x,grid.y,grid.H,x,y);
		end
		
		function zeta = interpZeta(run,x,y,t);
			zeta = zeros(size(x));
		end
		
		function mask = interpMask(run,x,y,t);
			mask(isin) = griddata(grid.x,grid.y,grid.mask,x,y);
			% note that t is not used here
		end
		
		function u = interpU(run,x,y,sigma,t);
			u0 = griddata(grid.xu3,grid.yu3,grid.sigmau3,run.F0.u,x,y,sigma);
			u1 = griddata(grid.xu3,grid.yu3,grid.sigmau3,run.F1.u,x,y,sigma);
			u = run.tinterp(t, u0, u1);
		end

		function v = interpV(run,x,y,sigma,t);
			v0 = griddata(grid.xu3,grid.yu3,grid.sigmau3,run.F0.v,x,y,sigma);
			v1 = griddata(grid.xu3,grid.yu3,grid.sigmau3,run.F1.v,x,y,sigma);
			v = run.tinterp(t, v0, v1);
		end

		function w = interpW(run,x,y,sigma,t);
			w0 = griddata(grid.xw3,grid.yw3,grid.sigmaw3,run.F0.w,x,y,sigma);
			w1 = griddata(grid.xw3,grid.yw3,grid.sigmaw3,run.F1.w,x,y,sigma);
			w = run.tinterp(t, w0, w1);
		end
		
		function Ks = interpKs(run,x,y,sigma,t);
			Ks = interpTracer(run,'Ks',x,y,sigma,t);
		end
		
		function c = interpTracer(run,name,x,y,sigma,t);
			if ndims(run.F0.(name))==2 % 2d
				c0 = griddata(grid.x,grid.y,run.F0.(name),x,y);
				c1 = griddata(grid.x,grid.y,run.F1.(name),x,y);
			else % 3d
				c0 = griddata(grid.x3,grid.y3,grid.sigma3,run.F0.(name),...
							  x,y,sigma);
				c1 = griddata(grid.x3,grid.y3,grid.sigma3,run.F1.(name),...
							  x,y,sigma);
			end
			c = run.tinterp(t, c0, c1);
		end
		
		
		function us = scaleU(run,u,x,y); % m/s -> deg lon per day
			us = u .* 86400 ./ 111325 ./ cos(y./180.*pi);
		end
		function vs = scaleV(run,v,x,y); % m/s -> deg lat per day
			vs = v .* 86400 ./ 111325;
		end
		
		
		function [x1,y1,active] = filterCoordinates(run,x,y);
			active = y > grid.ymin;
			% confine x coordinates to 0...360
			x1 = mod(x,360);
			% deal with the case where a point goes over the pole (y>90).
			% here's a bad way of dealing with it:
			y1 = min(y1,90);
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




