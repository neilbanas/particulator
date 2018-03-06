classdef modelRun_biomas2d < modelRun

	% for working with output from the Biomas model.
	% pre-averages over a specified depth range, so that all interpolations
	% are 2D.

	properties
		F0, F1					% two frames of in-memory storage
		year
		griddir = 'data/biomas/';
		dirname, basename
		localVars, fileVars		% lookup table for associating files with
								% standard variables
		tracerDims				% are the named variables 2D or 3D
		avg						% setup for depth-averaging
		pad						% setup for padding grid with strips at x=0,360
								% for wraparound interpolation
		si						% scatteredInterpolant objects for fields
								% that don't change (H,mask) 
	end
	
	methods
	
		function run = modelRun_biomas2d(dirname,year,griddir,depthRange);
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
			run.tracerDims = [3 3 3 3 2 2 2]; % 2D or 3D?
			
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
			grid.H(grid.mask) = -grid.zw(k(grid.mask));

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
			H3 = repmat(run.grid.H,[1 1 K]);
			run.avg.mask = zeros(I,J,K);
			run.avg.mask(:,:,run.avg.kRange(1):run.avg.kRange(2)) = 1;
			run.avg.mask(zw3 < -H3) = 0;
			run.avg.dz = dz3;
			run.avg.dz(run.avg.mask==0) = 0;
			run.avg.h = sum(run.avg.dz,3);
			
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
			
			% scatteredInterpolants
			run.si.H = scatteredInterpolant(...
				run.pad.x(:),run.pad.y(:),run.grid.H(run.pad.ind));
			run.si.mask = scatteredInterpolant(...
				run.pad.x(:),run.pad.y(:),double(run.grid.mask(run.pad.ind)));
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
					run.fileVars{strcmp('uv',run.localVars)} ...
					run.basename]);
			fseek(fid,framelength*(n-1),-1);
			A = reshape(fread(fid,I*J*K,'real*4'),[I J K]);
			run.F1.u = sum(A.*run.avg.dz,3) ./ run.avg.h;
			run.F1.u(~isfinite(run.F1.u)) = 0;
			fseek(fid,framelength*(365+n-1),-1);			
			A = reshape(fread(fid,I*J*K,'real*4'),[I J K]);
			run.F1.v = sum(A.*run.avg.dz,3) ./ run.avg.h;
			run.F1.v(~isfinite(run.F1.v)) = 0;
			fclose(fid);
				% perhaps could wrap these around the first dimension...
				% ([1:end 1],:)
				% ... and likewise for all other fields, in order to make
				% things interpolate correctly at x ~ 210 in the N Pacific.
			% Ks
			fid = fopen([run.dirname ...
					run.fileVars{strcmp('Ks',run.localVars)} ...
					run.basename]);
			fseek(fid,framelength*(n-1),-1);
			A = reshape(fread(fid,I*J*K,'real*4'),[I J K]);
			run.F1.Ks = sum(A.*run.avg.dz,3) ./ run.avg.h;
			run.F1.Ks(~isfinite(run.F1.Ks)) = 0;
			fclose(fid);
				% is Ks on the tracer grid or the w grid?
				% this version assumes tracer grid
			% everything else
			for i=1:length(tracers)
				j = find(strcmp(tracers{i},run.localVars));
				prefix = [run.fileVars{j}];
				if isempty(prefix), prefix = tracers{i}; end
				fid = fopen([run.dirname prefix run.basename]);
				if ~isempty(j) & run.tracerDims(j)==2 % 2D tracer
					fseek(fid,framelength/K*(n-1),-1);
					A = reshape(fread(fid,I*J,'real*4'),[I J]);
					run.F1.(tracers{i}) = A;
				else % 3D tracer
					fseek(fid,framelength*(n-1),-1);
					A = reshape(fread(fid,I*J*K,'real*4'),[I J K]);
					run.F1.(tracers{i}) = sum(A.*run.avg.dz,3) ./ run.avg.h;
					run.F1.(tracers{i})(~isfinite(run.F1.(tracers{i}))) = 0;
				end
				fclose(fid);
			end
			% create scatteredInterpolant objects for all fields
			run.F1.si.u = scatteredInterpolant(...
				run.pad.xu,run.pad.yu,run.F1.u(run.pad.indu));
			run.F1.si.v = scatteredInterpolant(...
				run.pad.xu,run.pad.yu,run.F1.v(run.pad.indu));
			run.F1.si.Ks = scatteredInterpolant(...
				run.pad.x,run.pad.y,run.F1.Ks(run.pad.ind));
			for i=1:length(tracers)
				run.F1.si.(tracers{i}) = scatteredInterpolant(...
					run.pad.x,run.pad.y,run.F1.(tracers{i})(run.pad.ind));
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
			u0 = run.F0.si.u(x,y);
			u1 = run.F1.si.u(x,y);
			u = run.tinterp(t, u0, u1);
		end

		function v = interpV(run,x,y,sigma,t);
			v0 = run.F0.si.v(x,y);
			v1 = run.F1.si.v(x,y);
			v = run.tinterp(t, v0, v1);
		end

		function w = interpW(run,x,y,sigma,t);
			w = 0;
		end
		
		function Ks = interpKs(run,x,y,sigma,t);
			Ks = interpTracer(run,'Ks',x,y,sigma,t);
		end
		
		function c = interpTracer(run,name,x,y,sigma,t);
			% note that sigma is ignored in all of these--2d interpolation only
			c0 = run.F0.si.(name)(x,y);
			c1 = run.F1.si.(name)(x,y);
			c = run.tinterp(t, c0, c1);
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




