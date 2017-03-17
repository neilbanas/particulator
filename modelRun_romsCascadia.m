classdef modelRun_romsCascadia < modelRun

	properties
		filename % list of files where output is stored 
		filestep % if more than one model time is stored per file,
				 % which index (for PNWTOX-era runs, this is just 1)
		F0, F1	 % two frames of in-memory storage
		outOfBoundsValue;
	end
	
	methods
	
		function run = modelRun_romsCascadia(dirname);
			run.outOfBoundsValue = 0; % set this to nan to have interp functions
								  	  % fail obviously, 0 to have them fail
								  	  % harmlessly
			
			% locate file list
			thefiles = dir([dirname 'ocean_his_*.nc']);
			run.filename = {thefiles.name};
			if isempty(run.filename)
				warning(['no files found at ' dirname 'ocean_his_*.nc']);
			end
			run.numFrames = length(run.filename);
			run.filestep = ones(run.numFrames,1);
			for i=1:run.numFrames
				run.filename{i} = [dirname run.filename{i}];
			end
			
			% look in first and last files to get times
			nc = netcdf.open(run.filename{1},'NOWRITE');
			t0 = netcdf.getVar(nc,netcdf.inqVarID(nc,'ocean_time'),'double');
			units = netcdf.getAtt(nc,netcdf.inqVarID(nc,'ocean_time'),'units');
			netcdf.close(nc);
			nc = netcdf.open(run.filename{end},'NOWRITE');
			t1 = netcdf.getVar(nc,netcdf.inqVarID(nc,'ocean_time'),'double');
			netcdf.close(nc);
			timeref = strrep(units, 'seconds since ', '');
			timeref = datenum(timeref,'yyyy-mm-dd HH:MM:SS');
			run.t = linspace(t0,t1,run.numFrames) ./ 86400 + timeref;
			run.loadedN = [nan nan];
			
			% load grid from first file
			nc = netcdf.open(run.filename{1},'NOWRITE');
			grid.lon = netcdf.getVar(nc,netcdf.inqVarID(nc,'lon_rho'),'double');
			grid.lat = netcdf.getVar(nc,netcdf.inqVarID(nc,'lat_rho'),'double');
			grid.lonu = netcdf.getVar(nc,netcdf.inqVarID(nc,'lon_u'),'double');
			grid.latu = netcdf.getVar(nc,netcdf.inqVarID(nc,'lat_u'),'double');
			grid.lonv = netcdf.getVar(nc,netcdf.inqVarID(nc,'lon_v'),'double');
			grid.latv = netcdf.getVar(nc,netcdf.inqVarID(nc,'lat_v'),'double');
			grid.DX = 1 ./ netcdf.getVar(nc,netcdf.inqVarID(nc,'pm'),'double');
			grid.DY = 1 ./ netcdf.getVar(nc,netcdf.inqVarID(nc,'pn'),'double');
			grid.cs = netcdf.getVar(nc,netcdf.inqVarID(nc,'Cs_r'),'double');
			grid.csw = netcdf.getVar(nc,netcdf.inqVarID(nc,'Cs_w'),'double');
			grid.mask = ...
					netcdf.getVar(nc,netcdf.inqVarID(nc,'mask_rho'),'double');
			grid.masku = ...
					netcdf.getVar(nc,netcdf.inqVarID(nc,'mask_u'),'double');
			grid.maskv = ...
					netcdf.getVar(nc,netcdf.inqVarID(nc,'mask_v'),'double');
			grid.H = ...
					netcdf.getVar(nc,netcdf.inqVarID(nc,'h'),'double');
			netcdf.close(nc);
			grid.Hu = interp2(grid.lat,grid.lon,grid.H,grid.latu,grid.lonu);
			grid.Hv = interp2(grid.lat,grid.lon,grid.H,grid.latv,grid.lonv);
			grid.bounds = [grid.lonu(1,[1 end])' grid.latv([1 end],1)];
			
			% construct 3d meshes for rho,u,v,w grids so this doesn't need to
			% be done for each interpolation
			grid.rho3.cs = ...
			
			run.grid = grid;
		end % constructor
		
		
		% reading from model files ---------------------------------------------
		
		
		function run = loadFrame(run,n,tracers);
			% !!! does not use filestep: assumes one frame per file !!!
			nc = netcdf.open(run.filename{n},'NOWRITE');
			run.F1.u = netcdf.getVar(nc,netcdf.inqVarID(nc,'u'),'double');
			run.F1.v = netcdf.getVar(nc,netcdf.inqVarID(nc,'v'),'double');
			run.F1.w = netcdf.getVar(nc,netcdf.inqVarID(nc,'w'),'double');
			run.F1.Ks = netcdf.getVar(nc,netcdf.inqVarID(nc,'AKs'),'double');
			run.F1.mask = ...
				netcdf.getVar(nc,netcdf.inqVarID(nc,'mask_rho'),'double');
			run.F1.zeta = netcdf.getVar(nc,netcdf.inqVarID(nc,'zeta'),'double');
			for i=1:length(tracers)
				run.F1.(tracers{i}) = ...
					netcdf.getVar(nc,netcdf.inqVarID(nc,tracers{i}),'double');
			end
			netcdf.close(nc);
			run.loadedN(2) = n;
		end
		
		function run = advanceTo(run,n,tracers);
			run.F0 = run.F1;
			run.loadFrame(n,tracers);
			run.loadedN = [run.loadedN(2) n];
		end
		
		
		% interpolating model variables ----------------------------------------


		function H = interpH(run,y,x);
			isin = run.in_xy_bounds(y,x);
			H = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				H(isin) = interp2(run.grid.lat, run.grid.lon, run.grid.H, ...
							      y(isin), x(isin));
			end
		end
		
		function zeta = interpZeta(run,t,y,x);
			isin = run.in_xy_bounds(y,x);
			zeta = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				zeta0 = interp2(run.grid.lat, run.grid.lon, run.F0.zeta, ...
								y(isin), x(isin));
				zeta1 = interp2(run.grid.lat, run.grid.lon, run.F1.zeta, ...
								y(isin), x(isin));
				zeta(isin) = run.tinterp(t, zeta0, zeta1);
			end
		end
		
		function mask = interpMask(run,t,y,x);
			isin = run.in_xy_bounds(y,x);
			mask = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				mask0 = interp2(run.grid.lat, run.grid.lon, run.F0.mask, ...
								y(isin), x(isin));
				mask1 = interp2(run.grid.lat, run.grid.lon, run.F1.mask, ...
								y(isin), x(isin));
				mask(isin) = run.tinterp(t, mask0, mask1);
			end
		end
		
		function u = interpU(run,t,sigma,y,x);
			isin = run.in_xy_bounds(y,x);
			sigma1 = max(min(sigma,0),-1);
			u = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				u0 = interpn(run.grid.cs, run.grid.latu, run.grid.lonu, ...
							 run.F0.u, sigma1(isin), y(isin), x(isin));
				u1 = interpn(run.grid.cs, run.grid.latu, run.grid.lonu, ...
							 run.F1.u, sigma1(isin), y(isin), x(isin));
				u(isin) = run.tinterp(t, u0, u1);
			end
		end

		function v = interpV(run,t,sigma,y,x);
			isin = run.in_xy_bounds(y,x);
			sigma1 = max(min(sigma,0),-1);
			v = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				v0 = interpn(run.grid.cs, run.grid.latv, run.grid.lonv, ...
							 run.F0.v, sigma1(isin), y(isin), x(isin));
				v1 = interpn(run.grid.cs, run.grid.latv, run.grid.lonv, ...
							 run.F1.v, sigma1(isin), y(isin), x(isin));
				v(isin) = run.tinterp(t, v0, v1);
			end
		end
		
		function w = interpW(t,sigma,y,x);
			isin = run.in_xy_bounds(y,x);
			sigma1 = max(min(sigma,0),-1);
			w = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				w0 = interpn(run.grid.csw, run.grid.lat, run.grid.lon, ...
							 run.F0.w, sigma1(isin), y(isin), x(isin));
				w1 = interpn(run.grid.csw, run.grid.lat, run.grid.lon, ...
							 run.F1.w, sigma1(isin), y(isin), x(isin));
				w = run.tinterp(t, w0, w1);
			end
		end
		
		function Ks = interpKs(t,sigma,y,x);
			isin = run.in_xy_bounds(y,x);
			sigma1 = max(min(sigma,0),-1);
			Ks = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				Ks0 = interpn(run.grid.csw, run.grid.lat, run.grid.lon, ...
							  run.F0.Ks, sigma1(isin), y(isin), x(isin));
				Ks1 = interpn(run.grid.csw, run.grid.lat, run.grid.lon, ...
							  run.F1.Ks, sigma1(isin), y(isin), x(isin));
				Ks(isin) = run.tinterp(t, Ks0, Ks1);
			end
		end
		
		function c = interpTracer(run,name,t,sigma,y,x);
			isin = run.in_xy_bounds(y,x);
			sigma1 = max(min(sigma,0),-1);
			c = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				if ndims(run.F0(name))==2 % 2d
					c0 = interp2(run.grid.lat, run.grid.lon, run.F0.(name), ...
								 y(isin), x(isin));
					c1 = interp2(run.grid.lat, run.grid.lon, run.F1.(name), ...
								 y(isin), x(isin));
					c(isin) = run.tinterp(t, c0, c1);				
				else % 3d
					c0 = interpn(run.grid.cs, run.grid.lat, run.grid.lon, ...
								 run.F0.(name), sigma1(isin), y(isin), x(isin));
					c1 = interpn(run.grid.cs, run.grid.lat, run.grid.lon, ...
								 run.F1.(name), sigma1(isin), y(isin), x(isin));
					c(isin) = run.tinterp(t, c0, c1);
				end
			end
		end
		
		
		function us = scaleU(run,u,y,x); % m/s -> deg lon per day
			us = u .* 86400 ./ 111325 ./ cos(y./180.*pi);
		end
		function vs = scaleV(run,v,y,x); % m/s -> deg lat per day
			vs = v .* 86400 ./ 111325;
		end
		
		
		function isin = in_xy_bounds(run,y,x);
			isin = x >= run.grid.bounds(1) & x <= run.grid.bounds(2) & ...
				   y <= run.grid.bounds(3) & y <= run.grid.bounds(4);
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




