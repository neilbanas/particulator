classdef modelRun_romsCascadia < modelRun

	% for working with output from MoSSea and Cascadia (PNWTOX) runs.
	% assumes a stretched Cartesian grid, and one output frame per file
	% (the _filestep_ property isn't used yet).
	%
	% unlike older code like 'post_tools' which was built around the snctools
	% for reading netcdfs, this code, using matlab's native netcdf library,
	% returns variables with dimensions [lon lat depth], in that order.
	% accordingly, calls to interp2 for 2D variables look like
	% 	interp2(y,x,H,yi,xi)
	% and calls to interpn for 3D variables look like
	%	interpn(x,y,sigma,S,xi,yi,sigmai)
	% interpolation in time is done by interpolating in (x,y,sigma) at each of 
	% the currently two loaded frames and then interpolating in 1D in time 
	% between the results.

	properties
		filename % list of files where output is stored 
		filestep % if more than one model time is stored per file,
				 % which index (for PNWTOX-era runs, this is just 1)
		F0, F1	 % two frames of in-memory storage
		outOfBoundsValue;
	end
	
	methods
	
		function run = modelRun_romsCascadia(dirname);
			run.outOfBoundsValue = 0; % set this to nan to have interp
									  % functions fail obviously, 0 to have
									  % them fail harmlessly
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
			grid.cs = [-1; grid.cs(:); 0];
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
			grid.Hu = interp2(grid.lat,grid.lon,grid.H,grid.lonu,grid.latu);
			grid.Hv = interp2(grid.lat,grid.lon,grid.H,grid.lonv,grid.latv);
			grid.bounds = [grid.lonu([1 end],1); grid.latv(1,[1 end])'];
			
			% construct 3d meshes for rho,u,v,w grids so this doesn't need to
			% be done for each interpolation
			Kw = length(grid.csw);
			K = length(grid.cs);
			[I,J] = size(grid.lon);
			[Iu,Ju] = size(grid.lonu);
			[Iv,Jv] = size(grid.lonv);
			grid.rho3.lon = repmat(reshape(grid.lon,[I J 1]),[1 1 K]);
			grid.rho3.lat = repmat(reshape(grid.lat,[I J 1]),[1 1 K]);
			grid.rho3.cs = repmat(reshape(grid.cs,[1 1 K]),[I J 1]);
			grid.u3.lon = repmat(reshape(grid.lonu,[Iu Ju 1]),[1 1 K]);
			grid.u3.lat = repmat(reshape(grid.latu,[Iu Ju 1]),[1 1 K]);
			grid.u3.cs = repmat(reshape(grid.cs,[1 1 K]),[Iu Ju 1]);
			grid.v3.lon = repmat(reshape(grid.lonv,[Iv Jv 1]),[1 1 K]);
			grid.v3.lat = repmat(reshape(grid.latv,[Iv Jv 1]),[1 1 K]);
			grid.v3.cs = repmat(reshape(grid.cs,[1 1 K]),[Iv Jv 1]);
			grid.w3.lon = repmat(reshape(grid.lon,[I J 1]),[1 1 Kw]);
			grid.w3.lat = repmat(reshape(grid.lat,[I J 1]),[1 1 Kw]);
			grid.w3.cs = repmat(reshape(grid.csw,[1 1 Kw]),[I J 1]);

			run.grid = grid;
		end % constructor
		
		
		% reading from model files ---------------------------------------------
		
		
		function run = loadFrame(run,n,tracers);
			% read from file
			% warning: does not use _filestep_: assumes one frame per file
			nc = netcdf.open(run.filename{n},'NOWRITE');
			run.F1.u = netcdf.getVar(nc,netcdf.inqVarID(nc,'u'),'double');
			run.F1.u = run.F1.u(:,:,[1 1:end end]);
			run.F1.v = netcdf.getVar(nc,netcdf.inqVarID(nc,'v'),'double');
			run.F1.v = run.F1.v(:,:,[1 1:end end]);
			run.F1.w = netcdf.getVar(nc,netcdf.inqVarID(nc,'w'),'double');
			run.F1.Ks = netcdf.getVar(nc,netcdf.inqVarID(nc,'AKs'),'double');
			run.F1.mask = ...
				netcdf.getVar(nc,netcdf.inqVarID(nc,'mask_rho'),'double');
			run.F1.zeta = netcdf.getVar(nc,netcdf.inqVarID(nc,'zeta'),'double');
			for i=1:length(tracers)
				run.F1.(tracers{i}) = ...
					netcdf.getVar(nc,netcdf.inqVarID(nc,tracers{i}),'double');
				run.F1.(tracers{i}) = run.F1.(tracers{i})(:,:,[1 1:end end]);
			end
			netcdf.close(nc);
			% replace netcdf out-of-range values (~1e37) with something more
			% useful. Velocity, diffusivity, and surface height should always 
			% be finite, so set these to 0. Tracers are set to nan, although it
			% might be better for this to be replaced with a nearest-neighbor
			% interpolation (otherwise, need to deal with this properly in
			% run.interpTracer()).
			run.F1.u(run.F1.u > 1e36) = 0;
			run.F1.v(run.F1.v > 1e36) = 0;
			run.F1.w(run.F1.w > 1e36) = 0;
			run.F1.Ks(run.F1.Ks > 1e36) = 0;
			run.F1.zeta(run.F1.zeta > 1e36) = 0;
			for i=1:length(tracers)
				run.F1.(tracers{i})(run.F1.(tracers{i}) > 1e36) = nan;
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
			isin = run.in_xy_bounds(x,y);
			H = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				H(isin) = interp2(run.grid.lat, run.grid.lon, run.grid.H, ...
							      y(isin), x(isin));
			end
		end
		
		function zeta = interpZeta(run,x,y,t);
			isin = run.in_xy_bounds(x,y);
			zeta = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				zeta0 = interp2(run.grid.lat, run.grid.lon, run.F0.zeta, ...
								y(isin), x(isin));
				zeta1 = interp2(run.grid.lat, run.grid.lon, run.F1.zeta, ...
								y(isin), x(isin));
				zeta(isin) = run.tinterp(t, zeta0, zeta1);
			end
		end
		
		function mask = interpMask(run,x,y,t);
			isin = run.in_xy_bounds(x,y);
			mask = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				mask0 = interp2(run.grid.lat, run.grid.lon, run.F0.mask, ...
								y(isin), x(isin));
				mask1 = interp2(run.grid.lat, run.grid.lon, run.F1.mask, ...
								y(isin), x(isin));
				mask(isin) = run.tinterp(t, mask0, mask1);
			end
		end
		
		function u = interpU(run,x,y,sigma,t);
			isin = run.in_xy_bounds(x,y);
			sigma1 = max(min(sigma,0),-1);
			u = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				u0 = interpn(...
				     run.grid.u3.lon, run.grid.u3.lat, run.grid.u3.cs, ...
					 run.F0.u, x(isin), y(isin), sigma1(isin));
				u1 = interpn(...
				     run.grid.u3.lon, run.grid.u3.lat, run.grid.u3.cs, ...
					 run.F1.u, x(isin), y(isin), sigma1(isin));
				u(isin) = run.tinterp(t, u0, u1);
			end
		end

		function v = interpV(run,x,y,sigma,t);
			isin = run.in_xy_bounds(x,y);
			sigma1 = max(min(sigma,0),-1);
			v = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				v0 = interpn(...
				     run.grid.v3.lon, run.grid.v3.lat, run.grid.v3.cs, ...
					 run.F0.v, x(isin), y(isin), sigma1(isin));
				v1 = interpn(...
				     run.grid.v3.lon, run.grid.v3.lat, run.grid.v3.cs, ...
					 run.F1.v, x(isin), y(isin), sigma1(isin));
				v(isin) = run.tinterp(t, v0, v1);
			end
		end
		
		function w = interpW(run,x,y,sigma,t);
			isin = run.in_xy_bounds(x,y);
			sigma1 = max(min(sigma,0),-1);
			w = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				w0 = interpn(...
				     run.grid.w3.lon, run.grid.w3.lat, run.grid.w3.cs, ...
					 run.F0.w, x(isin), y(isin), sigma1(isin));
				w1 = interpn(...
				     run.grid.w3.lon, run.grid.w3.lat, run.grid.w3.cs, ...
					 run.F1.w, x(isin), y(isin), sigma1(isin));
				w(isin) = run.tinterp(t, w0, w1);
			end
		end
		
		function Ks = interpKs(run,x,y,sigma,t);
			isin = run.in_xy_bounds(x,y);
			sigma1 = max(min(sigma,0),-1);
			Ks = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				Ks0 = interpn(...
				     run.grid.w3.lon, run.grid.w3.lat, run.grid.w3.cs, ...
					 run.F0.Ks, x(isin), y(isin), sigma1(isin));
				Ks1 = interpn(...
				     run.grid.w3.lon, run.grid.w3.lat, run.grid.w3.cs, ...
					 run.F1.Ks, x(isin), y(isin), sigma1(isin));
				Ks(isin) = run.tinterp(t, Ks0, Ks1);
			end
		end
		
		function c = interpTracer(run,name,x,y,sigma,t);
			% warning: either here or in run.loadFrame(), have to deal with the
			% case where we're interpolating between the last wet cell and the
			% land, since land values are allowed to be nan
			isin = run.in_xy_bounds(x,y);
			sigma1 = max(min(sigma,0),-1);
			c = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				if ndims(run.F0.(name))==2 % 2d
					c0 = interp2(run.grid.lat, run.grid.lon, run.F0.(name), ...
								 x(isin), y(isin));
					c1 = interp2(run.grid.lat, run.grid.lon, run.F1.(name), ...
								 x(isin), y(isin));
					c(isin) = run.tinterp(t, c0, c1);				
				else % 3d
					c0 = interpn(run.grid.rho3.lon, run.grid.rho3.lat, ...			
						 run.grid.rho3.cs, ...
						 run.F0.(name), x(isin), y(isin), sigma1(isin));
					c1 = interpn(run.grid.rho3.lon, run.grid.rho3.lat, ...
						 run.grid.rho3.cs, ...
						 run.F1.(name), x(isin), y(isin), sigma1(isin));
					c(isin) = run.tinterp(t, c0, c1);
				end
			end
		end
		
		
		function us = scaleU(run,u,x,y); % m/s -> deg lon per day
			us = u .* 86400 ./ 111325 ./ cos(y./180.*pi);
		end
		function vs = scaleV(run,v,x,y); % m/s -> deg lat per day
			vs = v .* 86400 ./ 111325;
		end
		function ws = scaleW(run,w);
			ws = w .* 86400; % m/s -> m/day
		end
		
		
		function isin = in_xy_bounds(run,x,y);
			isin = x >= run.grid.bounds(1) & x <= run.grid.bounds(2) & ...
				   y >= run.grid.bounds(3) & y <= run.grid.bounds(4);
		end		
		function [x1,y1,active] = filterCoordinates(run,x,y);
			active = run.in_xy_bounds(x,y);
			x1 = x;
			y1 = y;
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




