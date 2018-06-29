classdef modelRun_cartesianSigma < modelRun

	% for working with output from models using a staggered,stretched Cartesian
	% grid in the horizontal and sigma coordinates in the vertical.
	%
	% romsCascadia and nemoAMM7 are built on top of this.
	%
	% returns variables with dimensions [lon lat depth], in that order.
	% accordingly, calls to interp2 for 2D variables look like
	% 	interp2(y,x,H,yi,xi)
	% and calls to interpn for 3D variables look like
	%	interpn(x,y,sigma,S,xi,yi,sigmai)
	% interpolation in time is done by interpolating in (x,y,sigma) at each of 
	% the currently two loaded frames and then interpolating in 1D in time 
	% between the results.

	properties
		F0, F1	 			  % two frames of in-memory storage
		outOfBoundsValue = 0; % set this to nan to have interp functions
							  % fail obviously, 0 to have them fail harmlessly
	end
	
	methods
	
		function run = modelRun_cartesianSigma();
		end
		
		
		function grid = add_3d_meshes(run,grid0);
			grid = grid0;
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
		end
		
		
		% reading from model files ---------------------------------------------
		
		
		% this, like the constructor which reads the grid, needs to be
		% adapted for each model.
		function run = loadFrame(run,n,tracers);
		end
		
		
		function run = advanceTo(run,n,tracers);
			run.F0 = run.F1;
			run.loadedN(1) = run.loadedN(2);
			run.loadFrame(n,tracers);
		end
		
		
		% interpolating model variables ----------------------------------------


		function c = interp(run,name,x,y,sigma,t);
			isin = run.in_xy_bounds(x,y);
			c = run.outOfBoundsValue .* ones(size(x));
			if ~isempty(isin)
				if strcmpi(name,'H')
					c(isin) = interp2(run.grid.lat, run.grid.lon, run.grid.H,... 
							      y(isin), x(isin));
				elseif ndims(run.F0.(name))==2 % 2D
					% sigma ignored
					% warning: ubar and vbar (ROMS) won't be handled correctly
					c0 = interp2(run.grid.lat, run.grid.lon, run.F0.(name), ...
								 y(isin), x(isin));
					c1 = interp2(run.grid.lat, run.grid.lon, run.F1.(name), ...
								 y(isin), x(isin));
					c(isin) = run.tinterp(t, c0, c1);				
				else % 3D
					if strcmpi(name,'u')
						gr = 'u3';
					elseif strcmpi(name,'v')
						gr = 'v3';
					elseif strcmpi(name,'w') | strcmpi(name,'Ks')
						gr = 'w3';
					else
						gr = 'rho3';
					end
					sigma1 = max(min(sigma,0),-1);
					c0 = interpn(run.grid.(gr).lon, run.grid.(gr).lat, ...			
						 run.grid.(gr).cs, ...
						 run.F0.(name), x(isin), y(isin), sigma1(isin));
					c1 = interpn(run.grid.(gr).lon, run.grid.(gr).lat, ...
						 run.grid.(gr).cs, ...
						 run.F1.(name), x(isin), y(isin), sigma1(isin));
					c(isin) = run.tinterp(t, c0, c1);
				end
			end
		end
		
		
		function c = interpDepthAverage(run,name,x,y,zMinMax,t);
			if ndims(run.F0.(name)) < 3
				c = run.interp(name,x,y,[],t);
			else % 3D
				cpro = run.interpProfile(name,x,y,t);
				vax = run.verticalAxisForProfiles;
				vax = vax(:)';
				c = depthAverage(cpro,vax,zMinMax);
			end
		end
		
		
		function v_axis = verticalAxisForProfiles(run);
			v_axis = run.grid.cs;
		end
		
		
		function c = interpProfile(run,name,x,y,t);
			v_axis = run.verticalAxisForProfiles;
			% if x,y are [NP 1], c is [NP length(v_axis)].
			if ndims(run.F0.(name)) < 3
				c = run.interp(name,x,y,[],t);
			else % 3D
				if strcmpi(name,'u')
					gr = 'u3';
				elseif strcmpi(name,'v')
					gr = 'v3';
				elseif strcmpi(name,'w') | strcmpi(name,'Ks')
					gr = 'w3';
				else
					gr = 'rho3';
				end
				xx = repmat(x(:),[1 length(v_axis)]);
				yy = repmat(y(:),[1 length(v_axis)]);
				vv = repmat(v_axis(:)',[size(xx,1) 1]);
				c0 = interpn(run.grid.(gr).lon, run.grid.(gr).lat, ...			
					 run.grid.(gr).cs, run.F0.(name), xx, yy, vv);
				c1 = interpn(run.grid.(gr).lon, run.grid.(gr).lat, ...
					 run.grid.(gr).cs, run.F1.(name), xx, yy, vv);
				c = run.tinterp(t, c0, c1);
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




