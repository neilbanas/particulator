classdef modelRun_nemoAMM7 < modelRun_cartesianSigma

	% for working with output from the AMM7 model (NOC/PML).
	%
	% AMM7 uses a Cartesian grid with stretched vertical coordinates, and so
	% interpolation follows the same methods as romsCascadia, and is quite
	% different from modelRun_nemoGlobal, which uses a curvilinear grid and 
	% z coordinates.
	%
	% one difference from ROMS is that the direction of the vertical axis is
	% reversed; I haven't thought through whether this makes a difference.
	

	properties
		filenameTemplate
		
		vars					% lookup table for associating filenames with
								% standard variable names
	end
	
	
	methods
	
		function run = modelRun_nemoAMM7(basedir,basename);
			run.filenameTemplate = [basedir basename '_grid_[var].nc']
			% basename is like 'amm7_1d_19810101_19810131'

			%internal name, nc file, name inside nc file
			tab = {...
				'u',	'U',	'uoce'; ...		% m/s
				'v',	'V',	'voce'};
			run.vars.local = tab(:,1);
			run.vars.ncfile = tab(:,2);
			run.vars.ncname = tab(:,3);

			% read tracer grid + timebase
			% ...			
			
			% read u grid and timebase	
			ncname = strrep(run.filenameTemplate,'[var]','U');
			nc = netcdf.open(ncname,'NOWRITE');
			lon1 = netcdf.getVar(nc,netcdf.inqVarID(nc,'x_grid_U'),'double');			
			lat1 = netcdf.getVar(nc,netcdf.inqVarID(nc,'y_grid_U'),'double');
			[grid.latu,grid.lonu] = meshgrid(lat1,lon1);
			t = netcdf.getVar(nc,netcdf.inqVarID(nc,'time_counter'),'double');
			run.t = t./86400 + datenum(1950,1,1);
			run.numFrames = length(run.t);
			netcdf.close(nc);
			% read v grid	
			ncname = strrep(run.filenameTemplate,'[var]','V');
			nc = netcdf.open(ncname,'NOWRITE');
			lon1 = netcdf.getVar(nc,netcdf.inqVarID(nc,'x_grid_V'),'double');			
			lat1 = netcdf.getVar(nc,netcdf.inqVarID(nc,'y_grid_V'),'double');
			[grid.latv,grid.lonv] = meshgrid(lat1,lon1);
			netcdf.close(nc);
			% read tracer grid
			ncname = strrep(run.filenameTemplate,'[var]','T');
			nc = netcdf.open(ncname,'NOWRITE');
			lon1 = netcdf.getVar(nc,netcdf.inqVarID(nc,'x_grid_T'),'double');			
			lat1 = netcdf.getVar(nc,netcdf.inqVarID(nc,'y_grid_T'),'double');
			[grid.lat,grid.lon] = meshgrid(lat1,lon1);
			netcdf.close(nc);
			% read layer depths
			nc = netcdf.open([basedir 'mesh_zgr_e3.nc'],'NOWRITE');
			e3u = netcdf.getVar(nc,netcdf.inqVarID(nc,'e3u_0'),'double');
			e3v = netcdf.getVar(nc,netcdf.inqVarID(nc,'e3v_0'),'double');
			netcdf.close(nc);
			grid.Hu = sum(e3u,3); % total water depth
			grid.Hv = sum(e3v,3);
			grid.H = interp2(grid.latu,grid.lonu,grid.Hu,grid.lat,grid.lon);
			dcs = squeeze(e3u(2,2,:)./grid.Hu(2,2));
			grid.csw = [0; -cumsum(dcs)];
				% this is what would be, in ROMS, the vertical coordinates
				% for w, running all the way from -1 to the free surface.
				% w in AMM7 is not quite the same size.
				% if this is changed so that w can be read in for AMM7,
				% look through modelRun_cartesianSigma for places where we
				% are assuming that csw includes -1 and 0.
			grid.cs = 0.5.*(grid.csw(1:end-1)+grid.csw(2:end));
				 % scaled vertical coordinates (sigma coordinates, -1..0)
			grid.cs = [0; grid.cs(:); -1];
				% when we read in model fields, we'll fill in values at the
				% surface and bottom so that we can interpolate to all valid 
				% sigma values


			grid.bounds = [grid.lonu([1 end],1); grid.latv(1,[1 end])'];
			
			% don't have a real land mask at this point, so pretend it's
			% water everywhere
			grid.mask = ones(size(grid.lonu));
			grid.nearestWater = ...
				reshape(1:length(grid.mask(:)),size(grid.mask));
							
			grid = run.add_3d_meshes(grid);		
			run.grid = grid;				
		end % constructor
		
		
		% reading from model files ---------------------------------------------
		
		
		function c = read(run,localVarname,n);
			vi = strmatch(localVarname,run.vars.local,'exact');
			if ~isempty(vi) 
				ncfile = run.vars.ncfile{vi};
				ncname = strrep(run.filenameTemplate,'[var]',ncfile);
				nc = netcdf.open(ncname,'NOWRITE');
				disp(['reading ' run.vars.ncname{vi} ' from ' ncname]);
				varid = netcdf.inqVarID(nc,run.vars.ncname{vi});
				[~,~,dimids,~] = netcdf.inqVar(nc,varid);
				[I,J] = size(run.grid.lon);
				if length(dimids)==4 % 3D x time
					K = length(run.grid.cs) - 2;
						% not correct for Ks or w
					c = netcdf.getVar(nc,varid, [0 0 0 n-1], ...
									[I J K 1], 'double');
					c = c(:,:,[1 1:end end]);
						% extend values to the surface and bottom
				else % 2D x time
					c = netcdf.getVar(nc,varid, [0 0 n-1], ...
									[I J 1], 'double');
				end
			else
				warning(['don''t know what file contains ' localVarname '.']);
				c = [];
			end
			netcdf.close(nc);
			c(c>1e16) = nan;
		end
		
		
		function run = loadFrame(run,n,tracers);
			if nargin < 2, tracers = {}; end
			% load U, V
			run.F1.u = run.read('u',n);
			run.F1.v = run.read('v',n);
			% fill in zeros for Ks and w
			run.F1.Ks = zeros(size(run.grid.w3.lon));
			run.F1.w = zeros(size(run.grid.w3.lon));
			% same mask for all timesteps
			run.F1.mask = run.grid.mask;
			% no free surface variation
			run.F1.zeta = zeros(size(run.grid.H));
			% load tracers
			for m=1:length(tracers)
				run.F1.(tracers{m}) = run.read(tracers{m},n);
			end
			run.loadedN(2) = n;
		end
	
		
	end % methods

end % classdef




