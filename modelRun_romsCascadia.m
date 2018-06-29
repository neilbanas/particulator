classdef modelRun_romsCascadia < modelRun_cartesianSigma

	% for working with output from MoSSea and Cascadia (PNWTOX) runs.
	%
	% most of the code that used to be here has been moved to
	% modelRun_cartesianSigma.
	

	properties
		filename % list of files where output is stored 
		filestep % if more than one model time is stored per file,
				 % which index (for PNWTOX-era runs, this is just 1)
	end
	
	methods
	
		function run = modelRun_romsCascadia(dirname);
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
			
			grid = run.add_3d_meshes(grid);
			
			% find nearest wet cell for every dry cell, on the rho grid.
			% this makes it much faster to fill in land-masked tracer values
			% so that interpolation works correctly near the land edge.
			grid.nearestWater = reshape(1:length(grid.lon(:)), size(grid.lon));
			dry = grid.mask < 1;
			wet = ~dry;
			grid.nearestWater(dry) = griddata(...
				grid.lon(wet), grid.lat(wet), grid.nearestWater(wet), ...
				grid.lon(dry), grid.lat(dry), 'nearest');

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
			% be finite, so set these to 0.
			run.F1.u(run.F1.u > 1e36) = 0;
			run.F1.v(run.F1.v > 1e36) = 0;
			run.F1.w(run.F1.w > 1e36) = 0;
			run.F1.Ks(run.F1.Ks > 1e36) = 0;
			run.F1.zeta(run.F1.zeta > 1e36) = 0;
			for i=1:length(tracers)
				% bad tracer values could simply be set to nan...
%				run.F1.(tracers{i})(run.F1.(tracers{i}) > 1e36) = nan;
				% ... however, this causes interpolation problems when particles
				% are within one grid cell of the land. Filling in land values
				% with a nearest-neighbour interpolation once here means that
				% we don't have to do it repeatedly in run.interp().
				for k=1:size(run.F1.(tracers{i}),3)
					Ck = run.F1.(tracers{i})(:,:,k);
					dry = Ck > 1e36;
					Ck(dry) = Ck(run.grid.nearestWater(dry));
					run.F1.(tracers{i})(:,:,k) = Ck;
				end
			end
			% declare this frame loaded
			run.loadedN(2) = n;
		end
		

	end % methods

end % classdef




