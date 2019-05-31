classdef modelRun_nemoGlobal_zlevel < modelRun_nemoGlobal2d

	% for working with output from NOC Southampton NEMO/MEDUSA global model.
	% so far works on test files from Andy Yool, Nov 2017, from a 1/12 deg
	% hindcast
	%
	% picks a single z level near the mean of depthRange: never reads 3d fields into
	% memory. This is a solution to the problem that nemoGlobal2d was using > 50 GB of
	% swap space to load the pan-arctic region.	

	properties
	
		zIndex
		
	end
	
	
	methods
	
		function depthAveragingSetup(run,depthRange,dz3,mask3);
			run.avg = [];
			run.zIndex = find(abs(run.grid.z - mean(depthRange)) == ...
				min(abs(run.grid.z - mean(depthRange))));
			disp(['reading from grid.z(' num2str(run.zIndex) ') = ' ...
				num2str(run.grid.z(run.zIndex))]);
		end	
		
		
		
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
					c = netcdf.getVar(nc,varid, [run.box.start run.zIndex-1 0], ...
									[run.box.count 1 1], 'double');
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
			% don't load Ks by default: you can always include it in the list of tracers
			run.F1.Ks = zeros(size(run.F1.u));
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
					C = run.F1.(fields{i});
					warning off
					run.F1.si.(fields{i}) = scatteredInterpolant( ...
						run.grid.x(:), run.grid.y(:), C(:));
					warning on
				end
			end	
		end
	
		

		function v_axis = verticalAxisForProfiles(run);
			v_axis = [];
		end
		
		function c = interpProfile(run,name,x,y,t);
			c = [];
		end
		
		
		
	end % methods

end % classdef




