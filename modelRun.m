% base class definition for model runs.
% roms_cascadia, etc. should extend this class.

classdef modelRun < handle

	properties
		numFrames	% saved model output is identified by a series of indices 			
					% (1:numFrames)
					% even if this is arbitrary (as in an analytical flow 
					% field, set it so that a large fraction of one frame 
					% interval is a natural timestep to integrate over
					
		t			% times corresponding to (1:numFrames)
		
		loadedN		% t(loadedN) is the time range ([start end]) that is
					% saved in memory to be integrated within. We assume that
					% loadedN(2) - loadedN(1) = 1, so that
					% diff(t(loadedN)) * release.dt_per_DT is the timestep
					% for particle integration
					
		grid		% model grid (details completely model-dependent)
	end
	
	methods
		function run = modelRun;
		end
	
		function run = loadFrame(run,n,tracers);
		end
		function run = advanceTo(run,n,tracers);
		end
		
		function c = interp(run,name,x,y,sigma,t);
			% if name == 'H', callable as interp('H',x,y)
			% if a 2D variable, interp(name,x,y,[],t)				
		end
		function c = interpDepthAverage(run,name,x,y,zMinMax,t);
		end
		function c = interpProfile(run,name,x,y,t);
		end
		function v_axis = verticalAxisForProfiles(run);
		end

		function us = scaleU(run,u,x,y); % native units -> deg lon/day
		end
		function vs = scaleV(run,v,x,y); % native units -> deg lat/day
		end
		function ws = scaleW(run,w); % native units -> m/day
		end
		
		function [x1,y1,active] = filterCoordinates(run,x,y);
		end		
	end

end