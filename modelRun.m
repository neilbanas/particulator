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
	end
	
	methods
		function run = modelRun;
		end
	
		function run = loadFrame(run,n,tracers);
		end
		function run = advance(run,tracers);
		end
		
		function H = interpH(run,y,x);
		end
		function zeta = interpZeta(run,t,y,x);
		end
		function mask = interpMask(run,t,y,x);
		end
		
		function u = interpU(run,t,sigma,y,x);
		end
		function v = interpV(run,t,sigma,y,x);
		end
		function w = interpW(t,sigma,y,x);
		end
		function Ks = interpKs(t,sigma,y,x);
		end
		function c = interpTracer(run,name,t,sigma,y,x);
		end
		
		function us = scaleU(run,u,y,x);
		end
		function vs = scaleV(run,v,y,x);
		end
		
		function isin = in_xy_bounds(run,y,x);
		end		
	end

end