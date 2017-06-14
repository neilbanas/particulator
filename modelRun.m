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
		
		nativeSigma	% is the model on sigma as opposed to z levels?
					% used to speed things up if rel.sigmaTrapLevel or 
					% rel.zTrapLevel is set

		wScaleFactor % instead of methods like scaleU() and scaleV(), which
					% convert native u,v units to x,y units per day, assume
					% that the conversion factor is constant for
					% w, dKsdz, and wdiff
	end
	
	methods
		function run = modelRun;
		end
	
		function run = loadFrame(run,n,tracers);
		end
		function run = advanceTo(run,n,tracers);
		end
		
		function H = interpH(run,x,y);
		end
		function zeta = interpZeta(run,x,y,t);
		end
		function mask = interpMask(run,x,y,t);
		end
		
		function u = interpU(run,x,y,sigma,t);
		end
		function v = interpV(run,x,y,sigma,t);
		end
		function w = interpW(x,y,sigma,t);
		end
		function Ks = interpKs(x,y,sigma,t);
		end
		function c = interpTracer(run,name,x,y,sigma,t);
		end
		
		% these only need to be defined if nativeSigma = 0 
		function u = interpU_in_z(run,x,y,z,t);
		end
		function v = interpV_in_z(run,x,y,z,t);
		end
		function w = interpW_in_z(x,y,z,t);
		end
		function Ks = interpKs_in_z(x,y,z,t);
		end
		function c = interpTracer_in_z(run,name,x,y,z,t);
		end

		function us = scaleU(run,u,x,y);
		end
		function vs = scaleV(run,v,x,y);
		end
		
		function [x1,y1,active] = filterCoordinates(run,x,y);
		end		
	end

end