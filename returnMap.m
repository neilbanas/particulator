classdef returnMap < handle

	properties
		x,y,H	% grid
		t		% timebase
		
		M		% the map itself (lookup table of indices into x,y).
				% a particle that is at x(i),y(i) at time t(n) will be at
				% x(M(n,i,r)),y(M(n,i,r)) at time n+1, where r is a random
				% index that chooses among replicates
		c		% tracers and other scalar properties at (x,y,t) points
		
		xadv, yadv		% raw results of advection over each frame: allows one
						% to change choices about diffusion without having to
						% redo the whole integration

		maskingGrid		% the grid from the source model or some other means
						% of determining whether a point is valid when xy2ind
						% is called with 'strict' 
						% (fields = .x, .y, .mask)
		si, siMask	% scatteredInterpolant objects to speed up xy2ind
	end
	
	
	methods
	
		function map = returnMap(x,y,t);
			% a blank structure
			map.x = x;
			map.y = y;
			map.si = scatteredInterpolant(x(:),y(:),(1:length(x(:)))',...
						'nearest');
			map.t = t;
		end
		
		
		% map generation -------------------------------------------------------
		% (a big particle integration using par_integrate)
		
		function map = generate(map,run,varargin);
			% varargin is passed along to the particle release that is
			% performed for each time interval in map.t
			% it also contains values for 'KH' and 'Nreplicates'
			% 'tracers' is used twice, by both generate() and par_release()			
			% defaults
			opt.KH = 0;
			opt.Nreplicates = 1;
			opt.tracers = {};
			% overwrite defaults
			fields = {varargin{1:2:end}};
			vals = {varargin{2:2:end}};
			for i=1:length(fields)
				opt.(fields{i}) = vals{i};
			end
			% set the mask used by xy2ind(...'strict')
			map.setMask(run);
			% initialize the map M
			map.M = repmat(nan,...
						[length(map.t)-1 prod(size(map.x)) opt.Nreplicates]);
			% interpolate tracers for the first timestep
			sigma0 = nan; % haven't figured out how this should behave--
						  % so for now only covering cases where this is
						  % ignored anyway, like biomas2d
			map.H = run.interpH(map.x,map.y);
			n00 = floor(interp1(run.t, 1:run.numFrames, map.t(1)));
			run.loadFrame(n00,opt.tracers);
			run.advanceTo(n00,opt.tracers);
			for j=1:length(opt.tracers)
				map.c.(opt.tracers{j})(1,:) = run.interpTracer(...
					opt.tracers{j},map.x,map.y,sigma0,map.t(1));
			end
			% advance in time...
			for ni = 1:(length(map.t)-1)
				% perform a normal particle release in the modelRun _run_
				% for one interval in map.t
				rel = par_release('x0',map.x,'y0',map.y,'t0',map.t(ni),...
					 			  't1',map.t(ni+1),'sigma0',sigma0);
				fields = fieldnames(opt);
				for j=1:length(fields)
					rel.(fields{j}) = opt.(fields{j}); % pass along options
				end
				steps = par_integrate(rel,run);
				map.xadv(ni,:) = steps(end).x(:)';
				map.yadv(ni,:) = steps(end).y(:)';
				% the interpolation of tracers from _run_ onto map.x,map.y,map.t
				% is actually independent of the map generation. But while we 
				% have the timestep corresponding to map.t(ni) loaded in _run_
				% (from the end of par_integrate above) it's a good time to
				% look up tracer values...
				for j=1:length(opt.tracers)
					map.c.(opt.tracers{j})(ni+1,:) = run.interpTracer(...
						opt.tracers{j},map.x,map.y,sigma0,map.t(ni+1));
				end
			end
			% if we have no diffusion and no replicates, then
			% M = xy2ind(xadv,yadv) and we'd be done. But instead, add diffusion 
			% (this can be overwitten with a different diffusivity later)		
			map.addDiffusion(opt.KH,opt.Nreplicates);
		end
		
		
		function map = addDiffusion(map,KH,Nreplicates,doResampling);
			if nargin < 4, doResampling = 1; end
			if nargin < 3, Nreplicates = 1; end
			if nargin < 2, KH = 0; end
			if KH==0, Nreplicates = 1; end			
			% add a random normal distribution to xadv, yadv, scaled by the
			% diffusivity KH. Assumes KH is in m^2/s, x,y in lon,lat
			x = repmat(map.xadv,[1 1 Nreplicates]);
			y = repmat(map.yadv,[1 1 Nreplicates]);
			dt = map.t(2) - map.t(1); % assume it's a uniform timestep
			dt_secs = dt * 86400;
			udiff = sqrt(2*KH/dt_secs) .* randn(size(x)) .* ...
					86400 ./ 111325 ./ cos(y./180.*pi);
			vdiff = sqrt(2*KH/dt_secs) .* randn(size(x)) .* ...
					86400 ./ 111325;
			x = x + udiff .* dt;
			y = y + vdiff .* dt;			
			% figure out where in (map.x, map.y) these have ended up
			map.M = map.xy2ind(x,y,'strict');
			if doResampling
				% resample to replace the ones that ended at invalid locations
				[nbad,ibad,rbad] = ind2sub(size(map.M),find(~isfinite(map.M)));
				nbadu = unique(nbad);
				for ni=1:length(nbadu)
					n = nbadu(ni);
					ibadu = unique(ibad(nbad==n));
					for ii=1:length(ibadu)
						i = ibadu(ii);
						% each (n,i) that had bad values among its replicates
						rr = rbad(nbad==n & ibad==i);
						good = find(isfinite(map.M(n,i,:)));
						if length(good)>0
							j = randi(length(good),[length(rr) 1]);
							map.M(n,i,rr) = map.M(n,i,good(j));
						else
							% if no good values, make it a fixed pt
							map.M(n,i,rr) = i;
						end
					end
				end
			end
		end
		
		
		% coordinate conversions -----------------------------------------------
		
		function ind = xy2ind(map,x,y,howStrict);
			if nargin < 4, howStrict = 'not strict'; end
			if isempty(map.maskingGrid), howStrict = 'not strict'; end
			ind = map.si(x,y);
			if strcmpi(howStrict,'strict')
				m = map.siMask(x,y);
				ind(m < 0.5) = nan;
			end
		end
		
		function map = setMask(map,run);
			map.maskingGrid.x = run.grid.x;
			map.maskingGrid.y = run.grid.y;
			map.maskingGrid.mask = run.grid.mask;
			map.siMask = scatteredInterpolant(map.maskingGrid.x(:),...
							map.maskingGrid.y(:),...
							double(map.maskingGrid.mask(:)),'nearest');
		end

		function n = t2n(map,t);
			n = interp1(map.t(:),(1:length(map.t))',t,'nearest');
		end
		
		
		% integrating trajectories ---------------------------------------------
		% (applying the map; analogous to par_integrate followed by 
		% par_concatSteps)
		
		function P = integrate(map,x0,y0,t0,t1,tracers);
			if nargin < 6, tracers = fieldnames(map.c); end
			if nargin < 5, t1 = map.t(end); end
			if nargin < 4, t0 = map.t(1); end
			if nargin < 3, y0 = map.y; end
			if nargin < 2, x0 = map.x; end
			n0 = map.t2n(t0);
			n1 = map.t2n(t1);
			nn = n0:n1;
			numpar = prod(size(x0));
			numrep = size(map.M,3); 
			numtime = length(nn);
			ind = repmat(nan,[numtime numpar]);
			% repeatedly apply the map to fill out _ind_ (indices into x,y
			% over time) and then look up all other fields afterwards
			ind(1,:) = map.xy2ind(x0,y0);
			for ni = 2:numtime
				n = nn(ni);
				rep = randi(numrep,[1 numpar]);
				i3 = sub2ind(size(map.M), repmat(n-1,[1 numpar]), ...
							   ind(ni-1,:), rep);
					% linear index corresponding to M at timestep n-1, at
					% the positions of this particular particle set (ind),
					% for a random choice among the replicates
				ind(n,:) = map.M(i3);
			end
			% intialize a structure to hold all the output
			sz = size(x0);
			P.ind = reshape(ind,[numtime sz]);
			P.x = map.x(P.ind);
			P.y = map.y(P.ind);
			% look up all fields at (n, ind(n,...))
			n2 = repmat(nn(:),[1 numpar]);
			i2 = sub2ind([length(map.t) size(map.M,2)], n2, ind);
			i2 = reshape(i2,[numtime sz]);
			for j=1:length(tracers)
				P.(tracers{j}) = map.c.(tracers{j})(i2);
			end		
		end
		
		
	end
	
end
