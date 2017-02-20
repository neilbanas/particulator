function rel = par_release(varargin);

% r = par_release('param1',val1,'param2',val2,...);
%
% specification for a particle release experiment.

% defaults
rel.x0 = [];
rel.y0 = [];
rel.z0 = [];
rel.sigma0 = [];
rel.t0 = [];
rel.t1 = [];
rel.sigmaTrapLevel = [];
rel.zTrapLevel = [];
rel.diffusive = 1;
rel.tracers = {};
rel.dt_per_DT = 2; % internal timesteps per frame of saved model output

% overwrite with specified name-value pairs
fields = varargin{1:2:end};
vals = varargin{2:2:end};
for i=1:length(fields)
	rel.(fields{i}) = val{i};
end

% make sure everything is consistent. Anything that requires the model run
% itself (like harmonizing cs0 and z0) waits until par_integrate.m.
if ~isempty(rel.sigmaTrapLevel) & ~isempty(rel.zTrapLevel)
	warning(...
		'either sigmaTrapLevel or zTrapLevel should be empty. Keeping sigma.');
	rel.zTrapLevel = [];
end
if ~isempty(rel.sigmaTrapLevel) | ~isempty(rel.zTrapLevel)
	rel.diffusive = 0;
end