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
rel.verticalMode = '3D'; % 3D | zLevel | sigmaLevel | zAverage
rel.verticalLevel = []; % if mode is 3D, this is ignored;
						% if zLevel or sigmaLevel, this should be a scalar;
						% if zAverage, it should be [min max]
rel.verticalDiffusion = 1; % this is only allowed if verticalMode is 3D
rel.tracers = {}; % tracers to save as scalar time series
rel.profiles = {}; % tracers to save as complete vertical profiles
rel.Ninternal = 2; % internal timesteps per frame of saved model output
rel.parallel = 0; % serial calculation unless this is turned on
rel.verbose = 0; % quiet about status updates unless this is turned on


% overwrite with specified name-value pairs
fields = {varargin{1:2:end}};
vals = {varargin{2:2:end}};
for i=1:length(fields)
	rel.(fields{i}) = vals{i};
end


% make sure everything is consistent. Anything that requires the model run
% itself (like harmonizing sigma0 and z0) waits until par_integrate.m.
if strcmpi(rel.verticalMode,'3D')
	rel.verticalLevel = [];
else
	rel.verticalDiffusion = 0;
	if strcmpi(rel.verticalMode,'sigmaLevel')
		rel.verticalLevel = min(max(rel.verticalLevel,-1),0);
	elseif strcmpi(rel.verticalMode,'zAverage')
		rel.verticalLevel = sort(rel.verticalLevel);
	end
end
if length(rel.sigma0)==1
	rel.sigma0 = repmat(rel.sigma0,size(rel.x0));
	% could handle more special cases along these lines--this one seems to be
	% the most common
end

