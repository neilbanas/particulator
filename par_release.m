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
rel.verticalLevel = [];
	% if mode is 3D, this is ignored.
	% if zLevel or sigmaLevel, this should be a scalar, or else empty
	%    (in which case each particle is trapped at its initial sigma value).
	% if zAverage, it should be [min max]
rel.verticalDiffusion = 1; % this is only allowed if verticalMode is 3D
rel.horizIsotropicDiffusivity = 0; % if 0, non-diffusive; otherwise,
								   % horizontal diffusivity, in units that
								   % match the native velocity units:
								   % [velocity]^2 * sec
rel.horizShearDiffusion = 0; % if 1, adds a diffusivity in the direction
							 % of the depth-average flow proportional to
							 % the depth-average speed (only makes sense
							 % when using zAverage)
rel.tracers = {}; % tracers to save as scalar time series
rel.profiles = {}; % tracers to save as complete vertical profiles
rel.Ninternal = 2; % internal timesteps per frame of saved model output
rel.avoidLand = 1; % if 1, don't take steps that would place the particle
				   % at a point where the interpolated mask < landThreshhold
rel.landThreshhold = 0.5;
rel.verbose = 0; % quiet about status updates unless this is turned on


% overwrite with specified name-value pairs
fields = {varargin{1:2:end}};
vals = {varargin{2:2:end}};
for i=1:length(fields)
	rel.(fields{i}) = vals{i};
end


% make sure everything is consistent. Anything that requires the model run
% itself (like filling in sigma0 from z0 or vice versa) waits until 
% par_integrate.m.
if strcmpi(rel.verticalMode,'3D')
	rel.verticalLevel = [];
	rel.horizShearDiffusion = 0;
else
	rel.verticalDiffusion = 0;
	if strcmpi(rel.verticalMode,'sigmaLevel')
		if ~isempty(rel.verticalLevel)
			rel.sigma0 = min(max(rel.verticalLevel,-1),0);
		else
			rel.verticalLevel = min(max(rel.sigma0,-1),0);
		end
		rel.z0 = [];
		rel.horizShearDiffusion = 0;
	elseif strcmpi(rel.verticalMode,'zLevel')
		if ~isempty(rel.verticalLevel)
			rel.z0 = rel.verticalLevel;
		else
			rel.verticalLevel = rel.z0;
		end
		rel.sigma0 = [];
		rel.horizShearDiffusion = 0;
	elseif strcmpi(rel.verticalMode,'zAverage')
		rel.verticalLevel = sort(rel.verticalLevel);
		rel.z0 = mean(rel.verticalLevel);
		rel.sigma0 = [];
	else
		error(['don''t recognise the verticalMode ''' rel.verticalMode '''.']);
	end
end
if length(rel.sigma0)==1
	rel.sigma0 = repmat(rel.sigma0,size(rel.x0));
end
if length(rel.z0)==1
	rel.z0 = repmat(rel.z0,size(rel.x0));
end

