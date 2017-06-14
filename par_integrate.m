function steps = par_integrate(rel,run,basefilename);

% steps = par_integrate(rel,run);
% steps = par_integrate(rel,run,'basefilename');
%
% integrates the particle release experiment specified by _rel_ using the
% model run _run_.
% returns a vector of _step_ structures that can be easily broken up into 
% sequential files. If a basefilename is given, saves steps to
% basefilename1.mat, basefilename2.mat,....


if nargin<3, basefilename = ''; end
steps = [];
saveToVar = (nargout > 0);

% make sure flags in _rel_ are consistent and that undefined ones are set to
% sensible defaults
if isempty(rel.diffusive)
	rel.diffusive = 0;
end
if isempty(rel.verbose)
	rel.verbose = 0;
end
if ~isempty(rel.sigmaTrapLevel)
	rel.zTrapLevel = [];
	rel.diffusive = 0;
end
if ~isempty(rel.zTrapLevel)
	rel.diffusive = 0;
end

if ~run.nativeSigma
	if rel.diffusive
		warning('the dKsdz code hasn''t been rewritten to use the _in_z functions yet. Set rel.zTrapLevel and thus rel.diffusive = 0, use a sigma-coordinate model, or fill in this code.');
	end
	if isempty(rel.zTrapLevel)
		warning('so far the interp*_in_z() functions only handle the special case where rel.zTrapLevel is set. For the full 3D case, set rel.nativeSigma=1 to use the default functions, or write more code');
	end
end
		
if rel.verbose
	disp('integrating with particle release...');
	rel
	disp('...and model run...');
	run
end

% frame numbers (in the model run's terms) that span the particle integration
n00 = floor(min(interp1(run.t, 1:run.numFrames, rel.t0(:))));
n11 = ceil(max(interp1(run.t, 1:run.numFrames, rel.t1(:))));
nn = n00:n11;
if rel.verbose
	disp(['using frames ' num2str(n00) ' .. ' num2str(n11)]);
end

% load the initial frame
if (rel.verbose), disp('loading first frames'); end
run.loadFrame(n00,rel.tracers);
run.advanceTo(n00,rel.tracers);
% set up the initial set of particles--as if trajectories had previously been
% integrated to n00
s1.n = nn(1);
s1.t = run.t(nn(1));
s1.x = rel.x0;
s1.y = rel.y0;
if ~isempty(rel.sigma0) & isempty(rel.z0)
	s1.sigma = rel.sigma0;
	s1.z = [];
elseif isempty(rel.sigma0) & ~isempty(rel.z0)
	s1.z = rel.z0;
	s1.sigma = [];
elseif ~isempty(rel.sigma0) & ~isempty(rel.z0)
	warning('either rel.sigma0 or rel.z0 should be empty. keeping sigma0');
	s1.z = sigma2z(rel.sigma0,rel);
	s1.sigma = [];
else isempty(rel.sigma0) & isempty(rel.z0)
	if ~isempty(rel.zTrapLevel)
		s1.z = rel.zTrapLevel .* ones(size(s1.x));
		s1.sigma = [];
	elseif ~isempty(rel.sigmaTrapLevel)
		s1.sigma = rel.sigmaTrapLevel .* ones(size(s1.x));
		s1.z = [];
	else
		error('can''t find sigma0 or z0 or any other clues.');
	end
end
dt = (run.t(nn(2)) - run.t(nn(1))) / rel.Ninternal;
s1 = interpEverything(s1,dt,rel,run);
steps = saveStep(s1,1,saveToVar,steps,basefilename);


% main loop
for ni = 2:length(nn)
	run.advanceTo(nn(ni),rel.tracers);

	if rel.verbose
		disp(['step ' num2str(nn(ni)) ...
		      ', between model frames ' num2str(run.loadedN)]);
	end

	tt = run.t(run.loadedN); % model time range in memory
	Ninternal = rel.Ninternal;
	dt = diff(tt) / Ninternal; % internal timestep

	for m = 1:Ninternal
		if rel.verbose, disp('    .'); end
		s0 = s1;
		s1 = takeStep(s0,dt,rel,run);
		s1 = interpEverything(s1,dt,rel,run);
	end

	s1.n = run.loadedN(end);
	s1.t = repmat(tt(2),size(s1.t));
		% make sure particles are exactly at the time we think they're at
	steps = saveStep(s1,ni,saveToVar,steps,basefilename);
		% save to either files or memory, on the same timebase as the
		% model output itself
end



% ------------------------------------------------------------------------------
function s = interpEverything(s0,dt,rel,run);
% takes a set of particle positions s.x, s.y, s.z, s.t and interpolates
% cs, H, zeta, u, v, w, dksdz, wdiff, tracers, uScaled, vScaled, active
s = s0;

[s.x, s.y, s.active] = run.filterCoordinates(s.x, s.y);
s.active = s.active & (s0.t >= rel.t0);

s.H = run.interpH(s.x, s.y);
s.zeta = run.interpZeta(s.x, s.y, s.t);
s.mask = run.interpMask(s.x, s.y, s.t);
if ~isempty(rel.zTrapLevel) % z trapped
	s.z = rel.zTrapLevel;
	s.sigma = z2sigma(s.z, s.H, s.zeta);
elseif ~isempty(rel.sigmaTrapLevel) % sigma trapped
	s.sigma = rel.sigmaTrapLevel;
	s.z = sigma2z(s.sigma, s.H, s.zeta);
elseif isempty(s.z) % z not defined yet
	s.z = sigma2z(s.sigma, s.H, s.zeta);	
else % normal case
	s.sigma = z2sigma(s.z, s.H, s.zeta);
	s.z = sigma2z(s.sigma, s.H, s.zeta);	
end

if run.nativeSigma
	s.u = run.interpU(s.x, s.y, s.sigma, s.t);
	s.v = run.interpV(s.x, s.y, s.sigma, s.t);
	s.w = run.interpW(s.x, s.y, s.sigma, s.t);
	s.Ks = run.interpKs(s.x, s.y, s.sigma, s.t);
	for i=1:length(rel.tracers)
		s.(rel.tracers{i}) = run.interpTracer(rel.tracers{i}, ...
								s.x, s.y, s.sigma, s.t);
	end
else % z levels
	s.u = run.interpU_in_z(s.x, s.y, s.z, s.t);
	s.v = run.interpV_in_z(s.x, s.y, s.z, s.t);
	s.w = run.interpW_in_z(s.x, s.y, s.z, s.t);
	s.Ks = run.interpKs_in_z(s.x, s.y, s.z, s.t);
	for i=1:length(rel.tracers)
		s.(rel.tracers{i}) = run.interpTracer_in_z(rel.tracers{i}, ...
								s.x, s.y, s.z, s.t);
	end	
end

s.uScaled = run.scaleU(s.u, s.x, s.y);
s.vScaled = run.scaleV(s.v, s.x, s.y);

if rel.diffusive
	dt_secs = dt .* 86400; % assumes w is in (z units) per sec,
						   % Ks is in (z units)^2 per sec
	% diffusion gradient dKs/dz
	wdiff_approx = sqrt(2.*s.Ks./dt_secs);
	dsigma = wdiff_approx .* dt_secs ./ (s.H + s.zeta);
		% half-span to take gradient over--the scale of the next diffusive step
	sigmatop = min(s.sigma + dsigma,0);
	sigmabot = max(s.sigma - dsigma,-1);
	Kstop = run.interpKs(s.x, s.y, sigmatop, s.t);
	Ksbot = run.interpKs(s.x, s.y, sigmabot, s.t);
	s.dKsdz = (Kstop - Ksbot) ./ (sigmatop - sigmabot) ./ (s.H + s.zeta);
	% diffusion velocity wdiff
	sigma1 = z2sigma(s.z + 0.5.*s.dKsdz.* dt_secs, s.H, s.zeta);
	Ks1 = run.interpKs(s.x, s.y, sigma1, s.t);
	s.wdiff = sqrt(2.*Ks1./dt_secs) .* randn(size(Ks1));
else
	s.dKsdz = 0;
	s.wdiff = 0;
end


% ------------------------------------------------------------------------------
function s = z2sigma(z, H, zeta);
if nargin < 3, zeta = 0; end
s = (z - zeta) ./ (H + zeta);
s = min(max(s,-1),0);

function z = sigma2z(sigma, H, zeta);
if nargin < 3, zeta = 0; end
z = min(max(sigma,-1),0) .* (H + zeta) + zeta;


% ------------------------------------------------------------------------------
function s1 = takeStep(s0,dt,rel,run);
% the basic operation X1 = X0 + X*dt.
% midpoint method.
% fills in only x,y,z,t; other fields are calculated in interpEverything().
ac = double(s0.active); % when this is 0, x,y,z do not advance but t does
smid.x = s0.x + ac .* s0.uScaled .* 0.5 .* dt; % take half an advective step
smid.y = s0.y + ac .* s0.vScaled .* 0.5 .* dt;
smid.z = s0.z + ac .* (s0.w .* run.wScaleFactor) .* 0.5 .* dt;
smid.t = s0.t + 0.5 .* dt;
smid = interpEverything(smid,dt,rel,run); % calculate new advective velocities
s1.x = s0.x + ac .* smid.uScaled .* dt; % full step
s1.y = s0.y + ac .* smid.vScaled .* dt;
s1.z = s0.z + ac .* (smid.w + s0.wdiff + s0.dKsdz) .* run.wScaleFactor .* dt;
s1.t = s0.t + dt;


% ------------------------------------------------------------------------------
function steps1 = saveStep(step,i,saveToVar,steps,basefilename);
% saves _step_ either to the variable _steps_, a numbered file, or both.

if saveToVar
	if ~isempty(steps)
		steps1 = steps;
		steps1(i) = step;
	else
		steps1 = step;
	end
end
if ~isempty(basefilename)
	nstr = ['0000' num2str(i)];
	nstr = nstr(end-3:end);
	filename = [basefilename nstr '.mat'];
	save(filename,'step');
end


