function steps = par_integrate(rel,run,basefilename);

% steps = par_integrate(rel,run);
% steps = par_integrate(rel,run,'basefilename');
%
% integrates the particle release experiment specified by _rel_ using the
% model run _run_.
% returns a vector of _step_ structures that can be easily broken up into 
% sequential files. If a basefilename is given, saves steps to
% basefilename1.mat, basefilename2.mat,....

parallel = 0;

if nargin<2, basefilename = ''; end
steps = [];
saveToVar = nargout > 0;

% frame numbers (in the model run's terms) that span the particle integration
n00 = floor(min(interp1(run.t, 1:run.numFrames, rel.t0(:))));
n11 = ceil(max(interp1(run.t, 1:run.numFrames, rel.t1(:))));
nn = n00:n11;

% load the initial frame
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
elseif isempty(rel.sigma0) & ~isempty(rel.z0)
	s1.z = rel.z0;
elseif ~isempty(rel.sigma0) & ~isempty(rel.z0)
	warning('either sigma0 or z0 should be empty. keeping sigma0');
	s1.z = cs2z(rel.sigma0,rel);
else isempty(rel.sigma0) & isempty(rel.z0)
	if ~isempty(rel.zTrapLevel)
		s1.z = rel.zTrapLevel .* ones(size(s1.x));
	elseif ~isempty(rel.sigmaTrapLevel)
		s1.sigma = rel.sigmaTrapLevel .* ones(size(s1.x));
	else
		error('can''t find sigma0 or z0 or any other clues.');
	end
end
dt = (run.t(nn(2)) - run.t(nn(1))) / rel.dt_per_DT;
s1 = interpEverything(s1,dt,rel,run);
s1.active = isactive(s1);
saveStep(s1,1,saveToVar,steps,basefilename);

% main loop
for ni = 2:length(nn)
	run.advanceTo(nn(ni),rel.tracers);

	if verbose
		disp(nn(ni));
	end

	tt = run.t(run.loadedN); % model time range in memory
	Ninternal = rel.dt_per_DT;
	dt = diff(tt) / Ninternal; % internal timestep

	if parallel % --------------------------
		s1blocks = disassemble(s1);
		parfor j=1:length(s1blocks)
			s1blockj = s1blocks{j};
			for m = 1:Ninternal
				s0blockj = s1blockj;
				s1blockj = takeStep(s0blockj,dt,rel,run);
				s1blockj = interpEverything(s1blockj,dt,rel,run);
			end
			s1blocks{j} = s1blockj;	
		end
		s1 = reassemble(s1blocks);
		
	else % ---------------------------------
		for m = 1:Ninternal
			s0 = s1;
			s1 = takeStep(s0,dt,rel,run);
			s1 = interpEverything(s1,rel,run);
		end
	end
	
	s1.t = repmat(tt(2),size(s1.t));
		% make sure particles are exactly at the time we think they're at
	saveStep(s1,ni,saveToVar,steps,basefilename);
		% save to either files or memory, on the same timebase as the
		% model output itself
end



% ------------------------------------------------------------------------------
function s = interpEverything(s0,dt,rel,run);
% takes a set of particle positions s.x, s.y, s.z, s.t and interpolates
% cs, H, zeta, u, v, w, dksdz, wdiff, tracers, uScaled, vScaled, active
s = s0;
s.H = run.interpH(s.y, s.x);
s.zeta = run.interpZeta(s.t, s.y, s.x);
s.mask = run.interpMask(s.t, s.y, s.x);
if ~isempty(rel.zTrapLevel)
	s.z = rel.zTrapLevel;
	s.sigma = z2sigma(s.z, s.H, s.zeta);
elseif ~isempty(rel.sigmaTrapLevel)
	s.sigma = rel.sigmaTrapLevel;
	s.z = sigma2z(s.sigma, s.H, s.zeta);
else
end
s.u = run.interpU(s.t, s.sigma, s.y, s.x);
s.v = run.interpV(s.t, s.sigma, s.y, s.x);
s.w = run.interpW(s.t, s.sigma, s.y, s.x);
s.uScaled = run.scaleU(s.u, s.y, s.x);
s.vScaled = run.scaleV(s.v, s.y, s.x);
for i=1:length(rel.tracers)
	s.(rel.tracers{i}) = run.interpTracer(rel.tracers{i}, ...
							s.t, s.sigma, s.y, s.x);
end
if rel.diffusive
	% diffusion gradient dKs/dz
	dz = 1; % half-span to take gradient over. Not a good method, but preserving
			% particulator-java behaviour for now.
	dsigma = dz ./ (s.H + s.zeta);
	sigmatop = min(s.sigma + dsigma,0);
	sigmabot = max(s.sigma - dsigma,-1);
	Kstop = run.interpKs(s.t, sigmatop, s.y, s.x);
	Ksbot = run.interpKs(s.t, sigmabot, s.y, s.x);
	s.dKsdz = (Kstop - Ksbot) ./ (sigmatop - sigmabot) ./ (s.H + s.zeta);
	% diffusion velocity wdiff
	sigma1 = z2sigma(s.z + 0.5.*s.dKsdz.* dt, s.H, s.zeta);
	Ks1 = run.interpKs(s.t, sigma1, s.y, s.x);
	s.wdiff = sqrt(2.*Ks1./dt) .* randn(size(Ks1));
else
	s.dKsdz = 0;
	s.wdiff = 0;
end
s.active = run.in_xy_bounds(s.y, s.x);


% ------------------------------------------------------------------------------
function s = z2sigma(z, H, zeta);
if nargin < 3, zeta = 0; end
s = (z - zeta) ./ (H + zeta);
s = min(max(s,-1),0);

function z = sigma2z(sigma, H, zeta);
if nargin < 3, zeta = 0; end
z = min(max(sigma,-1),0) .* (H + zeta) + zeta;


% ------------------------------------------------------------------------------
function s1 = takeStep(s0,dt,rel);
% the basic operation X1 = X0 + X*dt.
% midpoint method.
% fills in only x,y,z,t; other fields are calculated in interpEverything().
smid.x = s0.x + s0.uScaled .* 0.5 .* dt; % take half an advective step
smid.y = s0.y + s0.vScaled .* 0.5 .* dt;
smid.z = s0.z + s0.w       .* 0.5 .* dt;
smid.t = s0.t + 0.5 .* dt;
smid = interpEverything(smid,rel); % calculate new advective velocities
s1.x = s0.x + smid.uScaled 					 .* dt; % full step
s1.y = s0.y + smid.vScaled                   .* dt;
s1.z = s0.z + (smid.w + s0.wdiff + s0.dksdz) .* dt;
s1.t = s0.t + dt;


% ------------------------------------------------------------------------------
function blocks = disassemble(s);
% splits the step _s_ into a number of blocks for parallel integration.
NP = length(s.x);
NB = gcp.NumWorkers;
lim = [0 (1:NB).*round(NP/NB)]; % block i runs from lim(i-1)+1 to lim(i)
fields = fieldnames(s);
for i=1:length(fields)
	if length(s.(fields{i})) == NP
		for j=1:NB
			blocks{j}.(fields{i}) = s.(fields{i})(lim(j-1)+1:lim(j));
		end
	else
		blocks{j}.(fields{i}) = s.(fields{i});
	end
end

function s = reassemble(blocks);
% reconcatenates _blocks_ into a single step _s_
NB = length(blocks);
fields = fieldnames(blocks{1});
s = blocks{1};
for j=2:length(blocks)
	NP = length(blocks{j}.x);
	for i=1:length(fields)
		if length(blocks{j}.(fields{i})) == NP
			s.(fields{i}) = concat(1,s.(fields{i})(:),blocks{j}.(fields{i})(:));
		else
			% ignore fields that aren't the same size as x; presumably these are
			% scalars and the value from block 1 is fine
		end
	end
end



% ------------------------------------------------------------------------------
function steps = saveStep(step,i,saveToVar,steps,basefilename);
% saves _step_ either to the variable _steps_, a numbered file, or both.

if saveToVar
	steps(i) = step;
end
if ~isempty(basefilename)
	nstr = ['0000' num2str(i)];
	nstr = nstr(end-3:end);
	filename = [basefilename nstr '.mat'];
	save(filename,fieldnames(step));
end
