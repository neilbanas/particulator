function steps = par_integrate(rel,run,basefilename)

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

% setup ------------------------------------------------------------------------
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
tracersAndProfiles = union(rel.tracers,rel.profiles);
run.loadFrame(n00,tracersAndProfiles);
run.advanceTo(n00,tracersAndProfiles);
% set up the initial set of particles--as if trajectories had previously been
% integrated to n00
s1.n = nn(1);
s1.t = run.t(nn(1));
s1.x = rel.x0;
s1.y = rel.y0;
% this is our first chance to harmonise the way sigma0 and z0 were specified
if strcmpi(rel.verticalMode,'3D')	
	if ~isempty(rel.sigma0)
		s1.sigma = rel.sigma0;
		s1.z = [];
		if ~isempty(rel.z0)
			warning('either sigma0 or z0 should be empty. Keeping sigma0');
		end
	elseif ~isempty(rel.z0)
		s1.z = rel.z0;
		s1.sigma = [];
	else
		error('can''t find sigma0 or z0 or any other clues.');
	end
	% for modes other than 3D, z and sigma are overwritten in interpEverything
	% so they don't need to be defined here
end
dt = (run.t(nn(2)) - run.t(nn(1))) / rel.Ninternal;
s1 = interpEverything(s1,dt,rel,run);
steps = saveStep(s1,1,saveToVar,steps,basefilename);

% number of input time steps per output time steps (for temporal averaging)
try
    dtIn_per_dtOut = run.dtOut;
    dtIn_per_dtOut(:,2) = run.dtOut(:,2)./run.dtIn(:,2);
catch
    dtIn_per_dtOut = ones(1,2);
end
if  any(abs(dtIn_per_dtOut(:,2)/round(dtIn_per_dtOut(:,2))-1)>1.e-6)
    i = find(abs(dtIn_per_dtOut(:,2)/round(dtIn_per_dtOut(:,2))-1)>1.e-6, 1, 'first');
    error('Input file #%i - Selected output time step (%.4f) no proper multiple of input time step (%.4f).', ...
          run.dtOut(i,:), run.dtIn(i,2));
end
dtIn_per_dtOut = round(dtIn_per_dtOut);

% main loop --------------------------------------------------------------------
nOff = 1; % averaging/storage index offset
tPassed = 0; % simulation time elapsed since last temporal averaging
nStep = 0;
eps = 1.e-6;
for ni = 2:length(nn)
	run.advanceTo(nn(ni),tracersAndProfiles);

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
    
    % increase elapsed time and step counter
    tPassed = tPassed + diff(tt);
    nStep = nStep + 1;
    steps = saveStep(s1,nOff+nStep,saveToVar,steps,basefilename);
		% save to either files or memory, on the same timebase as the
		% model output itself
    
    % average over time steps
    %  - do this with one step delay for correct calculation of first step
    %    of next averaging interval
    %  - updated 'steps' only contains initial conditions, step averages
    %    and first step of next averaging interval
    if  isfield(run, 'fileind')
        idt = dtIn_per_dtOut(:,1) == run.fileind(run.loadedN(1));
    else
        idt = 1;
    end
    if  (dtIn_per_dtOut(idt,2)>1) && ((tPassed-run.dtOut(idt,2)>eps) || (ni==length(nn)))
        if  ni==length(nn)
            nSteps = nStep;
        else
            nSteps = nStep - 1;
        end
        ia = nOff + (1:nSteps); % indeces of steps for averaging
        steps = averageSteps(steps,nOff+1,ia);
        nOff = nOff + 1; % increase storage index offset
        nStep = 1;       % reset step counter
        tPassed = tPassed - run.dtOut(idt,2); % reset elapsed time
    elseif (dtIn_per_dtOut(idt,2)==1)
        % no averaging => only increase storage index offset
        nOff = nOff + 1;
        nStep = 0;
    end
end

% in case of temporal averaging, remove last step (remnant of delayed averaging)
if dtIn_per_dtOut(idt,2)>1, steps = steps(1:end-1); end

% end of par_integrate

% ------------------------------------------------------------------------------
function s1 = takeStep(s0,dt,rel,run)
% the basic operation X1 = X0 + X*dt.
% midpoint method.
% fills in only x,y,z,t; other fields are calculated in interpEverything().
ac = double(s0.active); % when this is 0, x,y,z do not advance but t does
smid.x = s0.x + ac .* s0.uScaled .* 0.5 .* dt; % take half an advective step
smid.y = s0.y + ac .* s0.vScaled .* 0.5 .* dt;
smid.z = s0.z + ac .* s0.wScaled .* 0.5 .* dt;
smid.t = s0.t + 0.5 .* dt;
smid = interpEverything(smid,dt,rel,run); % calculate new advective velocities
s1.x = s0.x + ac .*  smid.uScaled .* dt; % full step
s1.y = s0.y + ac .*  smid.vScaled .* dt;
s1.z = s0.z + ac .* (smid.wScaled + s0.wdiff + s0.dKsdz) .* dt;
s1.t = s0.t + dt;
if rel.avoidLand
	% if the step is going to take the particle into a region where the
	% interpolated land mask is less than 0.5, don't take the step
	mask1 = run.interp('mask',s1.x,s1.y,[],s1.t);
	f = find(mask1 < rel.landThreshhold);
	if ~isempty(f)
		s1.x(f) = s0.x(f);
		s1.y(f) = s0.y(f);
		s1.z(f) = s0.z(f);
	end
end

% ------------------------------------------------------------------------------
function s = interpEverything(s0,dt,rel,run)
% takes a set of particle positions s.x, s.y, s.z, s.t and interpolates
% sigma, H, zeta, u, v, w, dksdz, wdiff, tracers, uScaled, vScaled, active
s = s0;

[s.x, s.y, s.active] = run.filterCoordinates(s.x, s.y);
s.active = s.active & (s0.t >= rel.t0);

s.H = run.interp('H',s.x, s.y);
s.zeta = run.interp('zeta',s.x, s.y, [], s.t);
s.mask = run.interp('mask',s.x, s.y, [], s.t);

if strcmpi(rel.verticalMode,'zLevel')
	s.z = repmat(rel.verticalLevel,size(s.x));
	s.sigma = z2sigma(s.z, s.H, s.zeta);
elseif strcmpi(rel.verticalMode,'sigmaLevel')
	s.sigma = repmat(rel.verticalLevel,size(s.x));
	s.z = sigma2z(s.sigma, s.H, s.zeta);
elseif strcmpi(rel.verticalMode,'zAverage')
	s.z = repmat(mean(rel.verticalLevel),size(s.x));
	s.sigma = z2sigma(s.z, s.H, s.zeta);
else % 3D
	if isempty(s.z) % z not defined yet, perhaps at the first step
		s.z = sigma2z(s.sigma, s.H, s.zeta);	
	else % normal case
		s.sigma = z2sigma(s.z, s.H, s.zeta);
		s.z = sigma2z(s.sigma, s.H, s.zeta);	
	end
end

if strcmpi(rel.verticalMode,'zAverage')
	s.u = run.interpDepthAverage('u', s.x, s.y, rel.verticalLevel, s.t);
	s.v = run.interpDepthAverage('v', s.x, s.y, rel.verticalLevel, s.t);
	s.w = run.interpDepthAverage('w', s.x, s.y, rel.verticalLevel, s.t);
else
	s.u = run.interp('u', s.x, s.y, s.sigma, s.t);
	s.v = run.interp('v', s.x, s.y, s.sigma, s.t);
	s.w = run.interp('w', s.x, s.y, s.sigma, s.t);
end
s.uScaled = run.scaleU(s.u, s.x, s.y);
s.vScaled = run.scaleV(s.v, s.x, s.y);
s.wScaled = run.scaleW(s.w);
if strcmpi(rel.verticalMode,'zAverage')
	s.Ks = run.interpDepthAverage('Ks', s.x, s.y, rel.verticalLevel, s.t);
	for i=1:length(rel.tracers)
		s.(rel.tracers{i}) = run.interpDepthAverage(rel.tracers{i}, ...
								s.x, s.y, rel.verticalLevel, s.t);
	end
else
	s.Ks = run.interp('Ks', s.x, s.y, s.sigma, s.t);
	for i=1:length(rel.tracers)
		s.(rel.tracers{i}) = run.interp(rel.tracers{i}, ...
								s.x, s.y, s.sigma, s.t);
	end
end

for i=1:length(rel.profiles)
	if i==1
		s.profiles.v_axis = run.verticalAxisForProfiles;
	end
	s.profiles.(rel.profiles{i}) = run.interpProfile(rel.profiles{i}, ...
									 s.x, s.y, s.t);
end

if rel.verticalDiffusion
	dt_secs = dt .* 86400;
		% advective velocity is stored as m/day, but this block of code works in
		% m/s and m^2/s, and converts to m/day at the end
	% diffusion gradient dKs/dz
	wdiff_approx = sqrt(2.*s.Ks./dt_secs);
	dsigma = wdiff_approx .* dt_secs ./ (s.H + s.zeta);
		% half-span to take gradient over--the scale of the next diffusive step
	sigmatop = min(s.sigma + dsigma,0);
	sigmabot = max(s.sigma - dsigma,-1);
	Kstop = run.interp('Ks', s.x, s.y, sigmatop, s.t);
	Ksbot = run.interp('Ks', s.x, s.y, sigmabot, s.t);
	s.dKsdz = (Kstop - Ksbot) ./ (sigmatop - sigmabot) ./ (s.H + s.zeta);
	% diffusion velocity wdiff
	sigma1 = z2sigma(s.z + 0.5.*s.dKsdz.* dt_secs, s.H, s.zeta);
	Ks1 = run.interp('Ks', s.x, s.y, sigma1, s.t);
	s.wdiff = sqrt(2.*Ks1./dt_secs) .* randn(size(Ks1));
	% now put both velocity terms in m/day
	s.dKsdz = s.dKsdz .* 86400;
	s.wdiff = s.wdiff .* 86400;
else
	s.dKsdz = 0;
	s.wdiff = 0;
end

if rel.horizIsotropicDiffusivity > 0
	dt_secs = dt .* 86400;
	udiffscale = sqrt(2.*rel.horizIsotropicDiffusivity./dt_secs); 
	s.u = s.u + udiffscale .* randn(size(s.u));
	s.v = s.v + udiffscale .* randn(size(s.u));
	s.uScaled = run.scaleU(s.u, s.x, s.y);
	s.vScaled = run.scaleV(s.v, s.x, s.y);
end

if rel.horizShearDiffusion
	% assume verticalMode = zAverage 
	ualong = sqrt(s.u.^2 + s.v.^2);
	rx = s.u./ualong; % unit vector in the along-flow direction is
	ry = s.v./ualong; % rx xhat + ry yhat
	ualongdiff = ualong ./ sqrt(3) .* randn(size(s.u));
	s.u = s.u + ualongdiff .* rx;
	s.v = s.v + ualongdiff .* ry;
	s.uScaled = run.scaleU(s.u, s.x, s.y);
	s.vScaled = run.scaleV(s.v, s.x, s.y);	
end



% ------------------------------------------------------------------------------
function s = z2sigma(z, H, zeta)
if nargin < 3, zeta = 0; end
s = (z - zeta) ./ (H + zeta);
s = min(max(s,-1),0);

function z = sigma2z(sigma, H, zeta)
if nargin < 3, zeta = 0; end
z = min(max(sigma,-1),0) .* (H + zeta) + zeta;



% ------------------------------------------------------------------------------
function steps1 = saveStep(step,i,saveToVar,steps,basefilename)
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



% -------------------------------------------------------------------------
function steps1 = averageSteps(steps0, nAve, iAve)
% temporally average selected steps
steps  = steps0(iAve);     % steps used for averaging
steps1 = steps0(1:nAve+1); % structure with size of updated structure
step = steps(1);           % structure for new average step
flds = fieldnames(step);   % field names of step structure
n = length(iAve);          % 
% average all fields
for i = 1:length(flds)
    for j = 2:n
        if  islogical(steps(1).(flds{i}))
            % convert logicals to doubles for arithmetics
            step.(flds{i}) = double(step.(flds{i})) + double(steps(j).(flds{i}));
        elseif isstruct(steps(1).(flds{i}))
            % treat fields of sub-structures
            sflds = fieldnames(steps(1).(flds{i}));
            for k = 1:length(sflds)
                step.(flds{i}).(sflds{k}) = step.(flds{i}).(sflds{k}) + ...
                                            steps(j).(flds{i}).(sflds{k});
            end
        else
            % numerical values (step doesn't have char fields)
            step.(flds{i}) = step.(flds{i}) + steps(j).(flds{i});
        end
    end
    if  islogical(steps(1).(flds{i}))
        % convert logical fields back to logical
        step.(flds{i}) = step.(flds{i}) > 0.5;
    elseif isstruct(steps(1).(flds{i}))
        % treat fields of sub-structures
        sflds = fieldnames(steps(1).(flds{i}));
        for k = 1:length(sflds)
            step.(flds{i}).(sflds{k}) = step.(flds{i}).(sflds{k})/n;
        end
    else
        % numerical values
        step.(flds{i}) = step.(flds{i})/n;
    end
end
step.n = nAve;       % n is set to new position in averaged structure
step.t = steps(n).t; % t is set to end of averaging period
steps1(nAve) = step; % update step structure with averaged step
steps1(nAve+1) = steps0(end); % keep first step of next averaging interval;
                              % needed for correct calculation of next step