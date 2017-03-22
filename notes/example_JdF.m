% an example particle release in the spirit of Giddings et al. JGR 2014,
% Figure 2

% locate a model run to use as source data
modelDir = '/Users/neil/Dropbox/util/flow samples/cascadia2004/';
run = modelRun_romsCascadia(modelDir);
disp(['found ' num2str(run.numFrames) ' model saves from ' ...
	  datestr(run.t(1)) ' to ' datestr(run.t(end))]);
disp(['the first file is ' run.filename{1}]);
							
% pick initial coordinates
% space
x0 = -125.8 : 0.02 : -125;
y0 = 48.2 : 0.02 : 48.8; % grid in the JdF Eddy region
sigma0 = [-0.75 0]; % surface and (fairly) deep
[x0,y0,sigma0] = ndgrid(x0(:),y0(:),sigma0(:)); % make a 3D grid out of those
H0 = run.interpH(x0,y0); % get bottom depth from the model run
x0 = x0(H0>50); % keep only points beyond the 50 m isobath
y0 = y0(H0>50);
sigma0 = sigma0(H0>50); 
% time
t0 = run.t(1) .* ones(size(x0)); % start at first model save
t1 = run.t(end) .* ones(size(x0));
t1 = min(t1, t0 + 10); % track for 20 days or until the last save

% set timestep
DT_saves = run.t(2) - run.t(1);
Ninternal = 24; % internal timesteps for particle integration per interval
			   % between saved frames in the model run
disp(['integrating ' num2str(length(x0(:))) ' particles with a timestep of ' ...
	  num2str(DT_saves/Ninternal) ' days']);


% this returns an object specifying the setup of a particle release,
% summarizing all the choices above
rel = par_release('x0',x0,'y0',y0,'sigma0',sigma0,'t0',t0,'t1',t1,...
				  'Ninternal',Ninternal,...
				  'tracers',{'salt','temp'});
rel.verbose = 1; % extra diagnostic info please

							
% the actual Lagrangian integration, saved in _steps_
% this can easily be saved to sequentiol files instead by adding a basename
% as a third argument
steps = par_integrate(rel,run);
		
							
% concatenates the fields in _steps_ variable by variable into a single 
% structure P that's easier to analyze
P = par_concatSteps(steps); 


% example figures
figure
fsurf = find(P.sigma(1,:)==0);
fdeep = find(P.sigma(1,:)==-0.75);
plot(P.x(:,fsurf),P.y(:,fsurf),'b'); % full trajectories for surface particles
hold on
plot(P.x(:,fdeep),P.y(:,fdeep),'k'); % trajectories for deep particles
hold on
plot(P.x(1,:),P.y(1,:),'g.'); % start locations
plot(P.x(end,:),P.y(end,:),'r.'); % end locations
contour(run.grid.lon,run.grid.lat,run.grid.mask,[0.5 0.5],'k'); % coastline

figure
subplot 311
plot(P.t,P.z(:,fsurf),'b',P.t,P.z(:,fdeep),'k'); % depth vs time
ylabel('depth (m)');
datetick('x','mmm dd');

subplot 312
plot(P.t,P.salt(:,fsurf),'b',P.t,P.salt(:,fdeep),'k'); % salinity vs time
ylabel('salinity (psu)');
datetick('x','mmm dd');
