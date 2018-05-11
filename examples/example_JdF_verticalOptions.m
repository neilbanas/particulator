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
x0 = -125.8 : 0.04 : -125;
y0 = 48.2 : 0.04 : 48.8; % grid in the JdF Eddy region
[x0,y0] = meshgrid(x0(:),y0(:));
H0 = run.interp('H',x0,y0); % get bottom depth from the model run
x0 = x0(H0>50); % keep only points beyond the 50 m isobath
y0 = y0(H0>50);
% time
t0 = run.t(1) .* ones(size(x0)); % start at first model save
t1 = run.t(end) .* ones(size(x0));
t1 = min(t1, t0 + 10); % track for 10 days or until the last save

% set timestep
DT_saves = run.t(2) - run.t(1);
Ninternal = 4; % internal timesteps for particle integration per interval
			   % between saved frames in the model run
disp(['integrating ' num2str(length(x0(:))) ' particles with a timestep of ' ...
	  num2str(DT_saves/Ninternal) ' days']);


% a few variants of surface-trapped particles:
% surface-trapped
rel_s0 = par_release('x0',x0,'y0',y0,'t0',t0,'t1',t1,'Ninternal',Ninternal,...
				  'tracers',{'salt','temp'},'verbose',1,...
				  'verticalMode','sigmaLevel','verticalLevel',0);
% depth-trapped to MSL, avoiding a 1-cell buffer zone around the land
rel_z0 = par_release('x0',x0,'y0',y0,'t0',t0,'t1',t1,'Ninternal',Ninternal,...
				  'tracers',{'salt','temp'},'verbose',1,...
				  'landThreshhold',1,...
				  'verticalMode','zLevel','verticalLevel',0);
% depth-trapped to MSL, not avoiding land
rel_z0_l = par_release('x0',x0,'y0',y0,'t0',t0,'t1',t1,'Ninternal',Ninternal,...
				  'tracers',{'salt','temp'},'verbose',1,...
				  'avoidLand',0,...
				  'verticalMode','zLevel','verticalLevel',0);

% slightly deeper:
% depth-trapped to MSL - 5m
rel_z5 = par_release('x0',x0,'y0',y0,'t0',t0,'t1',t1,'Ninternal',Ninternal,...
				  'tracers',{'salt','temp'},'verbose',1,...
				  'verticalMode','zLevel','verticalLevel',-5);
% depth average, -10..0
rel_zav = par_release('x0',x0,'y0',y0,'t0',t0,'t1',t1,'Ninternal',Ninternal,...
				  'tracers',{'salt','temp'},'verbose',1,...
				  'verticalMode','zLevel','verticalLevel',[-10 0]);

% don't bother actually tracking all of them--but track these two for
% comparison
Pz0 = par_concatSteps(par_integrate(rel_z0,run));
Pz0l = par_concatSteps(par_integrate(rel_z0_l,run));

figure
plot(Pz0.x(1,:),Pz0.y(1,:),'k.');
hold on
plot(Pz0.x(end,:),Pz0.y(end,:),'o');
plot(Pz0l.x(end,:),Pz0l.y(end,:),'o');
legend('initial','final, avoiding land','final, not avoiding land');




