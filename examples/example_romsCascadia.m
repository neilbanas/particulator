modelDir = '/Users/neil/Dropbox/util/flow samples/cascadia2004/';
run = modelRun_romsCascadia(modelDir);
							% _run_ identifies a model run to use as source data

x0 = run.grid.lon;			% pick initial coordinates--
y0 = run.grid.lat;			% here, all lon-lat points in the grid.
sigma0 = 0;					% at the surface
t0 = run.t(1);				% start at the start of the model run
t1 = t0 + 3;				% end three days later
rel = par_release('x0',x0,'y0',y0,'sigma0',sigma0,'t0',t0,'t1',t1,...
				  'Ninternal',4,'tracers',{'salt','temp'});
							% _rel_ = the setup of a particle release
rel.verbose = 1;

steps = par_integrate(rel,run);
							% _steps_ = the actual Lagrangian integration
							
P = par_concatSteps(steps); % concatenates _steps_ variable by variable into
							% a single structure P that's easier to analyze