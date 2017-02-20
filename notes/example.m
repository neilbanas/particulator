run = modelRun_romsCascadia('dirname/');
							% _run_ identifies a model run to use as source data

x0 = run.grid.lon;			% pick initial coordinates
y0 = run.grid.lat;
sigma0 = zeros(size(x0));
t0 = run.t(1) .* ones(size(x0));
t1 = run.t(end);
rel = par_release('x0',x0,'y0',y0,'sigma0',sigma0,'t0',t0,'t1',t1,...
				  'dt_per_DT',1,'tracers',{'salt','temp'});
							% _rel_ = the setup of a particle release

steps = par_integrate(rel,run);
							% _steps_ = the actual Lagrangian integration
							
P = par_concatSteps(steps); % concatenates _steps_ variable by variable into
							% a single structure P that's easier to analyze