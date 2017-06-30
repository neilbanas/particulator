% load a biomas model run
%modelDir = '/Volumes/sambal/biomas/2009/';
modelDir = '/Volumes/littlegray/biomas2009/';
gridDir = '/Users/neil/Dropbox/cmg_tools/particulator/data/biomas/';
depthRange = [0 0];
run = modelRun_biomas2d(modelDir,2009,gridDir,depthRange);

quick = 1;
% make a blank return map in the Bering-Chukchi region
if quick
	x00 = 120 : 0.25 : 240;
	y00 = 53 : 0.25 : 75;
	t0 = run.t(1) + (0:3:12); % short litle thing
else
	x00 = 120 : 0.25 : 240;
	y00 = 53 : 0.25 : 75;
	t0 = run.t(1) + (0:10:360); % a full year
end
[x0,y0] = meshgrid(x00,y00);
m = run.interpMask(x0,y0) > 0.5;
x0 = x0(m);
y0 = y0(m);
map = returnMap(x0,y0,t0);

% fill in the return map using the model run
map.generate(run,'KH',0,'Nreplicates',10,...
			 'Ninternal',2,'verbose',1,'diffusive',0,...
			 'tracers',{'ice','temp','swrad','diatom','flagel'});
save example_biomas_returnMap map
		
% use the return map to track particles	 
P = map.integrate();
