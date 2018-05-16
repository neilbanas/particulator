% defining a land mask on the u grid ("grid.dat.pop") empirically, as points
% where the annual max of the surface flow is 0 (in 2009)
%
% Neil Banas, May 2018 

modelDir = '/Users/neil/Dropbox/CAO/datasets/biomas-physics-2009/';
gridDir = '/Users/neil/Dropbox/particulator/data/biomas/';
depthRange = [-150 0];
run = modelRun_biomas2d(modelDir,2009,gridDir,depthRange);

nn = 1:length(run.t);
t = run.t(nn);
N = length(t);
for ni=1:length(nn)
	disp(ni);
	n = nn(ni);
	run.loadFrame(n,{});
	usurf(:,:,ni) = run.F1.u(:,:,1);
	vsurf(:,:,ni) = run.F1.v(:,:,1);
end
umag = sqrt(usurf.^2 + vsurf.^2);
umax = max(umag,[],3);

xu = run.grid.xu;
yu = run.grid.yu;
masku = double(umax>0);
% there are 2 rho points in Franz Josef Land where interpolating from this 
% velocity-based u mask onto the rho mask, one would think there was flow but 
% really there's land
% meanwhile, there are 1451 rho points on the coastal margins where the rho mask
% says there is water but the interpolated flow is 0
% results are similar interpolating the other way, from rho points to u points
save velocityBasedMask xu yu masku
