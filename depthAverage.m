function cavg = depthAverage(c,z,zr);
 
% cavg = depthAverage(c,z,zr);
%
% integrates c(z) from zr(1) to zr(2).
% c and z are both size [I K], where K is the z dimension. There's no assumption 
% that z is uniform in the K direction: each column is handled separately.

zr = sort(zr);
[I,K] = size(c);

% if zr extends beyond any of the individual column tops or bottoms, move it to
% the column limits
zmin = repmat(zr(1),[I 1]);
zmin = max(zmin,min(z,[],2));
zmax = repmat(zr(2),[I 1]);
zmax = min(zmax,max(z,[],2));

% interpolate c onto zmin and zmax, and sort to include these points in the
% c(z) mesh
ii = (1:I)';
iii = repmat(ii,[1 K]);
c_zmin = griddata(iii,z,c,ii,zmin);
c_zmax = griddata(iii,z,c,ii,zmax);
z = cat(2,z,zmin,zmax);
c = cat(2,c,c_zmin,c_zmax);
[z,sorti] = sort(z,2);
c = c(sub2ind(size(c),repmat(ii,[1 K+2]),sorti));

% points in between, and matching cell sizes
c2 = 0.5.*(c(:,1:end-1) + c(:,2:end));
z2 = 0.5.*(z(:,1:end-1) + z(:,2:end));
dz = diff(z,[],2);

% find the ones between zmin and zmax, and sum
isin = (z2 >= repmat(zmin,[1 K+1]) & z2 <= repmat(zmax,[1 K+1]));
cint = sum(isin.*c2.*dz,2);
DZ = sum(isin.*dz,2);

% divide by range of integration DZ to get the average; but if DZ=0,
% interpret this as a point interpolation, and see if omitting the dz
% weighting gives valid values
cavg = cint./DZ;
f = find(DZ==0);
cavg(f) = sum(isin(f,:).*c2(f,:),2) ./ sum(isin(f,:),2);