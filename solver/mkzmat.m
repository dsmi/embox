function Z=mkzmat(wg, mesh)
% Z=mkzmat(wg, xx, xy, yx, yy, dx, dy)
% Populates the impedance/reactions matrix
%
% wg     - shiedling parameters, see wgparams
% mesh   - meshed metal, see mkmesh
% 
 
nx=length(mesh.xi);
ny=length(mesh.yi);
n=nx+ny;

Z=zeros(n,n);

% Upper limits of the waveguide mode orders to be used when evaluating
% the matrix elements: T(E/M)[0..M-1][0..N-1] modes are used. 
maxm=wg.nx*wg.cnx;
maxn=wg.ny*wg.cny;

% angular frequency
freq=wg.freq;

a=wg.a; % x-size of the waveguide
b=wg.b; % y-size of the waveguide
h=wg.h; % height of the metallization above ground
c=wg.c; % upper ground height

% Layers stackup, from bottom to top along z
weps = wg.weps;
wmu  = wg.wmu;

% wavenumbers of the layers
k=freq*sqrt(weps.*wmu);

% wm(m,n)=m-1, wn(m,n)=n-1
[ wm, wn ] = ndgrid(0:maxm-1, 0:maxn-1);

% x and y wavenumbers of the waveguide mode
kx = wm*pi./a;
ky = wn*pi./b;

% cutoff wavenumber of the waveguide geometry
kc=sqrt(kx.^2+ky.^2);

% number of layers/sections
nl=length(weps);

% temporary, used below
k2=repmat(reshape(k.^2,1,1,nl),maxm,maxn);
kc2=repmat(kc.^2, [ 1 1 nl ]);

% waveguide layers/sections wavenumbers 
kz=sqrt(k2-kc2);

% Non-propagating wavenumbers for the loss-free case
kz2=sqrt(-k2+kc2)./j;
kz(find(real(k2)<kc2))=kz2(find(real(k2)<kc2));

% waveguide layers/sections propagation constants
gamma=j*kz;

% Characteristic admittances of the TE modes
Y0e=gamma./(j*freq*repmat(shiftdim(wmu(:), -2), maxm, maxn));

% Characteristic admittances of the TM modes
Y0m=(j*freq*repmat(shiftdim(weps(:), -2), maxm, maxn))./gamma;

% normalization coefficients for te and tm waveguide modes
[ Ne, Nm ] = wnorm(a, b, maxm, maxn);

% left-looking reflection coefficient at z=0
Gls0=wg.Gls0;
% left-looking reflection coefficient at z=h
Gls=Gls0*exp(-2*gamma(:,:,1)*h);

% right-looking reflection coefficient at z=c
Ggr0=wg.Ggr0; % metal at top
% right-looking reflection coefficient at z=h
Ggr=Ggr0*exp(2*gamma(:,:,2)*(h-c));

% TE and TM admittances at z->h+
Ye_gr=(1-Ggr)./(1+Ggr).*Y0e(:,:,2);
Ym_gr=(1-Ggr)./(1+Ggr).*Y0m(:,:,2);

% TE and TM admittances at z->h-
Ye_ls=-(1-Gls)./(1+Gls).*Y0e(:,:,1);
Ym_ls=-(1-Gls)./(1+Gls).*Y0m(:,:,1);

% mesh cell sizes
dx=wg.a/wg.nx;
dy=wg.b/wg.ny;

% Multiplier resulting from x-integration of the triangular basis function
Gdx_tri=gtri(dx,kx);

% Multiplier resulting from x-integration of the constant (rectangular) b.f.
Gdx_flat=gflat(dx,kx);

% Multiplier resulting from y-integration of the triangular b.f.
Gdy_tri=gtri(dy,ky);

% Multiplier resulting from y-integration of the constant (rectangular) b.f.
Gdy_flat=gflat(dy,ky);

% Admittance to be used when finding the mode voltages
Ye=Ye_gr-Ye_ls;
Ym=Ym_gr-Ym_ls;

% x-directed testing, x-directed source
Gxx=Gdx_tri.*Gdx_tri.*Gdy_flat.*Gdy_flat.*(-Ne.*Ne.*ky.*ky./Ye-Nm.*Nm.*kx.*kx./Ym);

[ cc, ss ] = myfft(Gxx);

[ mxi, nxi ] = ndgrid(mesh.xi, mesh.xi);
[ mxj, nxj ] = ndgrid(mesh.xj, mesh.xj);

idif = wrapidx((mxi-nxi)*wg.cnx/2+1, size(Gxx, 1));   % xt-xs
isum = wrapidx((mxi+nxi)*wg.cnx/2+1, size(Gxx, 1));   % xt+xs
jdif = wrapidx((mxj-nxj)*wg.cny/2+1, size(Gxx, 2));   % yt-ys
jsum = wrapidx((mxj+nxj+1)*wg.cny/2+1, size(Gxx, 2)); % yt+ys
% notice one is added to mxj+nxj when computing jsum - y-center of
% the x-directed basis is found as j*dy+dy/2, (j is the mesh index)
% the two halves added together give one

idif_jdif = sub2ind(size(cc), idif, jdif);
idif_jsum = sub2ind(size(cc), idif, jsum);
isum_jdif = sub2ind(size(cc), isum, jdif);
isum_jsum = sub2ind(size(cc), isum, jsum);

Zxx = (cc(idif_jdif) - cc(idif_jsum) + cc(isum_jdif) - cc(isum_jsum)) ./ 4;

% y-directed testing, x-directed source
Gyx=Gdx_tri.*Gdy_flat.*Gdx_flat.*Gdy_tri.*(Ne.*ky.*Ne.*kx./Ye-Nm.*kx.*Nm.*ky./Ym);

[ cc, ss ] = myfft(Gyx);

[ myi, nxi ] = ndgrid(mesh.yi, mesh.xi);
[ myj, nxj ] = ndgrid(mesh.yj, mesh.xj);

% y-position of x-directed bases (and x-position of y-directed) is the
% bottom/left edge (one with minimal x/y) - that is the reason for adding
% or subtracting wg.cnx/4 or wg.cny/4 (which corresponds to half-cell)
% Here for Zyx:
%  xt=myi*dx+dx/2 (y-directed testing)
%  ys=nxj*dy+dy/2 (x-directed source)
idif = wrapidx((myi-nxi)*wg.cnx/2 + wg.cnx/4 + 1, size(Gyx, 1)); % xt-xs
isum = wrapidx((myi+nxi)*wg.cnx/2 + wg.cnx/4 + 1, size(Gyx, 1)); % xt+xs
jdif = wrapidx((myj-nxj)*wg.cny/2 - wg.cny/4 + 1, size(Gyx, 2)); % yt-ys
jsum = wrapidx((myj+nxj)*wg.cny/2 + wg.cny/4 + 1, size(Gyx, 2)); % yt+ys

idif_jdif = sub2ind(size(ss), idif, jdif);
idif_jsum = sub2ind(size(ss), idif, jsum);
isum_jdif = sub2ind(size(ss), isum, jdif);
isum_jsum = sub2ind(size(ss), isum, jsum);

Zyx = (ss(isum_jsum) - ss(isum_jdif) + ss(idif_jsum) - ss(idif_jdif)) ./ 4;

% x-directed testing, y-directed source
Gxy=Gdx_flat.*Gdy_tri.*Gdx_tri.*Gdy_flat.*(Ne.*kx.*Ne.*ky./Ye-Nm.*ky.*Nm.*kx./Ym);

[ cc, ss ] = myfft(Gxy);

[ mxi, nyi ] = ndgrid(mesh.xi, mesh.yi);
[ mxj, nyj ] = ndgrid(mesh.xj, mesh.yj);

% y-position of x-directed bases (and x-position of y-directed) is the
% bottom/left edge (one with minimal x/y) - that is the reason for adding
% or subtracting wg.cnx/4 or wg.cny/4 (which corresponds to half-cell)
% Here for Zyx:
%  yt=mxj*dy+dy/2 (x-directed testing)
%  xs=nyi*dx+dx/2 (y-directed source)
idif = wrapidx((mxi-nyi)*wg.cnx/2 - wg.cnx/4 + 1, size(Gyx, 1)); % xt-xs
isum = wrapidx((mxi+nyi)*wg.cnx/2 + wg.cnx/4 + 1, size(Gyx, 1)); % xt+xs
jdif = wrapidx((mxj-nyj)*wg.cny/2 + wg.cny/4 + 1, size(Gyx, 2)); % yt-ys
jsum = wrapidx((mxj+nyj)*wg.cny/2 + wg.cny/4 + 1, size(Gyx, 2)); % yt+ys

idif_jdif = sub2ind(size(ss), idif, jdif);
idif_jsum = sub2ind(size(ss), idif, jsum);
isum_jdif = sub2ind(size(ss), isum, jdif);
isum_jsum = sub2ind(size(ss), isum, jsum);

Zxy = (-ss(idif_jsum) - ss(idif_jdif) + ss(isum_jsum) + ss(isum_jdif)) ./ 4;

% y-directed testing, y-directed source
Gyy=Gdx_flat.*Gdy_tri.*Gdx_flat.*Gdy_tri.*(-Ne.*kx.*Ne.*kx./Ye-Nm.*ky.*Nm.*ky./Ym);

[ cc, ss ] = myfft(Gyy);

[ myi, nyi ] = ndgrid(mesh.yi, mesh.yi);
[ myj, nyj ] = ndgrid(mesh.yj, mesh.yj);

% notice one is added to myi+nyi when computing isum - x-center of
% the y-directed basis is found as i*dx+dx/2, (i is the mesh index)
% the two halves added together give one
idif = wrapidx((myi-nyi)*wg.cnx/2+1, size(Gyy, 1));   % xt-xs
isum = wrapidx((myi+nyi+1)*wg.cnx/2+1, size(Gyy, 1)); % xt+xs
jdif = wrapidx((myj-nyj)*wg.cny/2+1, size(Gyy, 2));   % yt-ys
jsum = wrapidx((myj+nyj)*wg.cny/2+1, size(Gyy, 2));   % yt+ys

idif_jdif = sub2ind(size(cc), idif, jdif);
idif_jsum = sub2ind(size(cc), idif, jsum);
isum_jdif = sub2ind(size(cc), isum, jdif);
isum_jsum = sub2ind(size(cc), isum, jsum);

Zyy = (cc(idif_jdif) + cc(idif_jsum) - cc(isum_jdif) - cc(isum_jsum)) ./ 4;

% compose the entire matrix
Z = [ Zxx Zxy ; Zyx Zyy ];

% Identify segments which cross the waveguide boundary - the corresponding
% elements of the Z matrix need to be multiplied by 0.5
xmul=-(~mesh.xi | ~(mesh.xi-wg.nx))*0.5+1.0;
ymul=-(~mesh.yj | ~(mesh.yj-wg.ny))*0.5+1.0;
M=diag([ xmul(:) ; ymul(:) ]);

% scale the impedance matrix elements for the boundary-crossing segments
Z=M*Z*M;
