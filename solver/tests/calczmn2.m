function [ Z ] = calczmn2(wg, ti, tj, xdt, si, sj, xds)
% Z=calczmn2(wg, ti, tj, xdt, si, yj, xds)
%  Calculates element of the generalized impedance matrix. The summation is
%  done directly - not using FFT - so this function works quite slow and is
%  only used for tests.
%  The difference between calczmn and calczmn2 is how the summation of the
%  waveguide modes is organized - calczmn2 is a step towards the fft
%  accelerated matrix elements evaluation.
%
% wg     - shiedling parameters, see wgparams
% ti, tj - testing segment coordinates in mesh
% xdt    - 1 for x-directed testing segment, 0 for y-directed
% si, sj - source segment coordinates in mesh
% yds    - 1 for x-directed source segment, 0 for y-directed
%
% Returns the matrix element and mode voltages and admittances which can be
% (are supposed to be) used for some additional tests.
%

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

% coordinates of the source segment center
xs=si*dx;
ys=sj*dy;
if xds
   ys=ys+dy/2;
else
   xs=xs+dx/2;
end

% coordinates of the testing segment center
xt=ti*dx;
yt=tj*dy;
if xdt
   yt=yt+dy/2;
else
   xt=xt+dx/2;
end

% x-directed source
if xds
    if xdt
	% x-directed testing function
	Gxx=Gdx_tri.*Gdx_tri.*Gdy_flat.*Gdy_flat.*(-Ne.*Ne.*ky.*ky./Ye-Nm.*Nm.*kx.*kx./Ym);
	%pxx=cos(kx.*xs).*sin(ky.*ys).*cos(kx.*xt).*sin(ky.*yt);
	pxx=(cos(kx.*(xt-xs)).*cos(ky.*(yt-ys))-cos(kx.*(xt-xs)).*cos(ky.*(yt+ys)) ...
	    +cos(kx.*(xt+xs)).*cos(ky.*(yt-ys))-cos(kx.*(xt+xs)).*cos(ky.*(yt+ys)))./4;
	Z=sum(sum(Gxx.*pxx, 2), 1);
    else
	% y-directed testing function
	Gyx=Gdx_tri.*Gdy_flat.*Gdx_flat.*Gdy_tri.*(Ne.*ky.*Ne.*kx./Ye-Nm.*kx.*Nm.*ky./Ym);
	%pyx=cos(kx.*xs).*sin(ky.*ys).*sin(kx.*xt).*cos(ky.*yt);
	pyx=(sin(kx.*(xt+xs)).*sin(ky.*(yt+ys))-sin(kx.*(xt+xs)).*sin(ky.*(yt-ys))
             +sin(kx.*(xt-xs)).*sin(ky.*(yt+ys))-sin(kx.*(xt-xs)).*sin(ky.*(yt-ys)))./4;
	Z=sum(sum(Gyx.*pyx, 2), 1);
    end
else
    if xdt
	% x-directed testing function
	Gxy= Gdx_flat.*Gdy_tri.*Gdx_tri.*Gdy_flat.*(Ne.*kx.*Ne.*ky./Ye-Nm.*ky.*Nm.*kx./Ym);
	%pxy=sin(kx.*xs).*cos(ky.*ys).*cos(kx.*xt).*sin(ky.*yt);
	pxy=(-sin(kx.*(xt-xs)).*sin(ky.*(yt+ys))-sin(kx.*(xt-xs)).*sin(ky.*(yt-ys)) ...
	     +sin(kx.*(xt+xs)).*sin(ky.*(yt+ys))+sin(kx.*(xt+xs)).*sin(ky.*(yt-ys)))./4;
	Z=sum(sum( Gxy.*pxy, 2), 1);
    else
	% y-directed testing function
	Gyy=Gdx_flat.*Gdy_tri.*Gdx_flat.*Gdy_tri.*(-Ne.*kx.*Ne.*kx./Ye-Nm.*ky.*Nm.*ky./Ym);
	%pyy=sin(kx.*xs).*cos(ky.*ys).*sin(kx.*xt).*cos(ky.*yt);
	pyy=(cos(kx.*(xt-xs)).*cos(ky.*(yt-ys))+cos(kx.*(xt-xs)).*cos(ky.*(yt+ys)) ...
	    -cos(kx.*(xt+xs)).*cos(ky.*(yt-ys))-cos(kx.*(xt+xs)).*cos(ky.*(yt+ys)))./4;
	Z=sum(sum(Gyy.*pyy, 2), 1);
    end
end
