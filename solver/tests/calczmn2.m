function [ Z ] = calczmn2(wg, ti, tj, tl, ttype, si, sj, sl, stype)
% Z=calczmn2(wg, ti, tj, tl, ttype, si, sj, sl, stype)
%  Calculates element of the generalized impedance matrix. The summation is
%  done directly - not using FFT - so this function works quite slow and is
%  only used for tests.
%  The difference between calczmn and calczmn2 is how the summation of the
%  waveguide modes is organized - calczmn2 is a step towards the fft
%  accelerated matrix elements evaluation.
%
% wg     - shiedling parameters, see wgparams
% ti, tj - testing segment coordinates in mesh
% tl     - testing segment position in the layers stackup
% ttype  - 1 for x-directed testing segment, 0 for y-directed, 2 for via
% si, sj - source segment coordinates in mesh
% sl     - source segment layer
% stype  - 1 for x-directed source segment, 0 for y-directed, 2 for via
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
h=wg.h; % thickness of the layers

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

% layers/tlines endpoints coordinates
z=cumsum(h); 

% Prepare inputs for the tlines calculators - one for te and one for
% tm modes. The reshaping is needed because the calculator only allows
% one dimension for the tline parameters.
ztlc = reshape(repmat(shiftdim([ 0 ; z(:) ], -2), maxm, maxn), [], nl+1);
ktlc = reshape(gamma, [], nl);
tle=calc_tlines(ztlc, reshape(1./Y0e, [], nl), ktlc, wg.Gls0, wg.Ggr0);
tlm=calc_tlines(ztlc, reshape(1./Y0m, [], nl), ktlc, wg.Gls0, wg.Ggr0);

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

% coordinates of the source segment center
xs=si*dx;
ys=sj*dy;
if stype == 1      % x-directed
   ys=ys+dy/2;
elseif stype == 0  % y-directed
   xs=xs+dx/2;
else               % via
   xs=xs+dx/2;
   ys=ys+dy/2;
end

% coordinates of the testing segment center
xt=ti*dx;
yt=tj*dy;
if ttype == 1      % x-directed
   yt=yt+dy/2;
elseif ttype == 0  % y-directed
   xt=xt+dx/2;
else               % via
   xt=xt+dx/2;
   yt=yt+dy/2;
end

% Use the tlines calculator to find the admittance between source and
% observation points in the equivalent transmission lines
Ye = 1./reshape(calc_vi(tle, z(tl), tl, z(sl), sl), maxm, maxn);
Ym = 1./reshape(calc_vi(tlm, z(tl), tl, z(sl), sl), maxm, maxn);

% x-directed source
if stype == 1
    if ttype == 1
	% x-directed testing function
	Gxx=Gdx_tri.*Gdx_tri.*Gdy_flat.*Gdy_flat.*(-Ne.*Ne.*ky.*ky./Ye-Nm.*Nm.*kx.*kx./Ym);
	%pxx=cos(kx.*xs).*sin(ky.*ys).*cos(kx.*xt).*sin(ky.*yt);
	pxx=(cos(kx.*(xt-xs)).*cos(ky.*(yt-ys))-cos(kx.*(xt-xs)).*cos(ky.*(yt+ys)) ...
	    +cos(kx.*(xt+xs)).*cos(ky.*(yt-ys))-cos(kx.*(xt+xs)).*cos(ky.*(yt+ys)))./4;
	Z=sum(sum(Gxx.*pxx, 2), 1);
    elseif ttype == 0
	% y-directed testing function
	Gyx=Gdx_tri.*Gdy_flat.*Gdx_flat.*Gdy_tri.*(Ne.*ky.*Ne.*kx./Ye-Nm.*kx.*Nm.*ky./Ym);
	%pyx=cos(kx.*xs).*sin(ky.*ys).*sin(kx.*xt).*cos(ky.*yt);
	pyx=(sin(kx.*(xt+xs)).*sin(ky.*(yt+ys))-sin(kx.*(xt+xs)).*sin(ky.*(yt-ys)) ...
             +sin(kx.*(xt-xs)).*sin(ky.*(yt+ys))-sin(kx.*(xt-xs)).*sin(ky.*(yt-ys)))./4;
	Z=sum(sum(Gyx.*pyx, 2), 1);
    else
	% via testing function
	iii = reshape(calc_iii(tlm, tl, 0, 1, z(sl), sl), maxm, maxn);
	m=kc.*kc./(j*freq*weps(tl));
	Gvx = Nm.*kx.*Gdx_tri.*Gdy_flat.*Nm.*Gdx_flat.*Gdy_flat.*m.*iii;
	%pvx = cos(kx.*xs).*sin(ky.*ys).*sin(kx.*xt).*sin(ky.*yt);
	pvx = (sin(kx.*(xt+xs)).*cos(ky.*(yt-ys)) - sin(kx.*(xt+xs)).*cos(ky.*(yt+ys)) ...
              +sin(kx.*(xt-xs)).*cos(ky.*(yt-ys)) - sin(kx.*(xt-xs)).*cos(ky.*(yt+ys)))./4;
	Z = sum(sum(Gvx.*pvx, 2), 1);
    end
elseif stype == 0
    if ttype == 1
	% x-directed testing function
	Gxy= Gdx_flat.*Gdy_tri.*Gdx_tri.*Gdy_flat.*(Ne.*kx.*Ne.*ky./Ye-Nm.*ky.*Nm.*kx./Ym);
	%pxy=sin(kx.*xs).*cos(ky.*ys).*cos(kx.*xt).*sin(ky.*yt);
	pxy=(-sin(kx.*(xt-xs)).*sin(ky.*(yt+ys))-sin(kx.*(xt-xs)).*sin(ky.*(yt-ys)) ...
	     +sin(kx.*(xt+xs)).*sin(ky.*(yt+ys))+sin(kx.*(xt+xs)).*sin(ky.*(yt-ys)))./4;
	Z=sum(sum( Gxy.*pxy, 2), 1);
    elseif ttype == 0
	% y-directed testing function
	Gyy=Gdx_flat.*Gdy_tri.*Gdx_flat.*Gdy_tri.*(-Ne.*kx.*Ne.*kx./Ye-Nm.*ky.*Nm.*ky./Ym);
	%pyy=sin(kx.*xs).*cos(ky.*ys).*sin(kx.*xt).*cos(ky.*yt);
	pyy=(cos(kx.*(xt-xs)).*cos(ky.*(yt-ys))+cos(kx.*(xt-xs)).*cos(ky.*(yt+ys)) ...
	    -cos(kx.*(xt+xs)).*cos(ky.*(yt-ys))-cos(kx.*(xt+xs)).*cos(ky.*(yt+ys)))./4;
	Z=sum(sum(Gyy.*pyy, 2), 1);
    else
	% via testing function
	iii = reshape(calc_iii(tlm, tl, 0, 1, z(sl), sl), maxm, maxn);
	m = kc.*kc./(j*freq*weps(tl));
	Gvy = Nm.*Gdx_flat.*Gdy_flat.*Nm.*ky.*Gdx_flat.*Gdy_tri.*m.*iii;
	%pvy = sin(kx.*xs).*cos(ky.*ys).*sin(kx.*xt).*sin(ky.*yt);
	pvy = (cos(kx.*(xt-xs)).*sin(ky.*(yt+ys)) + cos(kx.*(xt-xs)).*sin(ky.*(yt-ys)) ...
	      -cos(kx.*(xt+xs)).*sin(ky.*(yt+ys)) - cos(kx.*(xt+xs)).*sin(ky.*(yt-ys)))./4;
	Z = sum(sum(Gvy.*pvy, 2), 1);
    end
else
    % via source
    if ttype == 1
	% x-directed testing function
	m = kc.*kc./(j*freq*weps(sl));
	vvd = reshape(calc_vvd(tlm, z(tl), tl, sl), maxm, maxn);
	Gxv = -m.*Nm.*Gdx_flat.*Gdy_flat.*Nm.*kx.*Gdx_tri.*Gdy_flat.*vvd;
	%pxv = cos(kx.*xt).*sin(ky.*yt).*sin(kx.*xs).*sin(ky.*ys);
	pxv = (sin(kx.*(xt+xs)).*cos(ky.*(yt-ys))-sin(kx.*(xt+xs)).*cos(ky.*(yt+ys))...
	      -sin(kx.*(xt-xs)).*cos(ky.*(yt-ys))+sin(kx.*(xt-xs)).*cos(ky.*(yt+ys)))./4;
	Z=sum(sum(Gxv.*pxv, 2), 1);
    elseif ttype == 0
	% y-directed testing function
	m = kc.*kc./(j*freq*weps(sl));
	vvd = reshape(calc_vvd(tlm, z(tl), tl, sl), maxm, maxn);
	Gyv = -Nm.*ky.*Gdx_flat.*Gdy_tri.*Nm.*Gdx_flat.*Gdy_flat.*m.*vvd;
	%pyv = sin(kx.*xs).*sin(ky.*ys).*sin(kx.*xt).*cos(ky.*yt);
	pyv = (cos(kx.*(xt-xs)).*sin(ky.*(yt+ys)) - cos(kx.*(xt-xs)).*sin(ky.*(yt-ys))...
	      -cos(kx.*(xt+xs)).*sin(ky.*(yt+ys)) + cos(kx.*(xt+xs)).*sin(ky.*(yt-ys)))./4;
	Z = sum(sum(Gyv.*pyv, 2), 1);
    else
	% via testing function
	iivd = calc_iivd(tlm, tl, sl);
	r = reshape(iivd, maxm, maxn);
	if tl == sl
	   r = r - h(sl).*Y0m(:,:,sl)./gamma(:,:,sl);
	end
	m = -kc.^4/(freq*weps(sl)*freq*weps(tl));
	Gvv = Nm.*Gdx_flat.*Gdy_flat.*Nm.*Gdx_flat.*Gdy_flat.*m.*r;
	%pvv = sin(kx.*xt).*sin(ky.*yt).*sin(kx.*xs).*sin(ky.*ys);
	pvv = (cos(kx.*(xt-xs)).*cos(ky.*(yt-ys))-cos(kx.*(xt-xs)).*cos(ky.*(yt+ys))...
	      -cos(kx.*(xt+xs)).*cos(ky.*(yt-ys))+cos(kx.*(xt+xs)).*cos(ky.*(yt+ys)))./4;
	Z = sum(sum(Gvv.*pvv, 2), 1);
    end
end
