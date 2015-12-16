function [ Z Ve Ye Vm Ym ] = calczmn(wg, ti, tj, tl, ttype, si, sj, sl, stype)
% Z=calczmn(wg, ti, tj, tl, ttype, si, sj, sl, stype)
%  Calculates element of the generalized impedance matrix. The summation is
%  done directly - not using FFT - so this function works quite slow and is
%  only used for tests.
%
% wg     - shiedling parameters, see wgparams
% ti, tj - testing segment coordinates in mesh
% tl     - testing segment position in the layers stackup
% ttype  - 1 for x-directed testing segment, 0 for y-directed, 2-4 for via
% si, sj - source segment coordinates in mesh
% sl     - source segment layer
% stype  - 1 for x-directed source segment, 0 for y-directed, 2-4 for via
%
% Via types are:
%  2 - via to the previous layer, 3 - via to the next layer, 4 - 'through'
%      via from the previous to the next layer
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

% The waveguide modes are described in terms of equivalent transmission
% lines, and here we compute the current sources in the transmission lines
% which corresponds to the source current in the waveguide.
if stype == 1
    % TE and TM mode current sources, x-directed current
    Ie=-Ne.*ky.*Gdx_tri.*Gdy_flat.*cos(kx.*xs).*sin(ky.*ys);
    Im= Nm.*kx.*Gdx_tri.*Gdy_flat.*cos(kx.*xs).*sin(ky.*ys);
elseif stype == 0
    % TE and TM mode current sources, y-directed current
    Ie=Ne.*kx.*Gdx_flat.*Gdy_tri.*sin(kx.*xs).*cos(ky.*ys);
    Im=Nm.*ky.*Gdx_flat.*Gdy_tri.*sin(kx.*xs).*cos(ky.*ys);
else               
    % Calculate the distributed voltage source intensity in the waveguide
    % equivalent transmission line due to the via current by evaluating
    % integral of the via current multiplied by the eigenfunction.
    % The eigenfunction for the tm modes is: Psi_m=nm*sin(kx*x)*sin(ky*y)
    m=kc.*kc./(j*freq*weps(sl));
    Vdm=m.*Nm.*Gdx_flat.*Gdy_flat.*sin(kx.*xs).*sin(ky.*ys);
end

% Only used when both source and testing are horizontal, otherwise zero
Ye = 0;
Ym = 0;
Ve = 0;
Vm = 0;

% horizontal-to-horizontal
if (stype == 0 || stype == 1) && (ttype == 0 || ttype == 1)

    % Use the tlines calculator to find the admittance between source and
    % observation points in the equivalent transmission lines
    ze = reshape(calc_vi(tle, z(tl), tl, z(sl), sl), maxm, maxn);
    zm = reshape(calc_vi(tlm, z(tl), tl, z(sl), sl), maxm, maxn);
    Ye = 1./ze;
    Ym = 1./zm;

    % Now find mode voltages from the mode currents
    Ve = Ie.*ze;
    Vm = Im.*zm;

% horizontal-to-via
elseif (stype == 0 || stype == 1) && (ttype >= 2)

    % Current src due to the horz current -> integral of current over the via
    iii = 0;

    % Connects to the previous layer - increasing: i = (z-zi)/d
    if (ttype == 2 || ttype == 4)
        iii = iii + calc_iii(tlm, tl  , 1./tlm.d(tl),   0, z(sl), sl);
    end

    % Connects to the next layer - decreasing: i = 1 - (z-z{i+i})/d{i+1}
    if (ttype == 3 || ttype == 4)
        iii = iii + calc_iii(tlm, tl+1, -1./tlm.d(tl+1), 1, z(sl), sl);
    end

    % Integral of the current over the via layer(s)
    IIm = Im.*reshape(iii, maxm, maxn);

% via-to-horizontal
elseif (stype >= 2) && (ttype == 0 || ttype == 1)

    % Modal V at the observation due to the via distributed source
    vvd = 0;

    % Connects to the previous layer - increasing (linear only)
    if (stype == 2 || stype == 4)
        vvd = vvd + calc_vvl(tlm, z(tl), tl, sl, 1./tlm.d(tl));
    end

    % Connects to the next layer - decreasing (both linear and constant parts)
    if (stype == 3 || stype == 4)
        vvd = vvd + calc_vvd(tlm, z(tl), tl, sl+1);
        vvd = vvd + calc_vvl(tlm, z(tl), tl, sl+1, -1./tlm.d(tl+1));
    end

    Vm = Vdm.*reshape(vvd, maxm, maxn);
    Ve = Vm*0;

else
% via-to-via

    % Integral of current over the observation segment due to the
    % via-induced voltage. Notice that we drop non-exponential terms from
    % the current integrals (by passing c=0)
    iivd = 0;

    % testing - prev, source - prev
    if (ttype == 2 || ttype == 4) && (stype == 2 || stype == 4)
        iivd = iivd + calc_iivl(tlm, tl, 1./tlm.d(tl), 0, -1, sl, 1./tlm.d(sl));
    end

    % testing - prev, source - next
    if (ttype == 2 || ttype == 4) && (stype == 3 || stype == 4)
        iivd = iivd + calc_iivd(tlm, tl, 1./tlm.d(tl), 0, -1, sl+1);
        iivd = iivd + calc_iivl(tlm, tl, 1./tlm.d(tl), 0, -1, sl+1, -1./tlm.d(sl+1));
    end

    % testing - next, source - prev
    if (ttype == 3 || ttype == 4) && (stype == 2 || stype == 4)
        iivd = iivd + calc_iivl(tlm, tl+1, -1./tlm.d(tl+1), 1, -1, sl, 1./tlm.d(sl));
    end

    % testing - next, source - next
    if (ttype == 3 || ttype == 4) && (stype == 3 || stype == 4)
        iivd = iivd + calc_iivd(tlm, tl+1, -1./tlm.d(tl+1), 1, -1, sl+1);
        iivd = iivd + calc_iivl(tlm, tl+1, -1./tlm.d(tl+1), 1, -1, sl+1, -1./tlm.d(sl+1));
    end

    IIm = Vdm.*reshape(iivd, maxm, maxn);

end

if ttype == 1
    % x-directed testing function
    Ze=sum(sum( Ve.*Ne.*ky.*Gdx_tri.*Gdy_flat.*cos(kx.*xt).*sin(ky.*yt), 2), 1);
    Zm=sum(sum(-Vm.*Nm.*kx.*Gdx_tri.*Gdy_flat.*cos(kx.*xt).*sin(ky.*yt), 2), 1);
elseif ttype == 0
    % y-directed testing function
    Ze=sum(sum(-Ve.*Ne.*kx.*Gdx_flat.*Gdy_tri.*sin(kx.*xt).*cos(ky.*yt), 2), 1);
    Zm=sum(sum(-Vm.*Nm.*ky.*Gdx_flat.*Gdy_tri.*sin(kx.*xt).*cos(ky.*yt), 2), 1);
else               
    % via
    m = kc.*kc./(j*freq*weps(tl));
    Zm = sum(sum(m.*IIm.*Nm.*Gdx_flat.*Gdy_flat.*sin(kx.*xt).*sin(ky.*yt), 2), 1);
    Ze = 0;
end

Z=Ze+Zm;
