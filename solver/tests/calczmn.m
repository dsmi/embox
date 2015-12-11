function [ Z Ve Ye Vm Ym ] = calczmn(wg, ti, tj, tl, ttype, si, sj, sl, stype)
% Z=calczmn(wg, ti, tj, tl, ttype, si, sj, sl, stype)
%  Calculates element of the generalized impedance matrix. The summation is
%  done directly - not using FFT - so this function works quite slow and is
%  only used for tests.
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

if stype == 0 || stype == 1 
    if ttype == 0 || ttype == 1 
	% Use the tlines calculator to find the admittance between source and
	% observation points in the equivalent transmission lines
	Ye = 1./reshape(calc_vi(tle, z(tl), tl, z(sl), sl), maxm, maxn);
	Ym = 1./reshape(calc_vi(tlm, z(tl), tl, z(sl), sl), maxm, maxn);

	% Now find mode voltages from the mode currents
	Ve = Ie./Ye;
	Vm = Im./Ym;
    else
	% The case when the testing function is via - we need to
	% find the integral of the current over the target layer
	IIm = Im.*reshape(calc_iii(tlm, tl, 0, 1, z(sl), sl), maxm, maxn);
    end
else
    if ttype == 0 || ttype == 1 
	% Shift the via-induced voltage to the observation segment position
	Vm = Vdm.*reshape(calc_vvd(tlm, z(tl), tl, sl), maxm, maxn);
	Ve = Vm*0;
    else
	% Integral of current over the observation segment due to the
	% via-induced voltage
	iivd=calc_iivd(tlm, tl, sl);
	IIm = Vdm.*reshape(iivd, maxm, maxn);
	% Via self-impedance needs to be handled separately
	if tl == sl
	    % In the case of distributed voltage source, the transmission
	    % line equations are
	    %  d2I/dz2 - YZI = -YS
	    %  d2V/dz2 - YZV = 0
            % Where
            %  Z is the series impedance per len
            %  Y is the shunt admittance per len
            %  S is the voltage source per len
            % And the solutions are
            %  V(z)=Vp*exp(-gamma*z) + Vm*exp(gamma*z)
            %  I(z)=Ip*exp(-gamma*z) + Im*exp(gamma*z)+S/Z
	    % where
	    %  gamma = sqrt(YZ) is a propagation constant
	    %  Vp/Ip = -Vm/Im = sqrt(Z/Y) = Z0 is a characteristic impedance
	    % But - notice that the solution for I has the constant term S/Z.
	    % We need to subtract it from the current when evaluating
	    % the self-impedance!
	    IIm = IIm - Vdm.*h(sl).*Y0m(:,:,sl)./gamma(:,:,sl);
	end
    end
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
