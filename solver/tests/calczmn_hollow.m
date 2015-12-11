function [ Z Ie Im Vdm tle tlm ] = calczmn_hollow(wg, ti, tj, tl, ttype, si, sj, sl, stype)
% [ Z Ie Im Vdm tle tlm ] = calczmn_hollow(wg, ti, tj, tl, ttype, si, sj, sl, stype)
%  Calculates element of the generalized impedance matrix. The summation is
%  done directly - not using FFT - so this function works quite slow and is
%  only used for tests.
%  Unlike calczmn this uses a hollow via with caps, result of some experiments.
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
z = [ 0 reshape(cumsum(h), 1, []) ]; 

% Prepare inputs for the tlines calculators - one for te and one for
% tm modes. The reshaping is needed because the calculator only allows
% one dimension for the tline parameters.
ztlc = reshape(repmat(shiftdim(z(:), -2), maxm, maxn), [], nl+1);
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

% Multipliers resulting from x- and y-integration of the 'out' basis function
% Notice that it also replaces sin with cos and cos with -sin! (see gout
% description for more details)
Gdx_out=gout(dx,kx);
Gdy_out=gout(dy,ky);

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

    % Via gives: current sources at the top of the layer (the same Ie and Im as
    % the horizontal segments), current sources at the bottom of the layer
    % (the same as Ie and Im but with negative sign) 
    % the distributed voltage source

    % Calculate the distributed voltage source intensity in the waveguide
    % equivalent transmission line due to the vertical via current by
    % evaluating integral of the via current multiplied by the eigenfunction.
    % The eigenfunction for the tm modes is: Psi_m=nm*sin(kx*x)*sin(ky*y)
    m = kc.*kc./(j*freq*weps(sl));
    % y-parallel edges (integration over x) and x-parallel edges (integration over y)
    % sin(ky.*(ys-dx/2)) + sin(ky.*(ys+dx/2)) = sin(ky.*ys).*2.*cos(ky.*dy./2)
    g = (Gdx_flat.*2.*cos(ky.*dy./2)+Gdy_flat.*2.*cos(kx.*dx./2));
    % Distributed voltage source due to vertical current in the source via
    Vdm = m.*Nm.*g.*sin(kx.*xs).*sin(ky.*ys);

    % TE and TM mode current sources due to the due to the bottom or top of
    % the via. For the top the signs need to be changed because the currents
    % are 'out' at the bottom and 'in' at the top
    Ie = Ne.*(ky.*Gdx_out.*Gdy_flat - kx.*Gdx_flat.*Gdy_out).*sin(kx.*xs).*sin(ky.*ys);
    Im = Nm.*(-kx.*Gdx_out.*Gdy_flat - ky.*Gdx_flat.*Gdy_out).*sin(kx.*xs).*sin(ky.*ys);

end

% horizontal-to-horizontal
if (stype == 0 || stype == 1) && (ttype == 0 || ttype == 1)

    % Use the tlines calculator to find the admittance between source and
    % observation points in the equivalent transmission lines
    Ye = 1./reshape(calc_vi(tle, z(tl+1), tl, z(sl+1), sl), maxm, maxn);
    Ym = 1./reshape(calc_vi(tlm, z(tl+1), tl, z(sl+1), sl), maxm, maxn);

    % Now find mode voltages from the mode currents
    Ve = Ie./Ye;
    Vm = Im./Ym;

% horizontal-to-via
elseif (stype == 0 || stype == 1) && (ttype == 2)

    % Integral of the current over the via layer
    IIm = Im.*reshape(calc_iii(tlm, tl, 0, 1, z(sl + 1), sl), maxm, maxn);

    % Mode voltages at the top of the via
    Ve_top = Ie.*reshape(calc_vi(tle, z(tl+1), tl, z(sl+1), sl), maxm, maxn);
    Vm_top = Im.*reshape(calc_vi(tlm, z(tl+1), tl, z(sl+1), sl), maxm, maxn);

    % Mode voltages at the bottom of the via
    Ve_bottom = Ie.*reshape(calc_vi(tle, z(tl), tl, z(sl+1), sl), maxm, maxn);
    Vm_bottom = Im.*reshape(calc_vi(tlm, z(tl), tl, z(sl+1), sl), maxm, maxn);

% via-to-horizontal
elseif (stype == 2) && (ttype == 0 || ttype == 1)

    % Modal V at the observation due to the via distributed source
    vvd = calc_vvd(tlm, z(tl+1), tl, sl);
    Vm = Vdm.*reshape(vvd, maxm, maxn);
    Ve = Vm*0;

    % Plus modal V at the observation due to the bottom horizontal segment.
    Ve = Ve + Ie.*reshape(calc_vi(tle, z(tl+1), tl, z(sl), sl), maxm, maxn);
    Vm = Vm + Im.*reshape(calc_vi(tlm, z(tl+1), tl, z(sl), sl), maxm, maxn);

    % Plus modal V at the observation due to the top horizontal segment of the
    % via. Notice minus signs - Ie and Im were calculated for the 'out' basis,
    % while the top is 'in'
    Ve = Ve - Ie.*reshape(calc_vi(tle, z(tl+1), tl, z(sl+1), sl), maxm, maxn);
    Vm = Vm - Im.*reshape(calc_vi(tlm, z(tl+1), tl, z(sl+1), sl), maxm, maxn);

else
% via-to-via!, we need to calculate V(e/m)_(top/bottom) and IIm

    % Integral of current over the observation segment due to the
    % via-induced voltage
    iivd = calc_iivd(tlm, tl, sl);
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

    % Integral over the vertical part of the observation via due to the bottom
    % part of the source via
    IIm = IIm + Im.*reshape(calc_iii(tlm, tl, 0, 1, z(sl), sl), maxm, maxn);

    % Integral over the vertical part of the observation via due to the top part
    % of the source - notice sign (top is 'in') and z(sl+1)
    IIm = IIm - Im.*reshape(calc_iii(tlm, tl, 0, 1, z(sl+1), sl), maxm, maxn);

    % source: vertical, observation: bottom
    Vm_bottom =           Vdm.*reshape(calc_vvd(tlm, z(tl), tl,        sl), maxm, maxn);
    Ve_bottom = Vm_bottom*0;

    % source: bottom, observation: bottom
    Ve_bottom = Ve_bottom + Ie.*reshape(calc_vi(tle, z(tl), tl, z(sl), sl), maxm, maxn);
    Vm_bottom = Vm_bottom + Im.*reshape(calc_vi(tlm, z(tl), tl, z(sl), sl), maxm, maxn);

    % source: top, observation: bottom
    % Notice minus signs - Ie and Im were calculated for the 'out' basis,
    % while the top is 'in'
    Ve_bottom = Ve_bottom - Ie.*reshape(calc_vi(tle, z(tl), tl, z(sl+1), sl), maxm, maxn);
    Vm_bottom = Vm_bottom - Im.*reshape(calc_vi(tlm, z(tl), tl, z(sl+1), sl), maxm, maxn);

    % source: vertical, observation: top
    Vm_top =        Vdm.*reshape(calc_vvd(tlm, z(tl+1), tl,        sl), maxm, maxn);
    Ve_top = Vm_top*0;

    % source: bottom, observation: top
    Ve_top = Ve_top + Ie.*reshape(calc_vi(tle, z(tl+1), tl, z(sl), sl), maxm, maxn);
    Vm_top = Vm_top + Im.*reshape(calc_vi(tlm, z(tl+1), tl, z(sl), sl), maxm, maxn);

    % source: top, observation: top
    % Notice minus signs - Ie and Im were calculated for the 'out' basis,
    % while the top is 'in'
    Ve_top = Ve_top - Ie.*reshape(calc_vi(tle, z(tl+1), tl, z(sl+1), sl), maxm, maxn);
    Vm_top = Vm_top - Im.*reshape(calc_vi(tlm, z(tl+1), tl, z(sl+1), sl), maxm, maxn);

end

if ttype == 1
    % x-directed testing function
    Ze = sum(sum( Ve.*Ne.*ky.*Gdx_tri.*Gdy_flat.*cos(kx.*xt).*sin(ky.*yt), 2), 1);
    Zm = sum(sum(-Vm.*Nm.*kx.*Gdx_tri.*Gdy_flat.*cos(kx.*xt).*sin(ky.*yt), 2), 1);
elseif ttype == 0
    % y-directed testing function
    Ze = sum(sum(-Ve.*Ne.*kx.*Gdx_flat.*Gdy_tri.*sin(kx.*xt).*cos(ky.*yt), 2), 1);
    Zm = sum(sum(-Vm.*Nm.*ky.*Gdx_flat.*Gdy_tri.*sin(kx.*xt).*cos(ky.*yt), 2), 1);
else

    % vertical segment of the via
    m = kc.*kc./(j*freq*weps(tl));
    % y-parallel edges (integration over x) and x-parallel edges (integration over y)
    % sin(ky.*(ys-dx/2)) + sin(ky.*(ys+dx/2)) = sin(ky.*ys).*2.*cos(ky.*dy./2)
    g = (Gdx_flat.*2.*cos(ky.*dy./2)+Gdy_flat.*2.*cos(kx.*dx./2));
    Zm = sum(sum(m.*IIm.*Nm.*g.*sin(kx.*xt).*sin(ky.*yt), 2), 1);
    Ze = Zm * 0;

    % bottom of the via
    me = -Ve_bottom.*Ne.*ky.*Gdx_out.*Gdy_flat + Ve_bottom.*Ne.*kx.*Gdx_flat.*Gdy_out;
    mm =  Vm_bottom.*Nm.*kx.*Gdx_out.*Gdy_flat + Vm_bottom.*Nm.*ky.*Gdx_flat.*Gdy_out;
    ss = sin(kx.*xt).*sin(ky.*yt);
    Ze = Ze + sum(sum(me.*ss, 2), 1);
    Zm = Zm + sum(sum(mm.*ss, 2), 1);

    % top of the via - the signs are changed because the top basis is 'in'
    me =  Ve_top.*Ne.*ky.*Gdx_out.*Gdy_flat - Ve_top.*Ne.*kx.*Gdx_flat.*Gdy_out;
    mm = -Vm_top.*Nm.*kx.*Gdx_out.*Gdy_flat - Vm_top.*Nm.*ky.*Gdx_flat.*Gdy_out;
    ss = sin(kx.*xt).*sin(ky.*yt);
    Ze = Ze + sum(sum(me.*ss, 2), 1);
    Zm = Zm + sum(sum(mm.*ss, 2), 1);

end

Z = Ze + Zm;
