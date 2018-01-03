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

% Squared horizontal wavenumber
kx2ky2 = kx.^2+ky.^2;

% cutoff wavenumber of the waveguide geometry
kc=sqrt(kx2ky2);

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
    % We replace Z-directed electric current with a loop of magnetic current.
    % Then instead of evaluating integral M*\int -h dl (where h is the mode
    % vector) encircling the electric current we calculate integral
    % M*\int -rot(H) da over the current crossection area in accordance with
    % the curl theorem.

    % Intensity of magnetic current in dx-by-dy loop so it gives the same
    % fields as the electric dipole of unit length and intensity dx*dy
    K = -(dx*dy)./(j*freq*weps(sl)*dx*dy);

    % Integrate -M*curl(h) as described above.
    Vdm = -K.*Nm.*kx2ky2.*Gdy_flat.*Gdx_flat.*sin(kx.*xs).*sin(ky.*ys);
    Vde = Vdm*0; % No TE modes from via
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
    Ye = 1./reshape(calc_vi(tle, z(tl), tl, z(sl), sl), maxm, maxn);
    Ym = 1./reshape(calc_vi(tlm, z(tl), tl, z(sl), sl), maxm, maxn);

    % Now find mode voltages from the mode currents
    Ve = Ie./Ye;
    Vm = Im./Ym;

elseif (stype == 0 || stype == 1) && (ttype == 2)
% horizontal-to-via
    % Integral of the current over the via layer
    IIm = Im.*reshape(calc_iii(tlm, tl, z(sl), sl), maxm, maxn);
elseif (stype == 2) && (ttype == 0 || ttype == 1)
% via-to-horizontal
    % Modal V at the observation due to the via distributed source
    Ve = Vde.*reshape(calc_vvd(tle, z(tl), tl, sl), maxm, maxn);
    Vm = Vdm.*reshape(calc_vvd(tlm, z(tl), tl, sl), maxm, maxn);
else
% via-to-via
    % Currents at observation from source voltages
    IIe = Vde.*reshape(calc_iivd(tle, tl, sl), maxm, maxn);
    IIm = Vdm.*reshape(calc_iivd(tlm, tl, sl), maxm, maxn);
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

    % Calculate reaction of field from the source segment (whatever it is) on the
    % via's Z-directed electric current. Notice that K here is the same (with different sign)
    % as the intensity of magnetic current loop in the via source calculation.
    K = (dx*dy)./(j*freq*weps(tl)*dx*dy);
    Zm = K*(sum(sum(IIm.*Nm.*kx2ky2.*Gdy_flat.*Gdx_flat.*sin(kx.*xt).*sin(ky.*yt), 2), 1));
    Ze = Zm*0;

    % Via self-impedance needs to be handled separately since we replaced z-directed
    % electric current with x-y loop of magnetic current.
    % If we write the equation
    %   \oint H \cdot dl = \int (j\omega\epsilon E + J) \cdot ds
    % for a contour encircling the via we can see that the left hand side of the equation
    % does not change once we replace J with M since fields outside the source region are
    % the same. Therefore difference between E-fields inside the via when J-source and
    % M-source is used is:
    %    \int E^J \cdot ds = \int E^M \cdot ds - \frac{1}{j\omega\epsilon} \int J \cdot ds
    % Therefore for the via self-reaction we subtract the last integral multiplied by the
    % via length.
    % Another way of viewing that is as follows. Equation for the curl of H can be written
    % as:
    %    \nabla \times H = J^t
    % where J^t is the total current:
    %    J^t = j\omega\epsilon E + J^i
    % which consists of displacement current and impressed current. When we replace electric
    % current with magnetic loop in this equation we replace the impressed electric current
    % with displacement one, and the latter is due to electric field which in turn is induced
    % by magnetic current loop. Therefore when evaluating the via self-reaction we need to
    % subtract this E due to magnetic current loop from the total electric field.
    % Two important aspects:
    % 1) We use not the 'wanted' value of the source current J which is unity (and wanted
    % value of E to create the corresponding displacement current is J / (j\omega\epsilon) )
    % but its value as it is expanded as a sum over the waveguide modes that we use to expand
    % fields and currents.
    % 2) We apply this correction not only to the self-reaction of the via but to all
    % reactions  between the vias on the same layer. Although the 'wanted' value of this
    % displacement current is zero everywhere outside of the source via the expanded value
    % is slightly nonzero outside of it and therefore can affect the neighboring vias as well.
    if tl == sl && ttype == stype
        Zm = Zm - h(sl)*sum(sum(Vdm.*Nm.*Gdx_flat.*Gdy_flat.*sin(kx.*xt).*sin(ky.*yt), 2), 1);
    end

end

Z=Ze+Zm;
