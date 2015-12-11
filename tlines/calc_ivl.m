function i = calc_ivl(tl,zobs,iobs,jsrc,lsrc)
% i = calc_ivl(tl,zobs,iobs,jsrc,lsrc)
%
% Part of the tlines calculator (see calc_tlines), given the linearly
% varying distributed voltage source with intensity given as S = (z-zj)*l
% (voltage unit per length unit) which spans a particular tline, calculate
% the current at some point.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   zobs - observation coordinate.
%   iobs - observation tline index.
%   jsrc - source tline index.
%   lsrc - linear source coefficient

if iobs<jsrc,
    % Find voltage at the left terminal of the source tline
    vt = calc_vvl(tl,tl.z(:,jsrc),jsrc,jsrc,lsrc);
    % Next, find voltage at the right terminal of the observation tline
    vr = vt.*prod(tl.Tls(:,iobs+1:jsrc-1),2);
    % Finally, find current at the observation
    i = vr.*calc_i_vterm(tl,zobs,iobs);
elseif iobs>jsrc,
    % Find voltage at the right terminal of the source tline
    vt = calc_vvl(tl,tl.z(:,jsrc+1),jsrc,jsrc,lsrc);
    % Next, find voltage at the left terminal of the observation tline
    vl = vt.*prod(tl.Tgr(:,jsrc+1:iobs-1),2);
    % Finally, find current at the observation
    i = vl.*calc_i_vleft(tl,zobs,iobs);
else

    % Voltages at the left and right terminal
    V1 = calc_vvl(tl, tl.z(:,iobs), iobs, jsrc, lsrc);
    V2 = calc_vvl(tl, tl.z(:,iobs+1), iobs, jsrc, lsrc);

    t1 = tl.t1(:,iobs);
    k = tl.k(:,iobs);
    k2 = k.^2;
    Y0 = tl.Y0(:,iobs);
    Z = tl.k(:,iobs).*tl.Z0(:,iobs); % per-length impedance

    % Forward and backward voltage waves - the z coordinate is shifted such that
    % z=0 at the left terminal
    Vm = (-V2 + V1.*t1 + lsrc./k2.*(t1-1))./(t1 - 1./t1);
    Vp = V1 - Vm + lsrc./k2;

    % Forward and backward current waves
    Ip = Y0.*Vp;
    Im = -Y0.*Vm;

    % Shifted z - distance from the beginning of the observation line
    z = zobs-tl.z(:,iobs);

    % Finally the current at the observation
    i = Ip.*exp(-k.*z) + Im.*exp(k.*z) + lsrc.*z./Z;

end
