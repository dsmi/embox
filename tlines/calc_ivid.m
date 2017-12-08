function iv = calc_ivid(tl,iobs,jsrc)
% iv = calc_ivid(tl,iobs,jsrc)
%
% Part of the tlines calculator (see calc_tlines), given the unit
% distributed current source (ampere unit per length unit) which
% spans a particular tline, calculate the integral of voltage over
% the length of this or another tline.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   iobs - observation tline index.
%   jsrc - source tline index.

if iobs<jsrc,
    % Find voltage at the left terminal of the source tline
    vt = calc_vid(tl,tl.z(:,jsrc),jsrc,jsrc);
    % Next, find voltage at the right terminal of the observation tline
    vt = vt.*prod(tl.Tls(:,iobs+1:jsrc-1),2);
    % Finally, find integral of voltage
    iv = vt.*calc_iv_vterm(tl,iobs);
else if iobs>jsrc
    % Find voltage at the right terminal of the source tline
    vt = calc_vid(tl,tl.z(:,jsrc+1),jsrc,jsrc);
    % Next, find voltage at the left terminal of the observation tline
    vt = vt.*prod(tl.Tgr(:,jsrc+1:iobs-1),2);
    % Finally, find integral of voltage
    iv = vt.*calc_iv_vleft(tl,iobs);
else
    n = iobs;
    k = tl.k(:,n);
    d = tl.d(:,n);
    t1 = tl.t1(:,n);
    t = tl.t(:,n);
    ex1 = 1 - 2*t1 + t;
    v1 = tl.Ggr(:,n).*ex1;
    v2 = tl.Gls(:,n).*ex1;
    ex3 = -t.*(2-t1-1./t1);
    v3 = tl.Gls(:,n).*tl.Ggr(:,n).*ex3;
    v4 = tl.Gls(:,n).*tl.Ggr(:,n).*ex3;
    ex0 = 2.*k.*d + 2*(t1 - 1); % Direct ray
    v = ex0 + (v1+v2+v3+v4)./(1-tl.Gls(:,n).*tl.Ggr(:,n).*tl.t(:,n));
    iv = v .* (tl.Z0(:,n)./(2*k.*k));
end
end
