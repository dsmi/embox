function iv = calc_iv_vterm(tl,iobs)
% iv = calc_iv_vterm(tl,iobs)
%
% Part of the tlines calculator (see calc_tlines), given unit voltage
% at the right terminal of a tline calculate the integral of voltage
% over the length of this tline.
% point of this tline.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   jsrc - source tline index.

denom = (1 + tl.Gls(:,iobs).*tl.t(:,iobs)).*tl.k(:,iobs);
ex1 = -exp(-tl.k(:,iobs).*(2*tl.d(:,iobs))) + tl.t1(:,iobs);
ex2 = 1 - tl.t1(:,iobs);
iv = (ex2 + tl.Gls(:,iobs).*ex1)./denom;
