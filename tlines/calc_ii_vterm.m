function i = calc_ii_vterm(tl,iobs)
% i = calc_ii_vterm(tl,iobs)
%
% Part of the tlines calculator (see calc_tlines), given unit voltage
% at the right terminal of a tline calculate the integral of current
% over the length of this tline.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   iobs - observation tline index.

denom = (1 + tl.Gls(:,iobs).*tl.t(:,iobs)).*tl.k(:,iobs);
ex1 = tl.t1(:,iobs)-tl.t(:,iobs);
ex2 = 1-tl.t1(:,iobs);
i = -tl.Y0(:,iobs).*(ex2 - tl.Gls(:,iobs).*ex1)./denom;
