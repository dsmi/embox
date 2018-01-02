function i = calc_ii_vleft(tl,iobs)
% i = calc_ii_vleft(tl,iobs)
%
% Part of the tlines calculator (see calc_tlines), given unit voltage
% at the left terminal of a tline calculate the integral of current
% over the length of this tline.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   iobs - observation tline index.

% right-looking reflection coefficient at the leftmost point
G = tl.Ggr(:,iobs).*tl.t(:,iobs);
ex1 = 1-tl.t1(:,iobs);
ex2 = tl.Ggr(:,iobs).*tl.t1(:,iobs) - G;
% forward- and backward-traveling waves
i = tl.Y0(:,iobs).*(ex1-ex2)./((1+G).*tl.k(:,iobs));
