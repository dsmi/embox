function iv = calc_iv_vleft(tl,iobs)
% iv = calc_iv_vleft(tl,iobs)
%
% Part of the tlines calculator (see calc_tlines), given unit voltage
% at the left terminal of a tline calculate the integral of voltage
% over the length of this tline.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   jsrc - source tline index.

% right-looking reflection coefficient at the leftmost point
G = tl.Ggr(:,iobs).*tl.t(:,iobs); 
k = tl.k(:,iobs);
% forward- and backward-traveling waves
ex1 = -tl.t1(:,iobs) + 1;
ex2 = 1./tl.t1(:,iobs) - 1;
iv = (ex1+ex2.*G)./(k.*(1+G));
