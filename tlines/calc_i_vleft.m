function i = calc_i_vleft(tl,zobs,iobs)
% i = calc_i_vleft(tl,zobs,iobs)
%
% Part of the tlines calculator (see calc_tlines), given unit voltage
% at the left terminal of a tline calculate current at an arbitrary
% point of this tline.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   zobs - observation coordinate.
%   iobs - observation tline index.

% right-looking reflection coefficient at the leftmost point
G = tl.Ggr(:,iobs).*tl.t(:,iobs); 
k = tl.k(:,iobs);
z = zobs - tl.z(:,iobs);
% forward- and backward-traveling waves
i = tl.Y0(:,iobs).*(exp(-k.*z)-exp(k.*z).*G)./(1+G);
