function v = calc_v_vleft(tl,zobs,iobs)
% v = calc_v_vleft(tl,zobs,iobs)
%
% Part of the tlines calculator (see calc_tlines), given unit voltage
% at the left terminal of a tline calculate voltage at an arbitrary
% point of this tline.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   zsrc - source coordinate.
%   jsrc - source tline index.

% right-looking reflection coefficient at the leftmost point
G = tl.Ggr(:,iobs).*tl.t(:,iobs); 
k = tl.k(:,iobs);
z = zobs - tl.z(:,iobs);
% forward- and backward-traveling waves
v = (exp(-k.*z)+exp(k.*(z - 2*tl.d(:,iobs))).*tl.Ggr(:,iobs))./(1+G);
