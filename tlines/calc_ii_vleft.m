function ii = calc_ii_vleft(tl,iobs,a,b)
% ii = calc_ii_vleft(tl,iobs,a,b)
%
% Part of the tlines calculator (see calc_tlines), given unit voltage
% at the left terminal of a tline calculate the integral of current
% multiplied by the linear function a*(z-zi)+b over the length of this
% tline.
% over the length of this tline.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   iobs - observation tline index.
%   a, b - coefficients of the linear function to multiply the current
%          by under the inegral sign

% right-looking reflection coefficient at the leftmost point
G = tl.Ggr(:,iobs).*tl.t(:,iobs);

Y0 = tl.Y0(:,iobs);
d = tl.d(:,iobs);
k = tl.k(:,iobs);
ak = a./k;

ex1 = -exp(-k.*d).*(a.*d + ak + b) + ak + b;
ex2 = exp(k.*d).*(a.*d - ak + b) + ak - b;
ii = (ex1 - G.*ex2).*Y0./(k.*(G+1));
