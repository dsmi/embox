function i = calc_iii(tl,iobs,a,b,zsrc,jsrc)
% i = calc_iii(tl,iobs,a,b,zsrc,jsrc)
%
% Part of the tlines calculator (see calc_tlines), given the unit
% shunt current source at an arbitrary point, calculate the integral
% of current over the length of another tline multiplied by the linear
% function a*(z-zi)+b.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   iobs - observation tline index.
%   a, b - coefficients of the linear function to multiply the current
%          by under the inegral sign
%   zsrc - source coordinate.
%   jsrc - source tline index.

if iobs<jsrc,
	% Find voltage at the left terminal of the source tline
	vt = calc_vterm_i(tl,zsrc,jsrc);
	% Next, find voltage at the right terminal of the observation tline
	vt = vt.*prod(tl.Tls(:,iobs+1:jsrc-1),2);
	% Finally, find current at the observation point
	i = vt.*calc_ii_vterm(tl,iobs,a,b);
elseif iobs>jsrc,
	% Find voltage at the right terminal of the source tline
	vt = calc_vi(tl,tl.z(:,jsrc+1),jsrc,zsrc,jsrc);
	% Next, find voltage at the left terminal of the observation tline
	vt = vt.*prod(tl.Tgr(:,jsrc+1:iobs-1),2);
	% Finally, find voltage at the observation point
	i = vt.*calc_ii_vleft(tl,iobs,a,b);
else
	n = iobs;
	k = tl.k(:,n);
	d = tl.d(:,n);
	z1 = tl.z(:,n);
	z2 = tl.z(:,n+1);
        ak = a./k;
	ex1 = exp(k.*(zsrc-z2)).*(a.*d-ak+b) + exp(k.*(zsrc-z2-d)).*(ak-b);
	i1 = -tl.Ggr(:,n).*ex1;
	ex2 = (ak+b).*exp(-k.*(zsrc-z1)) - exp(-k.*(zsrc-z1+d)).*(a.*d+ak+b);
	i2 = tl.Gls(:,n).*ex2;
	ex3 = (ak+b).*exp(k.*(zsrc+z1-2*z2)) - exp(k.*(zsrc-2*d-z2)).*(a.*d+ak+b);
	i3 = tl.Gls(:,n).*tl.Ggr(:,n).*ex3;
	ex4 = exp(-k.*(zsrc+2*d-z2)).*(a.*d-ak+b) + exp(-k.*(zsrc+2*d-z1)).*(ak-b);
	i4 = -tl.Gls(:,n).*tl.Ggr(:,n).*ex4;
        ex01 = exp(-k.*(zsrc-z1));
        ex02 = exp(-k.*(z2-zsrc));
        i0 = -(ak-b).*ex01 - (a.*d+ak+b).*ex02 + 2*ak;
	i = i0 + (i1+i2+i3+i4)./(1-tl.Gls(:,n).*tl.Ggr(:,n).*tl.t(:,n));
	i = i./(2*k);
end
