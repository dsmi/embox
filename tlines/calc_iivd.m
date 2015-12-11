function i = calc_iivd(tl,iobs,jsrc)
% i = calc_iivd(tl,iobs,jsrc)
%
% Part of the tlines calculator (see calc_tlines), given the unit
% distributed voltage source (voltage unit per length unit) which
% spans a particular tline, calculate the integral of current over
% the length of this or another tline.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   iobs - observation tline index.
%   jsrc - source tline index.

if iobs<jsrc,
	% Find voltage at the left terminal of the source tline
	vt = calc_vterm_vd(tl,jsrc);
	% Next, find voltage at the right terminal of the observation tline
	vt = vt.*prod(tl.Tls(:,iobs+1:jsrc-1),2);
	% Finally, find current at the observation point
	i = vt.*calc_ii_vterm(tl,iobs,0,1);
else if iobs>jsrc,
	% Find voltage at the right terminal of the source tline
	vt = calc_vvd(tl,tl.z(:,jsrc+1),jsrc,jsrc);
	% Next, find voltage at the left terminal of the observation tline
	vt = vt.*prod(tl.Tgr(:,jsrc+1:iobs-1),2);
	% Finally, find voltage at the observation point
	i = vt.*calc_ii_vleft(tl,iobs,0,1);
else
	n = iobs;
	k = tl.k(:,n);
	d = tl.d(:,n);
	t1 = tl.t1(:,n);
	t = tl.t(:,n);
	ex1 = (1-2*t1+t)./(k.*k);
	i1 = -tl.Ggr(:,n).*ex1;
	ex2 = (t-2*t1+1)./(k.*k);
	i2 = -tl.Gls(:,n).*ex2;
	ex3 = t.*(2-t1-1./t1)./(-k.*k);
	i3 = tl.Gls(:,n).*tl.Ggr(:,n).*ex3;
	ex4 = t.*(2-t1-1./t1)./(-k.*k);
	i4 = tl.Gls(:,n).*tl.Ggr(:,n).*ex4;
	ex0 = (2*d-(2-2*t1)./k)./k;
	i = ex0 + (i1+i2+i3+i4)./(1-tl.Gls(:,n).*tl.Ggr(:,n).*tl.t(:,n));
	i = i .* (tl.Y0(:,n)/2);
end
end
