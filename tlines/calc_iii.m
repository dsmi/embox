function i = calc_iii(tl,iobs,zsrc,jsrc)
% i = calc_iii(tl,iobs,zsrc,jsrc)
%
% Part of the tlines calculator (see calc_tlines), given the unit
% shunt current source at an arbitrary point, calculate the integral
% of current over the length of another tline.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   iobs - observation tline index.
%   zsrc - source coordinate.
%   jsrc - source tline index.

if iobs<jsrc,
	% Find voltage at the left terminal of the source tline
	vt = calc_vterm_i(tl,zsrc,jsrc);
	% Next, find voltage at the right terminal of the observation tline
	vt = vt.*prod(tl.Tls(:,iobs+1:jsrc-1),2);
	% Finally, find current at the observation point
	i = vt.*calc_ii_vterm(tl,iobs);
else if iobs>jsrc,
	% Find voltage at the right terminal of the source tline
	vt = calc_vi(tl,tl.z(:,jsrc+1),jsrc,zsrc,jsrc);
	% Next, find voltage at the left terminal of the observation tline
	vt = vt.*prod(tl.Tgr(:,jsrc+1:iobs-1),2);
	% Finally, find voltage at the observation point
	i = vt.*calc_ii_vleft(tl,iobs);
else
	n = iobs;
	k = tl.k(:,n);
	d = tl.d(:,n);
	z1 = tl.z(:,n);
	z2 = tl.z(:,n+1);
	ex1 = (exp(-k.*(z2-zsrc))-exp(-k.*(2.*z2-zsrc-z1)));
	i1 = -tl.Ggr(:,n).*ex1;
	ex2 = -(exp(k.*(2.*z1-zsrc-z2))-exp(k.*(z1-zsrc)));
	i2 = tl.Gls(:,n).*ex2;
	ex3 = -(exp(-k.*(2.*d-zsrc+z2))-exp(-k.*(2.*d-zsrc+z1)));
	i3 = tl.Gls(:,n).*tl.Ggr(:,n).*ex3;
	ex4 = (exp(-k.*(2.*d+zsrc-z2))-exp(-k.*(2.*d+zsrc-z1)));
	i4 = -tl.Gls(:,n).*tl.Ggr(:,n).*ex4;
	i0 = (exp(k.*(z1-zsrc)) - exp(-k.*(z2-zsrc)));
	i = i0 + (i1+i2+i3+i4)./(1-tl.Gls(:,n).*tl.Ggr(:,n).*tl.t(:,n));
	i = i./(2*k);
end
end
