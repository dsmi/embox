function v = calc_vvd(tl,zobs,iobs,jsrc)
% v = calc_vvd(tl,zobs,iobs,jsrc)
%
% Part of the tlines calculator (see calc_tlines), given the unit
% distributed voltage source (voltage unit per length unit) which
% spans a particular tline, calculate the voltage at some point.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   zobs - observation coordinate.
%   iobs - observation tline index.
%   jsrc - source tline index.

if iobs<jsrc,
	% Find voltage at the left terminal of the source tline
	vt = calc_vterm_vd(tl,jsrc);
	% Next, find voltage at the right terminal of the observation tline
	vt = vt.*prod(tl.Tls(:,iobs+1:jsrc-1),2);
	% Finally, find voltage at the observation point
	v = vt.*calc_v_vterm(tl,zobs,iobs);
else if iobs>jsrc,
	% Find voltage at the right terminal of the source tline
	vt = calc_vvd(tl,tl.z(:,jsrc+1),jsrc,jsrc);
	% Next, find voltage at the left terminal of the observation tline
	vt = vt.*prod(tl.Tgr(:,jsrc+1:iobs-1),2);
	% Finally, find voltage at the observation point
	v = vt.*calc_v_vleft(tl,zobs,iobs);
else
	n = iobs;
	k = tl.k(:,n);
	z1 = tl.z(:,n);
	z2 = tl.z(:,n+1);
	d = tl.d(:,n);
	ex1 = (exp(-k.*(z2-zobs))-exp(-k.*(2.*z2-zobs-z1)));
	v1 = tl.Ggr(:,n).*ex1;
	ex2 = -(exp(k.*(2.*z1-zobs-z2))-exp(k.*(2.*z1-zobs-z1)));
	v2 = -tl.Gls(:,n).*ex2;
	ex3 = (exp(-k.*(2.*d+zobs-z2))-exp(-k.*(2.*d+zobs-z1)));
	v3 = tl.Gls(:,n).*tl.Ggr(:,n).*ex3;
	ex4 = -(exp(-k.*(2.*d-zobs+z2))-exp(-k.*(2.*d-zobs+z1)));
	v4 = -tl.Gls(:,n).*tl.Ggr(:,n).*ex4;
	ex0 = (1-exp(-k.*(zobs-z1))) + (exp(k.*(zobs-z2))-1);
	v0 = ex0;
	v = v0 + (v1+v2+v3+v4)./(1-tl.Gls(:,n).*tl.Ggr(:,n).*tl.t(:,n));
	v = v./(2*k);
end
end
