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
	vt = calc_vvd(tl,tl.z(jsrc+1),jsrc,jsrc);
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
	iex1 = (exp(k.*z2)-exp(k.*z1))./k;
	ex1 = exp(-k.*2.*z2).*exp(k.*zobs).*iex1;
	v1 = tl.Ggr(:,n).*ex1;
	iex2 = (exp(-k.*z2)-exp(-k.*z1))./(-k);
	ex2 = exp(k.*2.*z1).*exp(-k.*zobs).*iex2;
	v2 = -tl.Gls(:,n).*ex2;
	ex3 = exp(-k.*2*d).*exp(-k.*zobs).*iex1;
	v3 = tl.Gls(:,n).*tl.Ggr(:,n).*ex3;
	ex4 = exp(-k.*2.*d).*exp(k.*zobs).*iex2;
	v4 = -tl.Gls(:,n).*tl.Ggr(:,n).*ex4;
	%% if zobs>zsrc,
	%%     ex0 = exp(-k.*zobs).*exp(k.*zsrc);
	%% else
	%%     ex0 = -exp(k.*zobs).*exp(-k.*zsrc);
	%% end
	% [ z1 zobs ] + [ zobs z2 ]
	ex0 = exp(-k.*zobs).*(exp(k.*zobs)-exp(k.*z1))./k - exp(k.*zobs).*(exp(-k.*z2)-exp(-k.*zobs))./(-k);
	v0 = ex0;
	v = v0 + (v1+v2+v3+v4)./(1-tl.Gls(:,n).*tl.Ggr(:,n).*tl.t(:,n));
	v = v/2;
end
end
