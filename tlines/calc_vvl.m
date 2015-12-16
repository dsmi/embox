function v = calc_vvl(tl,zobs,iobs,jsrc,lsrc)
% v = calc_vvl(tl,zobs,iobs,jsrc,lsrc)
%
% Part of the tlines calculator (see calc_tlines), given the linearly
% varying distributed voltage source with intensity given as S = (z-zj)*l
% (voltage unit per length unit) which spans a particular tline, calculate
% the voltage at some point.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   zobs - observation coordinate.
%   iobs - observation tline index.
%   jsrc - source tline index.
%   lsrc - linear source coefficient

if iobs<jsrc,
	% Find voltage at the left terminal of the source tline
	vt = calc_vvl(tl,tl.z(:,jsrc),jsrc,jsrc,lsrc);
	% Next, find voltage at the right terminal of the observation tline
	vt = vt.*prod(tl.Tls(:,iobs+1:jsrc-1),2);
	% Finally, find voltage at the observation point
	v = vt.*calc_v_vterm(tl,zobs,iobs);
else if iobs>jsrc,
	% Find voltage at the right terminal of the source tline
	vt = calc_vvl(tl,tl.z(:,jsrc+1),jsrc,jsrc,lsrc);
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
	ex1 = lsrc.*(exp(-k.*(z2-zobs)).*(k.*d-1)+exp(-k.*(z2-zobs+d)));
	v1 = tl.Ggr(:,n).*ex1;
	ex2 = lsrc.*exp(-k.*(zobs-z1+d)).*(-k.*d+exp(k.*d)-1);
	v2 = -tl.Gls(:,n).*ex2;
	ex3 = lsrc.*(exp(-k.*(zobs-z2+2*d)).*(k.*d-1)+exp(-k.*(zobs-z1+2*d)));
	v3 = tl.Gls(:,n).*tl.Ggr(:,n).*ex3;
	ex4 = lsrc.*(exp(k.*(zobs-z2-2*d)).*(-k.*d-1)+exp(k.*(zobs-z2-d)));
	v4 = -tl.Gls(:,n).*tl.Ggr(:,n).*ex4;
	ex0 = lsrc.*(exp(-k.*(zobs-z1)) + exp(-k.*(z2-zobs)).*(d.*k + 1) - 2);
	v0 = ex0;
	v = v0 + (v1+v2+v3+v4)./(1-tl.Gls(:,n).*tl.Ggr(:,n).*tl.t(:,n));
	v = v./(2*k.*k);
end
end
