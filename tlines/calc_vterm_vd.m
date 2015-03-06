function v = calc_vterm_vd(tl,jsrc)
% v = calc_vterm_vd(tl,jsrc)
%
% Part of the tlines calculator (see calc_tlines), given the unit
% distributed voltage source (voltage unit per length unit) which
% spans a particular tline, calculate the voltage at the left
% terminal of this tline.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   jsrc - source tline index. The voltage source spans the whole tline.

% Based on calc_vterm_v (see there for details how the expressions have beed
% obtained), the source is integrated over the tline
Yls=tl.Y0(:,jsrc).*(tl.Gls(:,jsrc)-1)./(tl.Gls(:,jsrc)+1);
Ygr=tl.Y0(:,jsrc-1).*(1-tl.Ggr(:,jsrc-1))./(1+tl.Ggr(:,jsrc-1));
M=Ygr./((Ygr-Yls).*(1-tl.Ggr(:,jsrc).*tl.t(:,jsrc)));

k = tl.k(:,jsrc);
d = tl.d(:,jsrc);
ex2 = 1-exp(-k.*d);
ex3 = exp(-k.*d)-exp(-k.*2.*d);

v = -M.*(ex2-tl.Ggr(:,jsrc).*ex3)./k;
