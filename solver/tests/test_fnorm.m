function test_fnorm
% test_fnorm
%
% Test the squared eigenfunction integral values returned by fnorm by
% evaluating integral numerically.
% 

a=0.2;
b=0.3;
maxm=3;
maxn=3;

[ ne, nm ] = wnorm(a, b, maxm, maxn);

% Eigenfunctions:
%   Psi_e=ne*cos(kx*x)*cos(ky*y)
%   Psi_m=nm*sin(kx*x)*sin(ky*y)

na=100;
nb=100;

dx=a/na;
dy=b/nb;

[ x, y ] = ndgrid(linspace(0,a-dx,na)+dx/2, linspace(0,b-dy,nb)+dy/2);

% wm(m,n)=m-1, wn(m,n)=n-1
[ wm, wn ] = ndgrid(0:maxm-1, 0:maxn-1);

% x and y wavenumbers of the waveguide mode
kx = wm*pi./a;
ky = wn*pi./b;

test_fe = zeros(maxm, maxn);
test_fm = zeros(maxm, maxn);

for m=1:(maxm),
    for n=1:(maxn),
	% integrage magnitude of psi_e^2
	psi_e2 = (cos(kx(m,n)*x).*cos(ky(m,n)*y)).^2;
	ipsi_e2 = sum(psi_e2(:))/(na*nb)*(a*b);
	test_fe(m,n) = ne(m,n)^2*ipsi_e2;
	% integrage magnitude of psi_m^2
	psi_m2 = (sin(kx(m,n)*x).*sin(ky(m,n)*y)).^2;
	ipsi_m2 = sum(psi_m2(:))/(na*nb)*(a*b);
	test_fm(m,n) = nm(m,n)^2*ipsi_m2;
    end
end

[ fe, fm ] = fnorm(a, b, maxm, maxn);

assertEquals(test_fe, fe, 1e-15);
assertEquals(test_fm, fm, 1e-15);
