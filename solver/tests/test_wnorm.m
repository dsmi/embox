function test_wnorm
% test_wnorm
%
% Test the normalization coefficients by evaluating integral of the mode
% vectors over the guide crossection numerically.
% 

a=0.2;
b=0.3;
maxm=3;
maxn=3;

[ ne, nm ] = wnorm(a, b, maxm, maxn);

% Eigenfunctions:
%   Psi_e=ne*cos(kx*x)*cos(ky*y)
%   Psi_m=nm*sin(kx*x)*sin(ky*y)
% Mode vectors (i,j,k are the basis)
%  ee = cross(k,grad(Psi_e)) = -i*dPsi_e/dy+j*dPsie_dx
%     = i*(ne*ky*cos(kx*x)*sin(ky*y))-j*(ne*kx*sin(kx*x)*cos(ky*y))
%  he = -grad(Psi_e) = -i*dPsi_e/dx-j*dPsi_e/dy
%     = i*(ne*kx*sin(kx*x)*cos(ky*y))+j*(ne*ky*cos(kx*x)*sin(ky*y))
%  em = -grad(Psi_m) = -i*(dPsi_m/dx)-j*dPsi_m/dy
%     = -i*(nm*kx*cos(kx*x)*sin(ky*y))-j*(nm*ky*sin(kx*x)*cos(ky*y))
%  hm = -cross(k,grad(Psi_m)) = i*dPsi_e/dy-j*dPsie_dx
%     = i*(nm*ky*sin(kx*x)*cos(ky*y))-j*(nm*kx*cos(kx*x)*sin(ky*y))

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

for m=1:(maxm),
    for n=((m==1)+1):(maxn),
	% integrage magnitude of ee
	eex = ne(m,n)*ky(m,n)*cos(kx(m,n)*x).*sin(ky(m,n)*y);
	eey = -ne(m,n)*kx(m,n)*sin(kx(m,n)*x).*cos(ky(m,n)*y);
	mee = eex.^2+eey.^2;
	imee = sum(mee(:))/(na*nb)*(a*b);
	assertEquals(1.0, imee, 1e-12);
	% integrage magnitude of he
	hex = ne(m,n)*kx(m,n)*sin(kx(m,n)*x).*cos(ky(m,n)*y);
	hey = ne(m,n)*ky(m,n)*cos(kx(m,n)*x).*sin(ky(m,n)*y);
	mhe = hex.^2+hey.^2;
	imhe = sum(mhe(:))/(na*nb)*(a*b);
	assertEquals(1.0, imhe, 1e-12);
	if m>1 && n>1,
	    % integrage magnitude of em
	    emx = nm(m,n)*kx(m,n)*cos(kx(m,n)*x).*sin(ky(m,n)*y);
	    emy = nm(m,n)*ky(m,n)*sin(kx(m,n)*x).*cos(ky(m,n)*y);
	    mem = emx.^2+emy.^2;
	    imem = sum(mem(:))/(na*nb)*(a*b);
	    assertEquals(1.0, imem, 1e-12);
	    % integrage magnitude of hm
	    hmx = nm(m,n)*ky(m,n)*sin(kx(m,n)*x).*cos(ky(m,n)*y);
	    hmy = nm(m,n)*kx(m,n)*cos(kx(m,n)*x).*sin(ky(m,n)*y);
	    mhm = hmx.^2+hmy.^2;
	    imhm = sum(mhm(:))/(na*nb)*(a*b);
	    assertEquals(1.0, imhm, 1e-12);
	end
    end
end
