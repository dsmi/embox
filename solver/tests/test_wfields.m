function test_wfields
% test_wfields
%
% The idea of this test is as follows. We call calczmn function giving it
% some basis function, and (along with the element of the z-matrix) obtain
% the mode voltages (or amplitudes or weitghts of the modes), which we then
% use to reconstruct the fields in the waveguide and compare the discontinuity
% of the fileds across the basis current with the expected one - which we know
% because we know the source current distribution.
% 


% angular frequency
freq=1.0e11;

a=0.01; % x-size of the waveguide
b=0.01; % y-size of the waveguide
h=5e-4; % height of the metallization above ground
c=1e-3; % height of the upper ground
nx=8;  % cells along x
ny=8;  % cells along y

% parametes of the enclosure to pass to calczmn
wg=wgparams(freq,a,b,h,c,nx,ny);

% to improve the accuracy
wg.cnx=32;
wg.cny=32;

% dimensions of the basis
dx=a/nx;
dy=b/ny;

% location of the source basis
si=3;
sj=4;
x0=dx*si;
y0=dy*sj+dy/2;

% try x-directed basis first
[ Z Ve Ye Vm Ym ]=calczmn(wg, si, sj, 1, si, sj, 1, dx, dy);

% Upper limits of the waveguide mode orders to be used when computing fields
% T(E/M)[0..M-1][0..N-1] modes are used. 
maxm=size(Ve, 1);
maxn=size(Ve, 2);

% normalization coefficients for te and tm waveguide modes
[ Ne, Nm ] = wnorm(a, b, maxm, maxn);

% wm(m,n)=m-1, wn(m,n)=n-1
[ wm, wn ] = ndgrid(0:maxm-1, 0:maxn-1);

% x and y wavenumbers of the waveguide mode
kx = wm*pi./a;
ky = wn*pi./b;

for x=linspace(x0-dx*1.5, x0+dx*1.5, 7),
    for y=linspace(y0-dy*0.6, y0+dy*0.6, 4),
	hex= Ne.*kx.*sin(kx.*x).*cos(ky.*y);
	hmx= Nm.*ky.*sin(kx.*x).*cos(ky.*y);
	dhx=sum(sum(hex.*Ve.*Ye+hmx.*Vm.*Ym, 2), 1);

	hey= Ne.*ky.*cos(kx.*x).*sin(ky.*y);
	hmy=-Nm.*kx.*cos(kx.*x).*sin(ky.*y);
	dhy=sum(sum(hey.*Ve.*Ye+hmy.*Vm.*Ym, 2), 1);

	jx=-dhy; % x-directed current from the field discontinuity
	jy=dhx; % y-directed current from the field discontinuity

	if x0-dx <= x && x < x0
	    fx = (x-x0)/dx+1;
	elseif x0 <= x && x < x0+dx
	    fx = (x0-x)/dx+1;
	else
	    fx = 0;
	end

	if y0-dy/2 <= y && y < y0+dy/2
	    fy = 1;
	else
	    fy = 0;
	end

	test_jx = fx*fy;

	assertEquals(test_jx, jx, 3e-2);
	assertEquals(0, jy, 1e-15);
    end
end

% location of the source basis
si=3;
sj=4;
x0=dx*si+dx/2;
y0=dy*sj;

% try y-directed basis
[ Z Ve Ye Vm Ym ]=calczmn(wg, si, sj, 1, si, sj, 0, dx, dy);

for x=linspace(x0-dx*0.6, x0+dx*0.6, 4),
    for y=linspace(y0-dy*1.5, y0+dy*1.5, 7),
	hex= Ne.*kx.*sin(kx.*x).*cos(ky.*y);
	hmx= Nm.*ky.*sin(kx.*x).*cos(ky.*y);
	dhx=sum(sum(hex.*Ve.*Ye+hmx.*Vm.*Ym, 2), 1);

	hey= Ne.*ky.*cos(kx.*x).*sin(ky.*y);
	hmy=-Nm.*kx.*cos(kx.*x).*sin(ky.*y);
	dhy=sum(sum(hey.*Ve.*Ye+hmy.*Vm.*Ym, 2), 1);

	jx=-dhy; % x-directed current from the field discontinuity
	jy=dhx; % y-directed current from the field discontinuity

	if x0-dx/2 <= x && x < x0+dx/2
	    fx = 1;
	else
	    fx = 0;
	end

	if y0-dy <= y && y < y0
	    fy = (y-y0)/dy+1;
	elseif y0 <= y && y < y0+dy
	    fy = (y0-y)/dy+1;
	else
	    fy = 0;
	end

	test_jy = fx*fy;

	assertEquals(0, jx, 1e-15);
	assertEquals(test_jy, jy, 3e-2);
   end
end
