function test_vfields
% test_vfields
%
% Similar to the test_wfields here we reconstruct the fields from the mode
% voltages and currents computed by calczmn_hollow and check if the field
% discontinuities match the source currents.
%

nx=8;   % cells along x
ny=8;   % cells along y
a=1e-2; % x-size of the waveguide
b=1e-2; % y-size of the waveguide

% Mesh cell size
dx = a/nx;
dy = b/ny;

% Number of the layers (number of the line segments is less by one)
nlay = 8;

% Thickness of the layers
hlay = dx;

% angular frequency
freq = 1e8;
wavelen = 2*pi/(freq * sqrt(eps0 * mu0));

% parametes of the enclosure to pass to mkzmat
h  = repmat(hlay, 1, nlay);
wg = wgparams(freq,a,b,h,nx,ny);
wg.cnx = 32; % to improve the accuracy
wg.cny = 32;
wg.Gls0 = 0; % no bottom ground
wg.Ggr0 = 0; % no top ground

% Center pixel coordinates
cnx = nx/2 + 1;
cny = ny/2 + 1;

% location of the source via
si = 3;
sj = 4;
sl = 4;
x0 = dx*si + dx/2;
y0 = dy*sj + dy/2;

% via
[ Z Ie Im Vdm tle tlm ] = calczmn_hollow(wg, si, sj, sl, 2, si, sj, sl, 2);

% Upper limits of the waveguide mode orders to be used when computing fields
% T(E/M)[0..M-1][0..N-1] modes are used. 
maxm=size(Ie, 1);
maxn=size(Ie, 2);

% normalization coefficients for te and tm waveguide modes
[ Ne, Nm ] = wnorm(a, b, maxm, maxn);

% wm(m,n)=m-1, wn(m,n)=n-1
[ wm, wn ] = ndgrid(0:maxm-1, 0:maxn-1);

% x and y wavenumbers of the waveguide mode
kx = wm*pi./a;
ky = wn*pi./b;

% First, check fields due to the bottom vertical part of the via
for x=linspace(x0-dx*0.6, x0+dx*0.6, 7),
    for y=linspace(y0-dy*0.6, y0+dy*0.6, 7),

	hex= Ne.*kx.*sin(kx.*x).*cos(ky.*y);
	hmx= Nm.*ky.*sin(kx.*x).*cos(ky.*y);
	dhx=sum(sum(hex.*Ie+hmx.*Im, 2), 1);

	hey= Ne.*ky.*cos(kx.*x).*sin(ky.*y);
	hmy=-Nm.*kx.*cos(kx.*x).*sin(ky.*y);
	dhy=sum(sum(hey.*Ie+hmy.*Im, 2), 1);

	jx = -dhy; % x-directed current from the field discontinuity
	jy = dhx; % y-directed current from the field discontinuity

        test_jx = fout(x, x0, dx)*fconst(y, y0, dy);
        test_jy = fconst(x, x0, dx)*fout(y, y0, dy);

	assertEquals(test_jx, jx, 6e-2);
	assertEquals(test_jy, jy, 6e-2);

    end
end
        
% layers/tlines endpoints coordinates
z = [ 0 reshape(cumsum(wg.h), 1, []) ]; 

% To obtain H-fields at the center of the vertical segment of the via we
% calculate the modal currents at the z-center of the via layer
Im = Vdm.*reshape(calc_ivd(tlm, (z(sl+1) + z(sl))./2, sl, sl), maxm, maxn);
Ie = Im*0;

% Helper functions which calculate H-field from the modal currents
hx = @(Ie, Im, x, y) sum(sum((Ne.*kx.*Ie + Nm.*ky.*Im).*sin(kx.*x).*cos(ky.*y), 2), 1);
hy = @(Ie, Im, x, y) sum(sum((Ne.*ky.*Ie - Nm.*kx.*Im).*cos(kx.*x).*sin(ky.*y), 2), 1);

%
% Vertical segment of the via is a cylinder of rectangular cros-section.
% Here we evaluate fields at the middle of each of the walls inside and
% outside, and the discontinuity of the H-field can be used to validate
% the via current:
%     Js = n x [ H1 - H2 ]
% Where Js is the surface current, n is the normal pointing towards (1),
% and H1 and H2 are the magnetic fields.
%
dhx = hx(Ie, Im, x0 - dx/2 - dx/20, y0) - hx(Ie, Im, x0 - dx/2 + dx/20, y0);
dhy = hy(Ie, Im, x0 - dx/2 - dx/20, y0) - hy(Ie, Im, x0 - dx/2 + dx/20, y0);
assertEquals(0.0, dhx, 1e-3)
assertEquals(-1.0, dhy, 0.15)

dhx = hx(Ie, Im, x0, y0 - dy/2 - dy/20) - hx(Ie, Im, x0, y0 - dy/2 + dy/20);
dhy = hy(Ie, Im, x0, y0 - dy/2 - dy/20) - hy(Ie, Im, x0, y0 - dy/2 + dy/20);
assertEquals(1.0, dhx, 0.15)
assertEquals(0.0, dhy, 1e-3)

dhx = hx(Ie, Im, x0 + dx/2 + dx/20, y0) - hx(Ie, Im, x0 + dx/2 - dx/20, y0);
dhy = hy(Ie, Im, x0 + dx/2 + dx/20, y0) - hy(Ie, Im, x0 + dx/2 - dx/20, y0);
assertEquals(0.0, dhx, 1e-3)
assertEquals(1.0, dhy, 0.15)

dhx = hx(Ie, Im, x0, y0 + dy/2 + dy/20) - hx(Ie, Im, x0, y0 + dy/2 - dy/20);
dhy = hy(Ie, Im, x0, y0 + dy/2 + dy/20) - hy(Ie, Im, x0, y0 + dy/2 - dy/20);
assertEquals(-1.0, dhx, 0.15)
assertEquals(0.0, dhy, 1e-3)

