function test_mkzsmat
% test_mkzsmat
%
% The mkzsmat is tested in the most straightforward way - by evaluating
% the scalar product integrals directly!
%

% angular frequency (not really needed here)
freq=1.0e11;

a=0.03; % x-size of the waveguide
b=0.01; % y-size of the waveguide
h=5e-4; % height of the metallization above ground
c=1e-3; % height of the upper ground
nx=3;   % cells along x
ny=4;   % cells along y

% parametes of the enclosure to pass to calczmn
wg=wgparams(freq,a,b,[h,c-h],nx,ny);

[ xi, xj ] = ndgrid(0:nx, 0:(ny-1));
xi = xi(:);
xj = xj(:);

[ yi, yj ] = ndgrid(0:nx, 0:ny);
yi = yi(:);
yj = yj(:);

layer = struct('xi', xi, 'xj', xj, 'yi', yi, 'yj', yj);
layer.conductivity = ccopper; % like it is make of copper

mesh=mkmesh(layer, 1);

Zs0 = (1+j)*sqrt(mu0*wg.freq/(2*layer.conductivity));

Zs=mkzsmat(wg, mesh);

% mesh cell sizes
dx=a/nx;
dy=b/ny;

nxx = length(xi);
Zxx = zeros(nxx);
for m=1:nxx,
    for n=1:nxx,
	if xj(m) == xj(n),
	    fp = @(x) ftri(x, dx*xi(m), dx)*ftri(x, dx*xi(n), dx);
	    Zxx(m,n) = Zs0*dy*quad(fp, 0, a);
	end
    end
end

nyy = length(yi);
Zyy = zeros(nyy);
for m=1:nyy,
    for n=1:nyy,
	if yi(m) == yi(n),
	    fp = @(x) ftri(x, dy*yj(m), dy)*ftri(x, dy*yj(n), dy);
	    Zyy(m,n) = Zs0*dx*quad(fp, 0, b);
	end
    end
end

% compose the entire matrix
Zxy = zeros(size(Zxx,1), size(Zyy,2));
Zyx = zeros(size(Zyy,1), size(Zxx,2));
Zst = [ Zxx Zxy ; Zyx Zyy ];

assertEquals(Zst, Zs, 1e-15);
