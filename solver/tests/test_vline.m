function test_vline
% vertical transmission line -- made of vias

nx=3;    % cells along x
ny=3;    % cells along y
a=1e-3;  % x-size of the waveguide
b=1e-3;  % y-size of the waveguide

% Mesh cell size
dx = a/nx;
dy = b/ny;

% Number of the layers and therefore number of the vias
nlay = 10;

% Thickness of the layers
hlay = dx;

% angular frequency
freq = 2e11;
wavelen = 2*pi/(freq * sqrt(eps0 * mu0));
k = freq*sqrt(eps0.*mu0);

% parametes of the enclosure to pass to mkzmat
h  = repmat(hlay, 1, nlay);
wg = wgparams(freq,a,b,h,nx,ny);
%% wg.cnx = 4;
%% wg.cny = 4;

% Center pixel coordinates
cnx = nx/2 + 1;
cny = ny/2 + 1;

% Make the mesh
clear mesh

% One via at the center
B0 = zeros(nx+2,ny+2);
BV = drawcir(B0, cnx, cny, 0.1);

for lidx=1:nlay
    layer = mklayer(B0, BV);
    layer.pos = lidx;
    mesh.layers(lidx) = layer;
end

% Ports
b1 = 1;
b2 = nlay;

ports = { b1'       b2'      };
portw = { b1'*0+dy  b2'*0-dy };

% Y-parameters of the line
Y = solvey(wg, mesh, ports, portw);

% ABCD parameters
A = y2a(Y);

% Characteristic impedance and delay from simulation
Z0s = sqrt(A(1,2)./A(2,1));
tds = acos(A(1,1))./freq;

c = 1/sqrt(mu0*eps0);
td = sum(h) / c; % expected delay
Z0 = 90; % expected Z0

assertEquals(Z0, Z0s, Z0*0.1);
assertEquals(td, tds, td*0.2);
