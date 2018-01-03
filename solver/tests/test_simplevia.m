function test_simplevia
% Line - via through ref plane - line
%

% Dimensions
a  = 1e-3;  % x-size of the enclosure/waveguide
b  = 1e-3;  % y-size of the enclosure/waveguide
tw = a/8;   % trace width
d1 = a/5;   % trace-to-plane separation
rh = a/7;   % via hole radius

% Mesh options
nx = 16;   % cells along x
ny = 16;   % cells along y
n1 = 1;    % number of layers between the trace and the plane

% Start meshing -- prepare bitmaps
cx = nx/2 + 1; % center pixel/cell coordinates
cy = ny/2 + 1; % not necessarily integer
dx = a/nx;
dy = a/ny;
B0 = zeros(nx+2, ny+2);
BB = drawcir(B0, cx, cy, tw/2/dy); % barrel
BT1 = drawline(BB, cx, cy, 0, cy, tw/dy); % trace1
BT2 = drawline(BB, cx, cy, nx+2, cy, tw/dy); % trace2
BH = 1.0 - drawcir(B0, cx, cy, rh/dx); % hole

% Populate mesh
clear mesh

layer = mklayer(BT1);
layer.pos = 1;
mesh.layers(1) = layer;

layer = mklayer(BH, BB);
layer.pos = 2;
mesh.layers(2) = layer;

layer = mklayer(BT2, BB);
layer.pos = 3;
mesh.layers(3) = layer;

% Thicknesses of layers
h = [ d1 d1 d1 ];

% Angular frequency
freq = 2e11;

% enclosure/waveguide parameters
wg      = wgparams(freq, a, b, h, nx, ny);
wg.weps = h*0 + eps0; 
wg.cnx  = 8; % want it fast
wg.cny  = 8;
wg.Gls0 = 0; % no bottom ground
wg.Ggr0 = 0; % no top ground

b1 = findbases(mesh, nx, ny, 0, 0, 0, 1, @(l) l == 1);
b2 = findbases(mesh, nx, ny, 1, 0, 1, 1, @(l) l == 3);

ports = { b1' b2'};
portw = { b1'*0-dy b2'*0+dy };

[ Y I Z ] = solvey(wg, mesh, ports, portw);

% ABCD parameters
A = y2a(Y);

% Characteristic impedance and delay from simulation
Z0s = sqrt(A(1,2)./A(2,1));
tds = acos(A(1,1))./freq;

c = 1/sqrt(mu0*eps0);
td = (a + 2*d1) / c; % expected delay

Z0 = 120; % expected impedance

assertEquals(Z0, Z0s, Z0*0.05);
assertEquals(td, tds, td*0.05);
