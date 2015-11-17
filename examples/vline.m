% vertical transmission line
% WORK IN PROGRESS

addpath(genpath([ pwd, '/..' ]));

nx=5;   % cells along x
ny=5;   % cells along y
a=1e-2; % x-size of the waveguide
b=1e-2; % y-size of the waveguide

% Mesh cell size
dx = a/nx;
dy = b/ny;

% Number of the layers (number of the line segments is less by one)
nlay = 10;

% Thickness of the layers
hlay = dx;

% angular frequency
freq = 1e7
wavelen = 2*pi/(freq * sqrt(eps0 * mu0))

% parametes of the enclosure to pass to mkzmat
h  = repmat(hlay, 1, nlay);
wg = wgparams(freq,a,b,h,nx,ny);
wg.cnx = 4; % to make the things a little faster
wg.cny = 4;
wg.Gls0 = 0; % no bottom ground
wg.Ggr0 = 0; % no top ground

% All zeros - we will draw here
B0 = zeros(nx+2,ny+2);

% Center pixel coordinates (not quite center because the # of pixels is even)
cnx = nx/2 + 1;%0.5;
cny = ny/2 + 1;%0.5;

% Half of the conductor offset
dy2 = 2;

% Line which connects the vertical conductors
B1 = drawline(B0, cnx, 0, cnx, cny, 1);

% Vias for the vertical conductors
B2 = drawline(B0, cnx, cny, cnx, cny + 1e-10, 1);

% Make the mesh
clear mesh

% Bottommost layer - no vias
layer = mklayer(B1, 0*B2);
layer = rmfield(layer, 'conductivity'); % perfect conductor
layer.pos = 1;
mesh.layers(1) = layer;

for lidx=2:nlay-1
    % only vias
    layer = mklayer(B2*0, B2);
    layer = rmfield(layer, 'conductivity'); % perfect conductor
    layer.pos = lidx;
    mesh.layers(lidx) = layer;
end

% topmost layer - both vias and connection
layer = mklayer(B1, B2);
layer = rmfield(layer, 'conductivity'); % perfect conductor
layer.pos = nlay;
mesh.layers(nlay) = layer;

% Ports
b1 = findbases(mesh, nx, ny, 0, 0, 1, 0, @(l) l == 1)
b2 = findbases(mesh, nx, ny, 0, 0, 1, 0, @(l) l == nlay)

I = zeros(1, 10000);
[ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I);
trimesh(Tri, X, Y, Z);
xlim([ 0      wg.a+3*(wg.a/nx) ])
ylim([ 0      wg.b+3*(wg.b/ny) ])
zlim([ -hlay  hlay * (nlay + 1) ])

ports = { b1'       b2' };
portw = { b1'*0-dy  b2'*0-dy };
%% [ Y I Z ]=solvey(wg, mesh, ports, portw);
%% Y

Z = mkzmat(wg, mesh);
%Z = mkzmat2(wg, mesh, @calczmn);

% Find the admittance matrix
N=size(Z,1);
np = length(ports);
portwi = cell2mat(portw);
[ faceidx fportidx ] = ports2subs(ports);
V = zeros(N,np);
V(sub2ind(size(V), faceidx, fportidx)) = portwi;
I = Z\V;
Y = -V.'*I;

A=y2a(Y);

% Characteristic impedance and delay from simulation
N=size(A,1)/2; % number of ports on each side
A12=A(1:N,N+1:end);
A21=A(N+1:end,1:N);
Z0s=sqrt(A12*inv(A21))
tds=acos(A(1,1))./freq

c = 1/sqrt(mu0*eps0);
td = (nlay - 1) * hlay / c
freq2 = 2*pi/(2*td) % 'half-wave' frequency
