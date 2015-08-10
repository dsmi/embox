% NEEDS MORE WORK
% vertical transmission line

addpath(genpath([ pwd, '/..' ]));

nx=32;  % cells along x
ny=32;  % cells along y
a=1e-2; % x-size of the waveguide
b=1e-2; % y-size of the waveguide

% Mesh cell size
dx = a/nx;
dy = b/ny;

% Number of the layers (and line segments)
nlay = 20;

% Thickness of the layers
hlay = dx

% angular frequency
freq = 1e8
wavelen = 2*pi/(freq * sqrt(eps0 * mu0))

% parametes of the enclosure to pass to mkzmat
h=repmat(hlay, 1, nlay);
wg=wgparams(freq,a,b,h,nx,ny);
wg.cnx = 4; % to make the things a little faster
wg.cny = 4;
wg.Ggr0 = 0; % no top ground

% Make the mesh
clear mesh
layer = mklayer(zeros(nx+2,ny+2));
layer.vi = [ nx/4 ]; % via - the tline segment
layer.vj = [ ny/2 ];
layer = rmfield(layer, 'conductivity'); % perfect conductor
for lidx=1:nlay
    layer.pos = lidx;
    mesh.layers(lidx) = layer;
end

ports = { [ 1  ] };
portw = { [ dy ] };

[ Y I ]=solvey(wg, mesh, ports, portw);
Z1 = 1/Y

k=freq*sqrt(eps0*mu0);
nu=sqrt(mu0/eps0);
gamma=j*k;
beta=k;
d=dx;        % wire radius
h=(nx/4)*dx; % separation
l=nlay*hlay; % length
Z0=nu/(2*pi)*log(4*h/d)
Zo=-j*Z0*cot(beta*l) % open line impedance

%% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I);
%% trimesh(Tri, X, Y, Z);
%% xlim([ 0 a+2*dx ])
%% ylim([ 0 b+2*dy ])
%% zlim([ 0 hlay*nlay ])
