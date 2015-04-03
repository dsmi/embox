%
% "Overclocked capacitor" - impedance of the waveguide formed by a pair of
% round plates. Compare the results against the analytical solution.
%

addpath(genpath([ pwd, '/..' ]));

nx=32;  % cells along x
ny=32;  % cells along y
a=1e-2; % x-size of the waveguide
b=1e-2; % y-size of the waveguide

% Mesh cell size
dx = a/nx;
dy = b/ny;

% angular frequency
freq=5.0e11;

% parameters of the 'capacitor' - separation, radius
d  = 1e-4; % separation
r2 = nx/3*dx; % outer radius
r1 = dx/2;
eps1 = eps0;

% parametes of the enclosure to pass to mkzmat
%h=[ dx*10 d dx*10 ];
h=[ d dx*10 ];
wg=wgparams(freq,a,b,h,nx,ny);
%% wg.cnx = 8; % to make the things a little faster
%% wg.cny = 8;
%wg.Gls0 = 0; % no bottom ground
%wg.Ggr0 = 0; % no top ground


B=zeros(nx+2,ny+2);
B=drawcir(B, nx/2+1+0.5, nx/2+1+0.5, r2/dx);

% imagesc(B)

%% % bottom layer - no via
%% clear mesh
%% layer = mklayer(B);
%% layer.pos = 1;
%% mesh.layers(1) = layer;

%% % top layer - with via (which will be a port)
%% layer.pos = 2;
%% layer.vi=[ nx/2 ];
%% layer.vj=[ ny/2 ];
%% mesh.layers(2) = layer;

% one and only layer - with via to the ground (which will be a port)
clear mesh
layer = mklayer(B);
layer.pos = 1;
layer.vi=[ nx/2 ];
layer.vj=[ ny/2 ];
mesh.layers(1) = layer;

%% [ Tri, X, Y, Z ] = mesh2tri(wg, mesh);
%% trimesh(Tri,X,Y,Z);
%% xlim([ 0 a+2*dx ])
%% ylim([ 0 b+2*dy ])
%% zlim([ 0 dx*8 ])

%Z=mkzmat(wg, mesh);

nxy = size(layer.xi) + size(layer.yi); % number of x- and y- directed
ports = { [ nxy + 1 ] }; % the very last b.f. is the via
portw = { [ dy      ] };

[ Y I ]=solvey(wg, mesh, ports, portw);
Z1 = 1/Y

function Z=rz(freq, eps1, r1, r2, d)

    % Parameters of the plane.
    Yplane = j*freq*eps1/d;
    Zplane = j*freq*mu0*d;
    k = sqrt(-Yplane*Zplane);

    % Analytical solution
    M = [ besselj(0, k*r1) bessely(0, k*r1) ; ...
	    besselj(0, k*r2) bessely(0, k*r2) ];
    %M = [ besselh(0,2,k*r1) besselh(0,1,k*r1) ; ...
    %      besselh(0,2,k*r2) besselh(0,1,k*r2) ]

    dM = [ -besselj(1, k*r1)*k*r1 -bessely(1, k*r1)*k*r1 ; ...
	   -besselj(1, k*r2)*k*r2 -bessely(1, k*r2)*k*r2 ];
    %dM = [ -besselh(1,2,k*r1)*k -besselh(1,1,k*r1)*k ; ...
    %       -besselh(1,2,k*r2)*k -besselh(1,1,k*r2)*k ]

    b = [ 1 ; 0 ];
    x = dM\b;

    Zplane = j*freq*mu0*d;
    V = M*x;
    I = -2*pi/Zplane;
    Z = V(1)/I;
end

Z0=rz(freq, eps1, r1, r2, d)

wavelen = 2*pi/(freq * sqrt(eps0 * mu0))

