
addpath(genpath([ pwd, '/..' ]));

% Dimensions
tw = 1e-3/15; % trace width
d1 = 2e-4;    % trace-to-plane separation
rh = 2e-4/2;  % via hole radius
a  = 1e-3; % x-size of the enclosure/waveguide
b  = 1e-3; % y-size of the enclosure/waveguide

% Mesh options
nx = 32;   % cells along x
ny = 32;   % cells along y
dx = a/nx;
dy = a/ny;

B = ones(nx+2, ny+2);

% Populate mesh
clear mesh

layer = mklayer(B);
layer.pos = 1;
mesh.layers(1) = layer;

% enclosure/waveguide parameters
wg      = wgparams(1, a, b, [ a/2 a/2 ], nx, ny);

%I = ones(10000,1);
I = [ (layer.xi+layer.xj*0)(:)-20 ; (layer.yi*0+layer.yj)(:)-20 ];
[ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I(:,1));
%trimesh(Tri, X, Y, Z);
trisurf(Tri, X, Y, Z, C, "facecolor", "interp");
xlim([ -dx  a+dx*2 ])
ylim([ -dy  b+dy*2 ])
zlim([ 0    a ])
xlabel('X')
ylabel('Y')
zlabel('Z')

% The improved colormap
MR=[0,0; 
    0.02,0.3; %this is the important extra point
    0.3,1;
    1,1];

MG=[0,0;
    0.3,0; 
    0.7,1;
    1,1];

MB=[0,0; 
    0.7,0;
    1,1];

hot2 = colormapRGBmatrices(500,MR,MG,MB);

colormap(hot2)
colorbar
