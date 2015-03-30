function test_calczmn3
% test_calczmn3
%
% One more test of calczmn2 (sligtly reorganized version of calczmn)
% Now with vias
% 

%
% test geometry - rectangular vertical loop - in y-z plane
%

nx=8;  % cells along x
ny=8;  % cells along y
a=1e-2; % x-size of the waveguide
b=1e-2; % y-size of the waveguide

% Mesh cell size
dx = a/nx;
dy = b/ny;

% angular frequency
freq=1.0e10;

nlay=nx; % number of layers (and tline sections)
layh=dx; % thickness of one layer

% parametes of the enclosure to pass to mkzmat
h=repmat(layh, 1, nlay);
wg=wgparams(freq,a,b,h,nx,ny);
wg.cnx = 4; % to make the things a little faster
wg.cny = 4;
%% wg.Gls0 = 0; % no bottom ground
%% wg.Ggr0 = 0; % no top ground


cx = nx/2 + 0.5; % center
cy = ny/2 + 0.5; 
n2 = 1; % half the width

% center coordinates for vias
ci = nx/2 - 1;
cj = ny/2 - 1;

B=zeros(nx+2,ny+2);
B=drawline(B, cx, cy-n2, cx, cy+n2, 1);
layer = mklayer(B);
layer = rmfield(layer, 'conductivity');
layer.pos = nlay/2 - n2;
mesh.layers(1) = layer;

for li=1:(n2*2)
    B=zeros(nx+2,ny+2);
    layer = mklayer(B);
    layer = rmfield(layer, 'conductivity');
    layer.vi = [ ci     ci    ];
    layer.vj = [ cj-n2  cj+n2 ];
    layer.pos = nlay/2 - n2 + li;
    mesh.layers(li + 1) = layer;
end

B=zeros(nx+2,ny+2);
B=drawline(B, cx, cy-n2, cx, cy+n2, 1);
layer = mklayer(B);
layer = rmfield(layer, 'conductivity');
layer.vi = [ ci     ci    ];
layer.vj = [ cj-n2  cj+n2 ];
layer.pos = nlay/2 + n2;
mesh.layers(n2*2 + 1) = layer;

%% Z1 = calczmn(wg, 7, 6, 7, 2, 7, 6, 6, 2)
%% Z2 = calczmn2(wg, 7, 6, 7, 2, 7, 6, 6, 2)
%% zdiff=Z1-Z2
%% Z3 = calczmn(wg, 7, 6, 6, 2, 7, 6, 6, 2)
%% Z4 = calczmn2(wg, 7, 6, 6, 2, 7, 6, 6, 2)
%% zdiff34=Z3-Z4

Zt=mkzmat2(wg, mesh, @calczmn);
Z=mkzmat2(wg, mesh, @calczmn2);

assertEquals(Zt, Z, 1e-17);
