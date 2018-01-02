function test_mkzmat2 
% test_mkzmat2
%
% Compare mkzmat results against the mkzmat2 which evaluates the sums
% of the waveguide modes directly instead of using fft
%
% This test uses the same geometry as test_calczmn3 - vertical loop
% in y-z plane
%

%
% test geometry - rectangular vertical loop - in y-z plane
%

nx=8;  % cells along x
ny=8;  % cells along y
a=1e-2; % x-size of the waveguide
b=2e-2; % y-size of the waveguide - intentionally different from a

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
wg.weps = eps0*(1:nlay); % vary the permittivity from layer to layer
wg.cnx = 8; % to make the things a little faster
wg.cny = 8;
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
layer.pos = nlay/2 - n2;
mesh.layers(1) = layer;

for li=1:(n2*2)
    B=zeros(nx+2,ny+2);
    layer = mklayer(B);
    layer.vi = [ ci     ci    ];
    layer.vj = [ cj-n2  cj+n2 ];
    layer.pos = nlay/2 - n2 + li;
    mesh.layers(li + 1) = layer;
end

B=zeros(nx+2,ny+2);
B=drawline(B, cx, cy-n2, cx, cy+n2, 1);
layer = mklayer(B);
layer.vi = [ ci     ci    ];
layer.vj = [ cj-n2  cj+n2 ];
layer.pos = nlay/2 + n2;
mesh.layers(n2*2 + 1) = layer;

Zt=mkzmat2(wg, mesh);
Z=mkzmat(wg, mesh);

assertEquals(Zt, Z, 1e-17);
