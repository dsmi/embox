function test_mkzmat
% test_mkzmat
%
% Compare mkzmat results against the mkzmat2 which evaluates the sums
% of the waveguide modes directly instead of using fft
%

% angular frequency
freq=1.0e11;

a=0.01; % x-size of the waveguide
b=0.01; % y-size of the waveguide
h=5e-4; % height of the metallization above ground
c=1e-3; % height of the upper ground
nx=4;  % cells along x
ny=4;  % cells along y

% parametes of the enclosure
wg=wgparams(freq,a,b,[h,c-h],nx,ny);
wg.cnx=4;
wg.cny=4;

xi=[ 1 2 3 ];
xj=[ 2 3 2 ];
yi=[ 3 1 3 2 ];
yj=[ 4 3 2 1 ];

% no vias
vi = ones(0,1);
vj = ones(0,1);

layer=struct('xi', xi, 'xj', xj, 'yi', yi, 'yj', yj, 'vi', vi, 'vj', vj);
layer.pos = 1; % position in the stackup

mesh.layers(1) = layer;
mesh.vias(1) = struct('vi', [], 'vj', [], 'pos', 1);

Zt=mkzmat2(wg, mesh);
Z=mkzmat(wg, mesh);

assertEquals(Zt, Z, 1e-15);
