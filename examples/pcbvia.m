
addpath(genpath([ pwd, '/..' ]));

% Dimensions
th = 3e-5; % trace thickness
tw = 1e-4; % trace width
d1 = 2e-4; % trace-to-plane separation
d2 = 3e-4; % plane-to-plane separation
rb = 1e-4/2; % via barrel radius
rp = 2e-4/2; % via pad radius
rh = 2e-4*1.1/2; % via hole radius
a  = 1e-3; % x-size of the enclosure/waveguide
b  = 1e-3; % y-size of the enclosure/waveguide

% Mesh options
nx = 64;   % cells along x
ny = 64;   % cells along y
n1 = 2;    % number of layers between the trace and the plane
n2 = 3;    % number of layers between the planes

% Start meshing -- prepare bitmaps
cx = nx/2 + 1; % center pixel/cell coordinates
cy = ny/2 + 1; % not necessarily integer
dx = a/nx;
dy = a/ny;
B0 = zeros(nx+2, ny+2);
BB = traceedges(drawcir(B0, cx, cy, rb/dx)); % barrel
BP = drawcir(B0, cx, cy, rp/dx); % pad
BT1 = drawline(BP, cx, cy, 0, cy, tw/dy); % trace1
BT2 = drawline(BP, cx, cy, nx+2, cy, tw/dy); % trace2
BH = 1.0 - drawcir(B0, cx, cy, rh/dx); % hole
BS = drawcir(B0, cx - nx/4, cy - ny/4, 1); % stitching
BS = drawcir(BS, cx - nx/4, cy + ny/4, 1); % stitching
BS = drawcir(BS, cx + nx/4, cy + ny/4, 1); % stitching
BS = drawcir(BS, cx + nx/4, cy - ny/4, 0.8); % stitching

% Start populating mesh
mesh.layers = struct( [] );

% Bottom trace
mesh.layers( end + 1 ) = mklayer(BT1, B0, 1);              % bottom surface
mesh.layers( end + 1 ) = mklayer(BT1, traceedges(BT1), 2); % top surface + vias

% Top trace
mesh.layers( end + 1 ) = mklayer(BT2, B0, 2 + n1 + n2 + n1); % bottom
mesh.layers( end + 1 ) = mklayer(BT2, traceedges(BT2), 2 + n1 + n2 + n1 + 1); % top + vias

% Barrel
for lp = 3:(3 + n1 + n2 + n1 - 1)
    mesh.layers( end + 1 ) = mklayer(BB*0, BB, lp);
end

% Bottom plane
mesh.layers( end + 1 ) = mklayer(BH, B0, 2 + n1);

% Top plane
mesh.layers( end + 1 ) = mklayer(BH, B0, 2 + n1 + n2);

% Stitching
for lp = (3 + n1):(3 + n1 + n2 - 1)
    mesh.layers( end + 1 ) = mklayer(BB*0, BS, lp);
end

% Merge layers one the same positions
mesh = mergelayers(mesh);

% Thicknesses of layers
ht = repmat(d1/n1, 1, n1); % trace-to-plane layers
hp = repmat(d2/n2, 1, n2); % plane-to-plane layers
h = [ th th ht hp ht th th ];

% Simulation results
Yf = [];

% Angular frequency
freqs = linspace(1e6,5e10,10)*2*pi;
freqs = 4e11;

% angular frequency
for freq = freqs

    freq
    wavelen = 2*pi/(freq * sqrt(eps0 * mu0))

    % enclosure/waveguide parameters
    wg      = wgparams(freq, a, b, h, nx, ny);
    wg.weps = h*0 + eps0; 
    wg.cnx  = 4;
    wg.cny  = 4;
    wg.Gls0 = 0; % no bottom ground
    wg.Ggr0 = 0; % no top ground

    b1 = findbases(mesh, nx, ny, 0, 0, 0, 1, @(l) l <= 2);
    b2 = findbases(mesh, nx, ny, 1, 0, 1, 1, @(l) l > 2 + n1 + n2);

    ports = { b1' b2'};
    portw = { b1'*0-dy b2'*0+dy };

    [ Y I Z ] = solvey(wg, mesh, ports, portw);

    Yf = cat(3, Yf, Y);

end

tswrite('pcbvia.y2p', freqs/(2*pi), Yf)

whos Z

% ABCD parameters
A = y2a(Y);

% Characteristic impedance and delay from simulation
Z0s = sqrt(A(1,2)./A(2,1))
tds = acos(A(1,1))./freq

c = 1/sqrt(mu0*eps0);
td = a / c % expected delay

% Calculate source voltages for the current plot
Yt = 1/100;
branches=[ 2 1 ; 3 1 ; 2 1 ; 3 1 ];
W = [ 0 ; 0 ; 1 ; 0 ];
K = 0*W;
Yc = zeros(size(branches, 1), size(branches, 1));
Yc(1:2,1:2) = Y;
Yc(3,3) = Yt;
Yc(4,4) = Yt;
[ F, V ] = solve(branches,Yc,W,K);

% excitation voltage vector
XV = zeros(size(Z,1), 1);
XV(ports{1}) = V(1);
XV(ports{2}) = V(2);

I = (Z\XV);

[ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I(:,1));
%trimesh(Tri, X, Y, Z);
trisurf(Tri, X, Y, Z, C);
xlim([ -dx  a+dx ])
ylim([ -dy  b+dy ])
zlim([ 0    a    ])

%% % The improved colormap
%% MR=[0,0; 
%%     0.02,0.3; %this is the important extra point
%%     0.3,1;
%%     1,1];

%% MG=[0,0;
%%     0.3,0; 
%%     0.7,1;
%%     1,1];

%% MB=[0,0; 
%%     0.7,0;
%%     1,1];

%% hot2 = colormapRGBmatrices(500,MR,MG,MB);
%% colormap(hot2)
%% colorbar
