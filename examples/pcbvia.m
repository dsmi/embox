
addpath(genpath([ pwd, '/..' ]));

% Dimensions
th = 3e-5; % trace thickness
tw = 1.22e-4; % trace width
d1 = 1e-4; % trace-to-plane separation
d2 = 3e-5; % plane-to-plane separation
rb = 0.267e-3/2; % via barrel radius
rp = 0.457e-3/2; % via pad radius
rh = 0.639e-3/2; % via hole radius
a  = tw*16; % x-size of the enclosure/waveguide
b  = tw*16; % y-size of the enclosure/waveguide

% dielectric
er = 4.05;  

% Mesh options
nx = 64;   % cells along x
ny = 64;   % cells along y
n0 = 3;    % number of layers in the trace
n1 = 3;    % number of layers between the trace and the plane
n2 = 1;    % number of layers between the planes

% Start meshing -- prepare bitmaps
cx = nx/2 + 1; % center pixel/cell coordinates
cy = ny/2 + 1; % not necessarily integer
dx = a/nx;
dy = a/ny;
b0 = zeros(nx+2, ny+2);
bb = traceedges( drawcir(b0, cx, cy, rb/dx) ); % barrel
bp = drawcir(b0, cx, cy, rp/dx); % pad
bt1 = drawline(bp, cx, cy, 0, cy, tw/dy); % trace1
bt2 = drawline(bp, cx, cy, nx+2, cy, tw/dy); % trace2
bh = 1.0 - drawcir(b0, cx, cy, rh/dx); % hole
bs = traceedges(bh); % stitching

function [ mesh h epsl ] = mklinemesh( freq, nx, ny, a, b, th, tw, d1, nl, er )

    % Start meshing -- prepare bitmaps
    cx = nx/2 + 1; % center pixel/cell coordinates
    cy = ny/2 + 1; % not necessarily integer
    dx = a/nx;
    dy = b/ny;
    b0 = zeros(nx+2, ny+2);
    bt = drawline(b0, 0, cy, nx+2, cy, tw/dy); % trace
    be = traceedges(bt);

    ht = repmat( th/max( nl, 1 ), 1, nl ); % trace layers
    h = [ d1 d1 ht tw ];

    epst = repmat( eps0, 1, nl ); % trace layers
    epsd = eps0*debye( er, 0.02, 1e9, freq/(2*pi) ); % dielectric
    epsc = eps0 - j*ccopper/freq; % permittivity for copper, with conductivity

    epsl = [ epsc epsd epst eps0 ];
    
    % Create the tline mesh
    mesh = mkhull( bt, be, 2, nl );

    %% (experiment) Meshed ground
    %% mesh.layers( end + 1 ) = mklayer( b0 + 1, b0, 1, ccopper );

end


function [ Y I ] = simline( freq, nx, ny, a, b, th, tw, d1, nl, er )

    % Create mesh for simulation
    [ mesh h epsl ] = mklinemesh( freq, nx, ny, a, b, th, tw, d1, nl, er );

    % parametes of the enclosure to pass to mkzmat
    wg = wgparams( freq, a, b, h, nx, ny );
    wg.weps = epsl; 
    wg.cnx  = 4;
    wg.cny  = 4;
    wg.Gls0 = 0; % no bottom ground
    wg.Ggr0 = 0; % no top ground

    % Identify ports
    b1 = findbases( mesh, nx, ny, 0, 0, 0, 1, @( l ) l > 1 );
    b2 = findbases( mesh, nx, ny, 1, 0, 1, 1, @( l ) l > 1 );

    % Mesh cell size
    dy = wg.b/ny;

    ports = { b1' b2'};
    portw = { b1'*0-dy b2'*0+dy };
    [ Y I ]=solvey(wg, mesh, ports, portw);

end

% Start populating the via mesh
mesh.layers = struct( [] );

% Bottom trace
btm = mkhull( bt1, traceedges( bt1 ), 1, n0 );
mesh.layers = [ mesh.layers btm.layers ];

% Middle layer pad
mpm = mkhull( bp, traceedges( bp ), 1+n0+n1, n2 );
mesh.layers = [ mesh.layers mpm.layers ];

% Top trace
ttm = mkhull( bt2, traceedges( bt2 ), 1 + n0 + n1*2 + n2, n0 );
mesh.layers = [ mesh.layers, ttm.layers ];

% Barrel
brm = mkhull( bb, traceedges( bb ), 1+n0, n1*2+n2 );
mesh.layers = [ mesh.layers , brm.layers ];

% Planes
plm = mkhull( bh, bs, 1+n0+n1, n2 );
mesh.layers = [ mesh.layers , plm.layers ];

% Merge layers one the same positions
mesh = mergelayers(mesh);

% Thicknesses of the layers
ht = repmat( th/n0, 1, n0 ); % trace layers
hs = repmat( d1/n1, 1, n1 ); % trace-to-plane layers
hp = repmat( d2/n2, 1, n2 ); % plane-to-plane layers
h = [ th ht hs hp hs ht ];

% angular frequencies
freqs = linspace( 1e6, 2e10, 21 )*2*pi;

% Simulation results
Yf   = [];
Yndf = [];
Yvf  = [];
Ylf  = [];

% angular frequency
for freq = freqs

    fprintf( 'Simulation %i of %i...\n', find(freq == freqs), length(freqs) )
    wavelen = 2*pi/(freq * sqrt(eps0 * er * mu0))

    % Determine length of the line used in deembedding -- we will use it
    % later remove the feeding lines from the via results
    nx1 = round( (a/2 - rp)/dx );
    nx2 = nx1*2;

    % Obtain the port discontinuity
    fsim1 = @( ) simline( freq, nx1, ny, nx1*dx, b, th, tw, d1, n0, er );
    fsim2 = @( ) simline( freq, nx2, ny, nx2*dx, b, th, tw, d1, n0, er );
    [ D Y1 Y2 ] = deembsims( fsim1, fsim2 );
    invD = inv(D); % To be used in ed-embedding

    % De-embedded parameters of the via feed line
    Al = invD*y2abcd(Y1)*invD;
    Yl = abcd2y(Al);

    % Dielectric permeabilities of the layers
    epsd = eps0*debye( er, 0.02, 1e9, freq/(2*pi) );
    et = repmat( eps0, 1, n0 );        % trace layers
    ed = repmat( epsd, 1, n1*2 + n2 ); % dielectric layers
    weps = [ eps0 et ed et ];

    % enclosure/waveguide parameters
    wg      = wgparams(freq, a, b, h, nx, ny);
    wg.weps = weps; 
    wg.cnx  = 4;
    wg.cny  = 4;
    wg.Gls0 = 0; % no bottom ground
    wg.Ggr0 = 0; % no top ground

    b1 = findbases(mesh, nx, ny, 0, 0, 0, 1, @(l) l <= n0 + 1 );
    b2 = findbases(mesh, nx, ny, 1, 0, 1, 1, @(l) l > n0 + n1*2 + n2 );

    ports = { b1' b2' };
    portw = { b1'*0-dy b2'*0+dy };

    fprintf('Running the via simulation...\n')
    Y1 = solvey( wg, mesh, ports, portw );

    % De-embedded Y1
    invD = inv(D);
    A = invD*y2abcd(Y1)*invD;
    Y = abcd2y(A);

    % Parameters of the via alone, without the feeding lines
    invAl = inv(Al);
    Av = invAl*A*invAl;
    Yv = abcd2y(Av);

    % Calculate via L and C, T-network equivalent
    Z = inv(Yv);
    Z1 = Z(1,1) - Z(2,1); % half-series
    Z2 = Z(2,1);          % shunt

    Ct = j*imag(1/Z2)/(j*freq)
    Lt2 = j*imag(Z1)/(j*freq)
    

    % Add y-parameters for this frequency to the overall results
    Yf = cat(3, Yf, Y);
    Yndf = cat(3, Yndf, Y1);
    Ylf = cat(3, Ylf, Yl);
    Yvf = cat(3, Yvf, Yv);

end

tswrite( 'via_64.y2p',      freqs/(2*pi), Yf   )
tswrite( 'via_64_nd.y2p',   freqs/(2*pi), Yndf )
tswrite( 'via_64_line.y2p', freqs/(2*pi), Ylf  )
tswrite( 'via_64_via.y2p',  freqs/(2*pi), Yvf  )

%% % For the drawing
%% wg = wgparams(1e9, a, b, h, nx, ny);

%% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, zeros(50000,1));
%% %trimesh(Tri, X, Y, Z);
%% trisurf(Tri, X, Y, Z, C);
%% xlim([ -dx  a+dx ])
%% ylim([ -dy  b+dy ])
%% zlim([ 0    a/3    ])

%% tswrite('pcbvia.y2p', freqs/(2*pi), Yf)

%% whos Z

%% % ABCD parameters
%% A = y2a(Y);

%% % Characteristic impedance and delay from simulation
%% Z0s = sqrt(A(1,2)./A(2,1))
%% tds = acos(A(1,1))./freq

%% c = 1/sqrt(mu0*eps0);
%% td = a / c % expected delay

%% % Calculate source voltages for the current plot
%% Yt = 1/100;
%% branches=[ 2 1 ; 3 1 ; 2 1 ; 3 1 ];
%% W = [ 0 ; 0 ; 1 ; 0 ];
%% K = 0*W;
%% Yc = zeros(size(branches, 1), size(branches, 1));
%% Yc(1:2,1:2) = Y;
%% Yc(3,3) = Yt;
%% Yc(4,4) = Yt;
%% [ F, V ] = solve(branches,Yc,W,K);

%% % excitation voltage vector
%% XV = zeros(size(Z,1), 1);
%% XV(ports{1}) = V(1);
%% XV(ports{2}) = V(2);

%% I = (Z\XV);

%% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I(:,1));
%% %trimesh(Tri, X, Y, Z);
%% trisurf(Tri, X, Y, Z, C);
%% xlim([ -dx  a+dx ])
%% ylim([ -dy  b+dy ])
%% zlim([ 0    a    ])

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
