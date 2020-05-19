
addpath(genpath([ pwd, '/..' ]));

% Dimensions
tw = 1.22e-4; % trace width
d1 = 1e-4; % trace-to-plane separation
d2 = 3e-5; % plane-to-plane separation
rb = 0.267e-3/2; % via barrel radius
rp = 0.457e-3/2; % via pad radius
rh = 0.639e-3/2; % via hole radius
a  = tw*16; % x-size of the enclosure/waveguide
b  = tw*16; % y-size of the enclosure/waveguide

% Mesh options
nx = 64;   % cells along x
ny = 64;   % cells along y
n0 = 3;    % number of layers in the trace
n1 = 3;    % number of layers between the trace and the plane
n2 = 1;    % number of layers between the planes
n3 = 6;    % number of layers in the 'thick' dielectric
%% n0 = n1 = n2 = n3 = 1; % Temp rough mesh

% Parameters of the via layers -- from bottom to top
%  layer thickness
lh = [60.96 91.44 30.48 116.84 30.48 914.4 30.48 116.84 30.48 91.44 60.96]*1e-6;
ll  = [  n0,   n1,   n2,    n1,   n2,   n3,   n0,    n1,   n2,   n1,   n0 ];
ler = [ 1.0, 4.05,  4.5,   4.5,  4.5,  4.5,  4.5,   4.5,  4.5, 4.05,  1.0 ];
llt = [ 0.0, 0.02, 0.02,  0.02, 0.02, 0.02, 0.02,  0.02, 0.02, 0.02,  0.0 ];

function [ mesh wg ] = mkmicrostrip( freq, nx, ny, a, b, th, w, h, nl, er )

    % Start meshing -- prepare bitmaps
    cx = nx/2 + 1; % center pixel/cell coordinates
    cy = ny/2 + 1; % not necessarily integer
    dx = a/nx;
    dy = b/ny;
    b0 = zeros(nx+2, ny+2);
    bt = drawline(b0, 0, cy, nx+2, cy, w/dy); % trace
    be = traceedges(bt);

    % Create the tline mesh
    mesh = mkhull( bt, be, 2, nl );

    epsd = eps0*debye( er, 0.02, 1e9, freq/(2*pi) ); % dielectric
    epsc = eps0 - j*ccopper/freq; % permittivity for copper, with conductivity

    h    = [ h    h    repmat( th/max( nl, 1 ), 1, nl )  h    h    ];
    weps = [ epsc epsd repmat( eps0, 1, nl )             eps0 eps0 ];

    % parametes of the enclosure to pass to mkzmat
    wg = wgparams( freq, a, b, h, nx, ny );
    wg.weps = weps;
    wg.cnx  = 4;
    wg.cny  = 4;
    wg.Gls0 = 0; % no bottom ground
    wg.Ggr0 = 0; % no top ground

end

function [ mesh wg ] = mkstripline( freq, nx, ny, a, b, th, w, h1, h2, nl, er )

    % Start meshing -- prepare bitmaps
    cx = nx/2 + 1; % center pixel/cell coordinates
    cy = ny/2 + 1; % not necessarily integer
    dx = a/nx;
    dy = b/ny;
    b0 = zeros(nx+2, ny+2);
    bt = drawline(b0, 0, cy, nx+2, cy, w/dy); % trace
    be = traceedges(bt);

    % Create the tline mesh
    mesh = mkhull( bt, be, 2, nl );

    epsd = eps0*debye( er, 0.02, 1e9, freq/(2*pi) ); % dielectric
    epsc = eps0 - j*ccopper/freq; % permittivity for copper, with conductivity

    h    = [ h1   h1    repmat( th/max( nl, 1 ), 1, nl )   h2    h2   ];
    weps = [ epsc epsd  repmat( eps0, 1, nl )              epsd  epsc ];

    % parametes of the enclosure to pass to mkzmat
    wg = wgparams( freq, a, b, h, nx, ny );
    wg.weps = weps;
    wg.cnx  = 4;
    wg.cny  = 4;
    wg.Gls0 = 0; % no bottom ground
    wg.Ggr0 = 0; % no top ground

end


function [ Y I ] = simline( fmesh )

    % Create mesh for simulation
    [ mesh wg ] = fmesh( );

    % Identify ports
    b1 = findbases( mesh, wg.nx, wg.ny, 0, 0, 0, 1, @( l ) l > 1 );
    b2 = findbases( mesh, wg.nx, wg.ny, 1, 0, 1, 1, @( l ) l > 1 );

    % Mesh cell size
    dy = wg.b/wg.ny;

    ports = { b1' b2'};
    portw = { b1'*0-dy b2'*0+dy };
    [ Y I ]=solvey(wg, mesh, ports, portw);

end

% Start meshing -- prepare bitmaps
cx = nx/2 + 1; % center pixel/cell coordinates
cy = ny/2 + 1; % not necessarily integer
dx = a/nx;
dy = a/ny;
b0 = zeros(nx+2, ny+2);
bb = traceedges( drawcir(b0, cx, cy, rb/dx) ); % barrel
bp = drawcir(b0, cx, cy, rp/dx); % pad
bt1 = drawline(bp, cx, cy, 0, cy, tw/dy);    % trace1,      0 -> xmax/2
bt2 = drawline(bp, cx, cy, nx+2, cy, tw/dy); % trace2, xmax/2 -> xmax
bh = 1.0 - drawcir(b0, cx, cy, rh/dx); % hole
bs = traceedges(bh); % stitching

% Bitmaps for each of the layer, bottom to top
lb{ 1 } = bp;                          % pad
le{ 1 } = traceedges( bp );
lb{ 2 } = bb;                          % barrel
le{ 2 } = traceedges( bb );
lb{ 3 } = double( bp | bh );           % pad+plane
le{ 3 } = traceedges( double( bp | bh ) );
lb{ 4 } = bb;                          % barrel
le{ 4 } = traceedges( bb );
lb{ 5 } = bp;                          % pad
le{ 5 } = traceedges( bp );
lb{ 6 } = bb;                          % barrel
le{ 6 } = traceedges( bb );
lb{ 7 } = bt2;                         % trace 2
le{ 7 } = traceedges( bt2 );
lb{ 8 } = bb;                          % barrel
le{ 8 } = traceedges( bb );
lb{ 9 } = double( bp | bh );           % pad+plane
le{ 9 } = traceedges( double( bp | bh ) );
lb{ 10 } = bb;                          % barrel
le{ 10 } = traceedges( bb );
lb{ 11 } = bt1;                         % trace 1
le{ 11 } = traceedges( bt1 );

% Start populating the mesh
mesh.layers = struct( [] );

h   = [ 1e-5 ]; % bottom simulator layer
wer = [ 1.0 ];
wlt = [ 0.0 ];
lastl = 1; % Last processed simulator layer
for lidx = 1:length(lh)
    numl = ll( lidx ); % number of the simulator layers
    laymesh = mkhull( lb{ lidx }, le{ lidx }, lastl, numl );
    lastl = lastl + numl;
    mesh.layers = [ mesh.layers laymesh.layers ];
    h   = [ h   repmat( lh( lidx )/numl, 1, numl ) ];
    wer = [ wer repmat( ler( lidx ),     1, numl ) ];
    wlt = [ wlt repmat( llt( lidx ),     1, numl ) ];
end

% Merge layers one the same positions
mesh = mergelayers(mesh);

% angular frequencies
freqs = linspace( 1e6, 2e10, 21 )*2*pi;

% Simulation results
Yf    = [];
Yndf  = [];
Yvf   = [];
Yl1f  = [];
Yl2f  = [];

% angular frequency
for freq = freqs

    fprintf( 'Simulation %i of %i...\n', find(freq == freqs), length(freqs) )
    wavelen = 2*pi/(freq * sqrt(eps0 * ler(3) * mu0))

    % Determine length of the line used in deembedding -- we will use it
    % later remove the feeding lines from the via results
    nx1 = round( (a/2 - rp)/dx );
    nx2 = nx1*2;
    a1 = nx1*dx; % length to use in the deembedding
    a2 = nx2*dx; % 2*l for the deembedding

    % Parameters of the microstrip for the deembedding
    h1 = lh(1); % thickness of the trace
    h2 = lh(2); % trace-to-plane separation
    er = ler(2); % dielectric

    % Obtain the port discontinuity -- first line, microstrip
    fmesh1 = @( ) mkmicrostrip( freq, nx1, ny, a1, b, h1, tw, h2, n0, er );
    fmesh2 = @( ) mkmicrostrip( freq, nx2, ny, a2, b, h1, tw, h2, n0, er );
    [ D1 Y1 Y2 ] = deembsims( @( ) simline( fmesh1 ), @( ) simline( fmesh2 ) );

    % De-embedded parameters of the first via feed line
    Al1 = inv(D1)*y2abcd(Y1)*inv(D1);
    Yl1 = abcd2y(Al1);

    % Parameters of the stripline for the deembedding
    h1 = lh(5); % thickness of the trace
    h2 = lh(4); % trace-to-plane separation
    h3 = sum(lh(6:8)); % trace-to-plane separation
    er = ler(4); % dielectric

    % Obtain the port discontinuity -- second line, stripline
    fmesh1 = @( ) mkstripline( freq, nx1, ny, a1, b, h1, tw, h2, h3, n0, er );
    fmesh2 = @( ) mkstripline( freq, nx2, ny, a2, b, h1, tw, h2, h3, n0, er );
    [ D2 Y1 Y2 ] = deembsims( @( ) simline( fmesh1 ), @( ) simline( fmesh2 ) );

    % De-embedded parameters of the first via feed line
    Al2 = inv(D2)*y2abcd(Y1)*inv(D2);
    Yl2 = abcd2y(Al2);

    % Dielectric permeabilities of the layers
    weps = eps0*debye( wer, wlt, 1.0e9, freq/(2*pi) );

    % enclosure/waveguide parameters
    wg      = wgparams(freq, a, b, h, nx, ny);
    wg.weps = weps; 
    wg.cnx  = 4;
    wg.cny  = 4;
    wg.Gls0 = 0; % no bottom ground
    wg.Ggr0 = 0; % no top ground

    b1 = findbases( mesh, nx, ny, 0, 0, 0, 1, @(l) l > sum( ll )-n0 );

    l0 = sum(ll(1:6));   % starting layer of the second trace
    l1 = sum(ll(1:7))+1; % ending layer of the second trace
    b2 = findbases( mesh, wg.nx, wg.ny, 1, 0, 1, 1, @( l ) l > l0 & l <= l1 );

    ports = { b1' b2' };
    portw = { b1'*0-dy b2'*0+dy };

    fprintf('Running the via simulation...\n')
    [ Y1 I ] = solvey( wg, mesh, ports, portw );

    % De-embedded Y1
    A = inv(D1)*y2abcd(Y1)*inv(D2);
    Y = abcd2y(A);

    % Parameters of the via alone, without the feeding lines
    Av = inv(Al1)*A*inv(Al2);
    Yv = abcd2y(Av);

    % Calculate via L and C, T-network equivalent
    Z = inv(Yv);
    Z1 = Z(1,1) - Z(2,1); % half-series
    Z2 = Z(2,1);          % shunt
    Z3 = Z(2,2) - Z(2,1); % half-series

    Ct = j*imag(1/Z2)/(j*freq)
    Lt1 = j*imag(Z1)/(j*freq)
    Lt2 = j*imag(Z3)/(j*freq)

    % Add y-parameters for this frequency to the overall results
    Yf = cat(3, Yf, Y);
    Yndf = cat(3, Yndf, Y1);
    Yl1f = cat(3, Yl1f, Yl1);
    Yl2f = cat(3, Yl2f, Yl2);
    Yvf = cat(3, Yvf, Yv);

end

tswrite( 'via2_64.y2p',      freqs/(2*pi), Yf   )
tswrite( 'via2_64_nd.y2p',   freqs/(2*pi), Yndf )
tswrite( 'via2_64_line1.y2p', freqs/(2*pi), Yl1f  )
tswrite( 'via2_64_line2.y2p', freqs/(2*pi), Yl2f  )
tswrite( 'via2_64_via.y2p',  freqs/(2*pi), Yvf  )

%% % For the drawing
%% wg = wgparams(1e9, a, b, h, nx, ny);

%% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, zeros(50000,1));
%% %trimesh(Tri, X, Y, Z);
%% trisurf(Tri, X, Y, Z, C);
%% xlim([ -dx  a+dx ])
%% ylim([ -dy  b+dy ])
%% zlim([ 0    sum(lh)*1.2    ])

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
