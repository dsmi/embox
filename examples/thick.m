% x-directed transmission line of nonzero thickness represented by
% by a pair of layers, deembedded by simulating line of length L
% and then 2L

addpath(genpath([ pwd, '/..' ]));

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

% Dimensions
th = 3e-5;    % trace thickness
tw = 1.22e-4; % trace width
d1 = 1e-4;    % trace-to-plane separation
a  = tw*16;   % x-size of the enclosure/waveguide
b  = tw*16;   % y-size of the enclosure/waveguide

% dielectric
er = 4.05;  

% Mesh options
nx = 96; % cells along x
ny = 96; % cells along y
nl = 3;  % number of layers in the trace

% Mesh cell size
dx = a/nx;
dy = b/ny;

%% % Create mesh for plotting
%% [ mesh h epsl ] = mklinemesh( 1e9, nx, ny, a, b, th, tw, d1, nl, er );

%% % For the plotting
%% wg = wgparams( 1e9, a, b, h, nx, ny );

%% [ Tri, X, Y, Z, C ] = mesh2tri( wg, mesh, zeros( 50000,1 ) );
%% trimesh(Tri, X, Y, Z);
%% %% trisurf(Tri, X, Y, Z, C);
%% xlim([ -dx  a+dx ])
%% ylim([ -dy  b+dy ])
%% zlim([ 0   a ])

% angular frequencies
freqs = linspace(1e6,2e10,10)*2*pi;
%% freqs = 1e8;

% Simulation results
Yf = [];
Yndf = [];

function [ D, Y1, Y2 ] = deembsims( fsim1, fsim2 )

    fprintf('Running deembedding-l simulation...\n')
    Y1 = fsim1( );

    fprintf('Running deembedding-l*2 simulation...\n')
    Y2 = fsim2( );

    A1 = y2abcd(Y1);
    A2 = y2abcd(Y2);

    % Obtain the port discontinuity
    [ D dtol ] = splitd( A1*inv(A2)*A1 );
    fprintf('Deembedding accuracy is %.8e\n', dtol)
    
end

% angular frequency
for freq = freqs

    fprintf( 'Simulation %i of %i...\n', find(freq == freqs), length(freqs) )

    % Obtain the port discontinuity
    fsim1 = @( ) simline( freq, nx/8, ny, a/8, b, th, tw, d1, nl, er );
    fsim2 = @( ) simline( freq, nx/4, ny, a/4, b, th, tw, d1, nl, er );
    D = deembsims( fsim1, fsim2 );

    fprintf('Running the entire line simulation...\n')
    Y1 = simline( freq, nx, ny, a, b, th, tw, d1, nl, er );

    % De-embedded Y1
    invD = inv(D);
    A = invD*y2abcd(Y1)*invD;
    Y = abcd2y(A);

    % Add y-parameters for this frequency to the overall results
    Yf = cat(3, Yf, Y);
    Yndf = cat(3, Yndf, Y1);

    % Characteristic impedance and delay from simulation
    N=size(A,1)/2; % number of ports on each side
    A12=A(1:N,N+1:end);
    A21=A(N+1:end,1:N);
    Z0s=sqrt(A12*inv(A21))
    tds=acos(A(1,1))./freq

end

tswrite( 'tline_96.y2p', freqs/(2*pi), Yf )
tswrite( 'tline_96_nd.y2p', freqs/(2*pi), Yndf )

%% wavelen = 2*pi/(freq * sqrt(eps0 * mu0))

%% l=a;   % length

%% C = 6.2415337e-011;
%% L = 1.7826551e-007;
%% Ro = 2.7583998e-001;
%% Rs = 2.1258788e-004;
%% %% conductivity = 5.8e7; % Copper
%% %% skin_depth = sqrt(2./(freq*mu0*conductivity));
%% %% Rs = 1./(w*skin_depth*conductivity);
%% Zo = j*freq*L + Ro + Rs*(1+j)*sqrt(freq/(2*pi)); % impedance per len
%% Yo = j*freq*C; % admittance per len
%% Z0 = sqrt(Zo/Yo)

%% gamma=sqrt(Zo*Yo); % propagation constant
%% vp=j*freq/gamma;
%% td=l/vp;

%% Z11=Z0/tanh(gamma.*l);
%% Z12=Z0/sinh(gamma.*l);
%% Ztl=[ Z11 Z12 ; Z12 Z11 ];
%% Ytl=inv(Ztl)
