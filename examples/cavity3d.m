%
% Rectangular cavity with two ports on the narrow side, to be combined
% with the stripline models -- power aware si experiments.
%


addpath(genpath([ pwd, '/..' ]));

% The geometry
inch2meter = 2.54e-2;
mil2meter = 1.0e-3*inch2meter;
d = (5.0 + 0.675 + 5.0)*mil2meter; % plane-to-plane separation
w = 100.0*mil2meter;  % cavity width
l = 200.0*mil2meter;  % cavity length
lw = 4.0*mil2meter;   % stripline (and therefore the port) width
fw = 4.0*mil2meter;   % feed lines width

% Mesh settings
cny = 50;
cnx = cny*2;
xo = 12; % edge-to-enclosure doubled
yo = 12;

% Mesh cell size
dx = l/cnx;
dy = w/cny;

% Dielectric params
lt = 0.02;
er0 = 4.3;
fr = 1e9;

% angular frequencies
%freqs = 1e9*2*pi;
freqs = linspace(1e7, 4e10, 100)*2*pi;

function Y = simfeed(freq, h, weps, nx, ny, dx, dy, fw)

    % parametes of the enclosure to pass to mkzmat
    a       = dx * nx;
    b       = dy * ny;
    wg      = wgparams(freq, a, b, h, nx, ny);
    wg.weps = weps; 
    wg.Ggr0 = 0; % no top ground
    wg.Gls0 = 0; % no bottom ground
    wg.cnx  = 4; % to speed up things
    wg.cny  = 4;

    % Canvas
    B0 = zeros(nx+2, ny+2);
    BL = drawline(B0, 0, (ny+2)/2, nx+2, (ny+2)/2, fw/dy); % line

    % Make mesh from the layers
    mesh.layers = struct( [] );
    mesh.layers(end + 1) = mklayer(BL, 0*BL, 1, ccopper); % Lower wall feed
    mesh.layers(end + 1) = mklayer(BL, 0*BL, 2, ccopper); % Upper wall feed

    % Identify ports
    b1 = findbases(mesh, nx, ny, 0, 0, 0, 1, @(l) l == 1);
    b2 = findbases(mesh, nx, ny, 1, 0, 1, 1, @(l) l == 1);

    ports = { b1' b2' };
    portw = { b1'*0-dy b2'*0+dy };
    [ Y I ] = solvey(wg, mesh, ports, portw);
end

% Feed line ABCD matrix and feed-to-wall discontinuity
function [ Ad A1 R ] = calcfeeda(freq, h, weps, dx, dy, cny, xo, yo, fw)

    % double-l
    nx = xo;
    ny = cny + yo;
    Y2 = simfeed(freq, h, weps, nx, ny, dx, dy, fw);

    % signle-l
    nx = xo / 2;
    ny = cny + yo;
    Y1 = simfeed(freq, h, weps, nx, ny, dx, dy, fw);

    A1 = y2abcd(Y1);
    A2 = y2abcd(Y2);

    Add = A1*inv(A2)*A1; % double port discontinuity
    N = size(Add,1)/2;  % number of ports on each side
    A = Add(1:N,1:N);
    B = Add(1:N,N+1:end);
    C = Add(N+1:end,1:N);
    D = Add(N+1:end,N+1:end);

    Yd = C*0.5; % port discontinuity admittance
    Ad  = [ A B ; Yd D ];
    %Ad = splitd(A1*inv(A2)*A1);

    
    % Estimate accuracy of deembedding.
    % Expected values of the shunt ABCD
    % are: A = 1; B = 0; C = 1/Zs; D = 1.0
    R = abs( A - ones(N, N) ) + abs( B ) + abs( D - ones(N, N) );
end


% Simulation results
Yf = [];
Ydf = [];
Ylf = [];

for ifr = 1:length(freqs)

    freq = freqs(ifr);
    wavelen = 2*pi/(freq * sqrt(eps0 * mu0));
    fprintf( 'Solving for frequency %.8e, step %i of %i, wlen %.8e\n', freq, ifr, length(freqs), wavelen );

    % Layers stack
    h    = [ d d d ];
    er   = debye(er0, lt, fr, freq/(2*pi));
    weps = [ eps0 eps0*er eps0 ];

    % Feed line params for deembedding
    [ Ad A1 R ] = calcfeeda(freq, h, weps, dx, dy, cny, xo, yo, fw);

    % parametes of the enclosure to pass to mkzmat
    nx      = cnx + xo;
    ny      = cny + yo;
    a       = dx * nx;
    b       = dy * ny;
    wg      = wgparams(freq, a, b, h, nx, ny);
    wg.weps = weps; 
    wg.Ggr0 = 0; % no top ground
    wg.Gls0 = 0; % no bottom ground
    wg.cnx  = 4; % to speed up things
    wg.cny  = 4;

    % Canvas
    B0 = zeros(nx+2, ny+2);

    % Cavity upper/lower wall
    lxb = 1.5 + xo/2; % beginning of the cavity and line
    lxe = nx + 0.5 - xo/2; % end of the cavity and line
    yc = (ny + 2)/2;
    B1 = drawline(B0, lxb, yc, lxe, yc, cny);
    B2 = drawline(B1, 0, (ny+2)/2, nx+2, (ny+2)/2, fw/dy); % feed line

    % Make mesh from the layers
    mesh.layers = struct( [] );
    mesh.layers(end + 1) = mklayer(B2, 0*B2, 1, ccopper); % Lower wall + feed
    mesh.layers(end + 1) = mklayer(B2, 0*B2, 2, ccopper); % Upper wall + feed

    % Identify ports
    b1 = findbases(mesh, nx, ny, 0, 0, 0, 1, @(l) l == 1);
    b2 = findbases(mesh, nx, ny, 1, 0, 1, 1, @(l) l == 1);

    ports = { b1' b2' };
    portw = { b1'*0-dy b2'*0+dy };
    [ Y I ] = solvey(wg, mesh, ports, portw);

    % De-embedding
    invAd = inv(Ad); % inverted discontinuity
    invF = inv(invAd * A1 * invAd); % inverted feed

    Yd = abcd2y( invF * invAd * y2abcd(Y) * invAd * invF );

    Yf = cat(3, Yf, Y);
    Ydf = cat(3, Ydf, Yd);
    Ylf = cat(3, Ylf, abcd2y(inv(invF)));

end

resultsName = 'cavity_3d_2port';
tswrite([ resultsName, '.y2p' ], freqs/(2*pi), Yf, 'Y', 50);
tswrite([ resultsName, '_deemb.y2p' ], freqs/(2*pi), Ydf, 'Y', 50);
tswrite([ resultsName, '_feed.y2p' ], freqs/(2*pi), Ylf, 'Y', 50);

%% I = zeros(100000, 1);
%% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I(:,1));
%% %trisurf(Tri, X, Y, Z, C);
%% trimesh(Tri, X, Y, Z);
%% xlim([ -(wg.a/wg.nx)  wg.a+2*(wg.a/wg.nx) ])
%% ylim([ -(wg.b/wg.ny)  wg.b+2*(wg.b/wg.ny) ])
%% zlim([ 0 3*d ])

