% x-directed stripline line of nonzero thickness
% Deembedded by simulating line of length L
% and then 2L

addpath(genpath([ pwd, '/..' ]));

mil2meter = 2.54e-5;

% Dimensions
lnpar.w = 8.0*mil2meter;   % stripline width
lnpar.l = 200.0*mil2meter; % stripline length
lnpar.t = 0.675*mil2meter; % stripline thickness
lnpar.d1 = 5*mil2meter;  % line to lower plane
lnpar.d2 = 5*mil2meter;  % line to upper plane

% Mesh options
lnpar.nl = 4;       % number of the metal layers in the stripline mesh
lnpar.nx = 100;     % cells along x
lnpar.ny = 50;      % cells along y
lnpar.a  = lnpar.l;        % x-size of the waveguide
lnpar.b  = (lnpar.w*25/4); % y-size of the waveguide

function [ Y I wg mesh ports portw ] = simline(freq, lnpar)

    % Layers stack
    nls = lnpar.nl - 1; % number of 'steps' between metal layers
    h = [ lnpar.d1 lnpar.d1 repmat(lnpar.t/nls, 1, nls) lnpar.d2 lnpar.d2 ];

    % Layer of copper at the bottom is the ground - it is not an ideal conductor
    epsc = eps0 - j*ccopper/freq; % permittivity for copper
    epsd = eps0 * debye2(4.3, 0.02, 1e9*2*pi, freq); % dielectric
    weps = [ epsc epsd repmat(epsd, 1, nls) epsd epsc ];

    % parametes of the enclosure to pass to mkzmat
    nx      = lnpar.nx;
    ny      = lnpar.ny;
    wg      = wgparams(freq, lnpar.a, lnpar.b, h, nx, ny);
    wg.weps = weps; 
    wg.Ggr0 = 0; % no top ground (we have layer of copper)
    wg.Gls0 = 0; % no bottom ground (we have layer of copper)
    wg.cnx  = 4; % to speed up things
    wg.cny  = 4;

    % Mesh cell size
    dy = wg.b/lnpar.ny;

    w = lnpar.w;
    nl = lnpar.nl;
    l = lnpar.l;

    % Make mesh from the layers
    mesh.layers = struct( [] );

    % trace mesh
    for li=1:lnpar.nl
        % x-directed trace
        B0 = zeros(nx+2, ny+2);
        BL = drawline(B0, 0, (ny+2)/2, nx+2, (ny+2)/2, w/dy);
        BV = B0;
        if li > 1
            BV = traceedges(BL); % vias
        end
        if li > 1 && li < nl
            BL = traceedges(BL); % only edges at intermediate layers
        end
        mesh.layers(end + 1) = mklayer(BL, BV, li + 1, ccopper);
    end
    
    % Identify ports
    b1 = findbases(mesh, nx, ny, 0, 0, 0, 1);
    b2 = findbases(mesh, nx, ny, 1, 0, 1, 1);

    ports = { b1' b2'};
    portw = { b1'*0-dy b2'*0+dy };
    [ Y I ]=solvey(wg, mesh, ports, portw);
end

% Obtain port discontinuity matrix via simulation of the calibration
% standards of single and double lengths
function [ Ad R Y1 Y2 Yd wg mesh ports portw ] = calcporta(freq, lnpar)

    % signle-l
    lnpar1 = lnpar;
    lnpar1.a  = lnpar.a/2;
    lnpar1.nx = lnpar.nx/2;
    Y1 = simline(freq, lnpar1);

    % double-l
    lnpar2 = lnpar;
    [ Y2 I2 wg mesh ports portw ] = simline(freq, lnpar2);

    A1=y2abcd(Y1);
    A2=y2abcd(Y2);

    Add = A1*inv(A2)*A1; % double port discontinuity
    N = size(Add,1)/2;  % number of ports on each side
    A = Add(1:N,1:N);
    B = Add(1:N,N+1:end);
    C = Add(N+1:end,1:N);
    D = Add(N+1:end,N+1:end);

    Yd = C*0.5; % port discontinuity admittance
    Ad  = [ A B ; Yd D ]; 
    
    % Estimate accuracy of deembedding.
    % Expected values of the shunt ABCD
    % are: A = 1; B = 0; C = 1/Zs; D = 1.0
    R = abs( A - ones(N, N) ) + abs( B ) + abs( D - ones(N, N) );
end

function writesp(fileName, freq, Y, t, r0)
  tswrite(fileName, freq/(2*pi), renorms(y2s(Y), [ 1 1 ], [ 50 50 ]), 'S', 50)
end

% angular frequencies
freqs = linspace(1e7, 4e10, 100)*2*pi;
%freqs = 3e10*2*pi;

% Simulation results
Y1f = [];
Y2f = [];
Yndf = [];
Yf = [];

for ifr = 1:length(freqs)

    %% step = find(freqs == freq);
    %% steps = length(freqs);
    freq = freqs(ifr);
    wavelen = 2*pi/(freq * sqrt(eps0 * mu0))
    fprintf( 'Solving for frequency %.8e, step %i of %i\n', freq, ifr, length(freqs) );

    % deembedding simulations
    [ Ad R Y1 Y2 Yd wg mesh ports portw ] = calcporta(freq, lnpar);

    fprintf( 'De-embdeeing accuracy %.8e\n', sum(sum(R)) );
    
    % Non-deembedded Y (we reuse the deembedding simulation result)
    Ynd = Y2;

    And = y2abcd(Ynd);

    invAd = inv(Ad);
    A = invAd * And * invAd;

    Y = abcd2y(A);

    % Add y-parameters for this frequency to the overall results
    Yf = cat(3, Yf, Y);

    % Capture y-parameters of the de-embedding standards and no-deembedded results
    Yndf = cat(3, Yndf, Ynd);
    Y1f = cat(3, Y1f, Y1);
    Y2f = cat(3, Y2f, Y2);

    % Characteristic impedance and delay from simulation
    N=size(A,1)/2; % number of ports on each side
    A12=A(1:N,N+1:end);
    A21=A(N+1:end,1:N);
    Z0s=sqrt(A12*inv(A21))
    tds=acos(A(1,1))./freq

    resultsName = 'stripline_3d';
    writesp([ resultsName, '.s2p' ], freqs, Yf);
    writesp([ resultsName, '_nd.s2p' ], freqs, Yndf);
    writesp([ resultsName, '_l1.s2p' ], freqs, Y1f);
    writesp([ resultsName, '_l2.s2p' ], freqs, Y2f);

end

%% I = zeros(20000, 1);

%% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I(:,1));
%% %trisurf(Tri, X, Y, Z, C);
%% trimesh(Tri, X, Y, Z);
%% xlim([ -(wg.a/wg.nx)  wg.a+2*(wg.a/wg.nx) ])
%% ylim([ -(wg.b/wg.ny)  wg.b+2*(wg.b/wg.ny) ])
%% zlim([ 0 (lnpar.d1*2 + lnpar.t + lnpar.d2*2) ])


