% Via through cavity decomposition 
% Feed line simulation

addpath(genpath([ pwd, '/..' ]));

mil2meter = 2.54e-5;

% Dimensions
lnpar.lw = 8.0*mil2meter;     % microstrip width
lnpar.l  = 120.0*mil2meter;   % length of the simulation area
lnpar.w  = 120.0*mil2meter;   % width of the simulation area
lnpar.cw = 80.0*mil2meter;    % cavity width and length 
lnpar.hr = 10.0*mil2meter;    % via hole radius

% Stackup
lnpar.h = 0.675*mil2meter; % thickness
lnpar.d = (5.00)*mil2meter;

% Mesh options
lnpar.nx = 120/2;      % cells along x
lnpar.ny = 120/2;      % cells along y
lnpar.a  = lnpar.l; % x-size of the waveguide
lnpar.b  = lnpar.w; % y-size of the waveguide


% line above ground
% d  -- line-to-ground
% ld -- length divisor
function [ wg mesh ports portw ] = mklinemesh(freq, lnpar, d, ld)

    % Layers stack
    h = [ d d lnpar.h(end) ]; % ground, space, line

    % Dielectric/metal permittivities
    epsc = eps0 - j*ccopper/freq; % permittivity for copper
    epsd = eps0 * debye2(4.3, 0.02, 1e9*2*pi, freq); % dielectric
    weps = [ epsc epsd eps0 ];

    % parametes of the enclosure to pass to mkzmat
    nx      = lnpar.nx/ld;
    ny      = lnpar.ny;
    wg      = wgparams(freq, lnpar.a/ld, lnpar.b, h, nx, ny);
    wg.weps = weps; 
    wg.Ggr0 = 0; % no top ground
    wg.Gls0 = 0; % no bottom ground
    wg.cnx  = 4; % to speed up things
    wg.cny  = 4;

    % Mesh cell size
    dx = wg.a/lnpar.nx;
    dy = wg.b/lnpar.ny;

    % Bitmaps
    B0 = zeros(nx+2, ny+2); % canvas
    xc = (nx+2)/2;
    yc = (ny+2)/2;

    % Bottom line
    BL = drawline(B0, 0, yc, nx+2, yc, lnpar.lw/dy);
    BV = traceedges(BL); % line vias

    % Make mesh from the layers
    mesh.layers = struct( [] );
    mesh.layers(end + 1) = mklayer(BL, BV*0, 2, ccopper); % top line
    mesh.layers(end + 1) = mklayer(BL, BV, 3, ccopper);   % top line

    % Merge layers one the same positions
    mesh = mergelayers(mesh);

    % Identify ports
    b1 = findbases(mesh, nx, ny, 0, 0, 0, 1);
    b2 = findbases(mesh, nx, ny, 1, 0, 1, 1);

    ports = { b1' b2'};
    portw = { b1'*0-dy b2'*0+dy };

end

% Obtain port discontinuity matrix via simulation of the calibration
% standards of single and double lengths
function [ Ad R Y1 Y2 Yd wg mesh ports portw ] = calcporta(freq, lnpar)

    % signle-l
    [ wg mesh ports portw ] = mklinemesh(freq, lnpar, lnpar.d, 4);
    Y1 = solvey(wg, mesh, ports, portw);

    % double-l
    [ wg mesh ports portw ] = mklinemesh(freq, lnpar, lnpar.d, 2);
    Y2 = solvey(wg, mesh, ports, portw);

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
freqs = linspace(1e7, 4e10, 120)*2*pi;

% Simulation results
Y1f = [];
Y2f = [];
Yndf = [];
Yf = [];

for ifr = 1:length(freqs)

    freq = freqs(ifr);
    wavelen = 2*pi/(freq * sqrt(4.3 * eps0 * mu0));
    fprintf( 'Solving for frequency %.8e, wlen %.8e, step %i of %i\n', freq, wavelen, ifr, length(freqs) );

    % deembedding simulations
    [ Ad R Y1 Y2 Yd wg mesh ports portw ] = calcporta(freq, lnpar);

    fprintf( 'De-embdeeing accuracy %.8e\n', sum(sum(R)) );
    
    % This is a line simulation!
    Ynd = Y1;

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

    resultsName = 'via_through_cavity_line_5_30';
    writesp([ resultsName, '.s2p' ], freqs, Yf);
    writesp([ resultsName, '_nd.s2p' ], freqs, Yndf);
    writesp([ resultsName, '_l1.s2p' ], freqs, Y1f);
    writesp([ resultsName, '_l2.s2p' ], freqs, Y2f);

end

