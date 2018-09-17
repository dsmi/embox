% Via through cavity

addpath(genpath([ pwd, '/..' ]));
% nodal solver for the currents plot
addpath('d:/octave/nodal/nodal');


mil2meter = 2.54e-5;

% Dimensions
lnpar.lw = 8.0*mil2meter;     % microstrip width
lnpar.l  = 120.0*mil2meter;   % length of the simulation area
lnpar.w  = 120.0*mil2meter;   % width of the simulation area
lnpar.cw = 80.0*mil2meter;    % cavity width and length 
lnpar.hr = 10.0*mil2meter;    % via hole radius

% Stackup
lnpar.h(10) = 0.675*mil2meter;  % upper line
lnpar.h(9) = 5.00*mil2meter;    % between cavity and upper line
lnpar.h(8) = 1.35*mil2meter;    % cavity upper copper
lnpar.h(7) = 5.00*mil2meter;   % dielectric, in-cavity
lnpar.h(6) = 1.35*mil2meter;  % plane
lnpar.h(5) = 5.00*mil2meter;  % dielectric
lnpar.h(4) = 1.35*mil2meter;  % dielectric
lnpar.h(3) = 5.00*mil2meter;  % dielectric
lnpar.h(2) = 0.675*mil2meter; % bottom line
lnpar.h(1) = 1.00*mil2meter;  % empty layer at the bottom

% Mesh options
lnpar.nx = 120/2;      % cells along x
lnpar.ny = 120/2;      % cells along y
lnpar.a  = lnpar.l; % x-size of the waveguide
lnpar.b  = lnpar.w; % y-size of the waveguide

function [ wg mesh ports portw ] = mkviamesh(freq, lnpar)

    % Layers stack
    h = lnpar.h;

    % Dielectric/metal permittivities
    epsc = eps0 - j*ccopper/freq; % permittivity for copper
    epsd = eps0 * debye2(4.3, 0.02, 1e9*2*pi, freq); % dielectric
    weps = [ eps0 eps0 repmat(epsd, 1, numel(h) - 3) eps0 ];

    % parametes of the enclosure to pass to mkzmat
    nx      = lnpar.nx;
    ny      = lnpar.ny;
    wg      = wgparams(freq, lnpar.a, lnpar.b, h, nx, ny);
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

    % Via
    BC = drawcir(B0, xc, yc, lnpar.lw/dy/2);
    BV = traceedges(BC);
    BH = drawcir(B0, xc, yc, lnpar.hr/dy);

    % Bottom line
    BL0 = drawline(BC, 0, yc, xc, yc, lnpar.lw/dy);
    BV0 = traceedges(BL0); % line vias

    % Top line
    BL1 = drawline(BC, xc, yc, nx+2, yc, lnpar.lw/dy);
    BV1 = traceedges(BL1); % line vias

    % Ground -- filled with hole
    GR = B0 + 1 - BH;
    GV = traceedges(GR); % ground vias

    % Cavity
    chp = lnpar.cw/dx;
    cwp = lnpar.cw/dy;
    CA = drawline(B0, xc - chp/2, yc, xc + chp/2, yc, cwp) - BH;
    CV = traceedges(CA); % cavity vias

    % Make mesh from the layers
    mesh.layers = struct( [] );
    mesh.layers(end + 1) = mklayer(BL0, BV0*0, 1, ccopper); % Bottom line
    mesh.layers(end + 1) = mklayer(BL0, BV0, 2, ccopper);   % Bottom line
    mesh.layers(end + 1) = mklayer(BV*0, BV, 3, ccopper);   % via
    mesh.layers(end + 1) = mklayer(BV*0, BV, 4, ccopper);   % via
    mesh.layers(end + 1) = mklayer(BV*0, BV, 5, ccopper);   % via
    mesh.layers(end + 1) = mklayer(BV*0, BV, 6, ccopper);   % via
    mesh.layers(end + 1) = mklayer(BV*0, BV, 7, ccopper);   % via
    mesh.layers(end + 1) = mklayer(BV*0, BV, 8, ccopper);   % via
    mesh.layers(end + 1) = mklayer(BV*0, BV, 9, ccopper);   % via
    mesh.layers(end + 1) = mklayer(BL1, BV1*0, 9, ccopper); % top line
    mesh.layers(end + 1) = mklayer(BL1, BV1, 10, ccopper);   % top line
    mesh.layers(end + 1) = mklayer(GR, GV*0, 5, ccopper);   % ground
    mesh.layers(end + 1) = mklayer(GR, GV, 6, ccopper);   % ground
    mesh.layers(end + 1) = mklayer(CA, CV*0, 7, ccopper);   % cavity
    mesh.layers(end + 1) = mklayer(CA, CV, 8, ccopper);   % cavity

    % Merge layers one the same positions
    mesh = mergelayers(mesh);

    % Identify ports
    b1 = findbases(mesh, nx, ny, 0, 0, 0, 1, @(l) l <= 2);
    b2 = findbases(mesh, nx, ny, 1, 0, 1, 1, @(l) l >= 9);

    ports = { b1' b2'};
    portw = { b1'*0-dy b2'*0+dy };

end

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
    [ wg mesh ports portw ] = mklinemesh(freq, lnpar, sum(lnpar.h(3:5)), 2);
    Y1 = solvey(wg, mesh, ports, portw);

    % double-l
    [ wg mesh ports portw ] = mklinemesh(freq, lnpar, sum(lnpar.h(3:5)), 1);
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
%freqs = 1e8*2*pi;

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
    
    % simulate the via
    [ wg mesh ports portw ] = mkviamesh(freq, lnpar);
    Ynd = solvey(wg, mesh, ports, portw);
    %% Ynd = Y2;

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

    resultsName = 'via_through_cavity_3d';
    writesp([ resultsName, '.s2p' ], freqs, Yf);
    writesp([ resultsName, '_nd.s2p' ], freqs, Yndf);
    writesp([ resultsName, '_l1.s2p' ], freqs, Y1f);
    writesp([ resultsName, '_l2.s2p' ], freqs, Y2f);

end

%% [ wg mesh ] = mklinemesh(1e9, lnpar, sum(lnpar.h(3:5)), 1);
%% I = zeros(25000,1);

%% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I(:,1));
%% graphics_toolkit('fltk')
%% %% trisurf(Tri, X, Y, Z, C);
%% %% shading interp
%% trimesh(Tri, X, Y, Z);
%% xlim([ -(wg.a/wg.nx)  wg.a+2*(wg.a/wg.nx) ])
%% ylim([ -(wg.b/wg.ny)  wg.b+2*(wg.b/wg.ny) ])
%% %zlim([ -(wg.b/wg.ny)  wg.b+2*(wg.b/wg.ny) ])
%% %zlim([ 0 ( lnpar.h1 + lnpar.h2 + lnpar.h3 + lnpar.h4 + lnpar.h5 ) ])
%% zlim([ 0 lnpar.l/4 ])


%% % Calculate source voltages for the current plot
%% %Yt = 1/Z0s
%% Yt = 1/81.3 % termination admittance for de-embedded DUT
%% Yd
%% Ys = Yt - Yd % non-deembeded source and termination admittance
%% branches=[ 2 1 ; 3 1 ; 2 1 ; 3 1 ];
%% W = [ 0 ; 0 ; 1 ; 0 ];
%% K = 0*W;
%% Yc = zeros(size(branches, 1), size(branches, 1));
%% Yc(1:2,1:2) = Y;
%% Yc(3,3) = Ys;
%% Yc(4,4) = Ys;
%% [ F, V ] = solve(branches,Yc,W,K);

%% % excitation voltage vector
%% XV = zeros(size(Z,1), 1);
%% XV(ports{1}) = V(1);
%% XV(ports{2}) = V(2);

%% I = (Z\XV);

