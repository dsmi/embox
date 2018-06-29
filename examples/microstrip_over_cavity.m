% x-directed transmission line of nonzero thickness represented by
% by two or more metal layers, deembedded by simulating line of length L
% and then 2L

addpath(genpath([ pwd, '/..' ]));
% nodal solver for the currents plot
%% addpath('d:/octave/nodal/nodal');


mil2meter = 2.54e-5;

% Dimensions
lnpar.w = 8.0*mil2meter;     % microstrip width
lnpar.l = 240.0*mil2meter;   % microstrip length
lnpar.h5 = 0.675*mil2meter;  % line copper
lnpar.h4 = 5.00*mil2meter;   % dielectric between cavity and line
lnpar.h3 = 1.35*mil2meter;   % cavity upper copper
lnpar.h2 = 5.00*mil2meter;   % cavity dielectric
lnpar.h1 = 1.35*mil2meter;   % bottom copper
lnpar.cv = 1;                % 0 -- no cavity, 1 -- cavity, 2 -- grounded
lnpar.cw = 80.0*mil2meter;   % cavity width and length 

% Mesh options
lnpar.nl = 4;       % number of the metal layers in the microstrip mesh
lnpar.nx = 240/2;     % cells along x -- line is x-directed
lnpar.ny = 240/2;     % cells along y
lnpar.a  = lnpar.l; % x-size of the waveguide
lnpar.b  = lnpar.l; % y-size of the waveguide

function [ wg mesh ports portw ] = mymkmesh(freq, lnpar)

    % Layers stack
    nls = lnpar.nl - 1; % number of 'steps' between metal layers
    h5 = lnpar.h5;
    h = [ lnpar.h1 lnpar.h2 lnpar.h3 lnpar.h4 repmat(h5/nls, 1, nls) ];

    % Dielectric/metal permittivities
    epsc = eps0 - j*ccopper/freq; % permittivity for copper
    epsd = eps0 * debye2(4.3, 0.02, 1e9*2*pi, freq); % dielectric
    weps = [ epsc epsd epsd epsd repmat(eps0, 1, nls) ];

    % parametes of the enclosure to pass to mkzmat
    nx      = lnpar.nx;
    ny      = lnpar.ny;
    wg      = wgparams(freq, lnpar.a, lnpar.b, h, nx, ny);
    wg.weps = weps; 
    wg.Ggr0 = 0; % no top ground
    wg.Gls0 = 0; % no bottom ground (we have layer of copper)
    wg.cnx  = 4; % to speed up things
    wg.cny  = 4;

    % Mesh cell size
    dx = wg.a/lnpar.nx;
    dy = wg.b/lnpar.ny;

    % Bitmaps
    B0 = zeros(nx+2, ny+2); % canvas
    BL = drawline(B0, 0, (ny+2)/2, nx+2, (ny+2)/2, lnpar.w/dy); % line
    BV = traceedges(BL); % line vias

    % Cavity bitmaps
    xc = (nx+2)/2;
    yc = (ny+2)/2;
    chp = lnpar.cw/dx;
    cwp = lnpar.cw/dy;
    BC = drawline(B0, xc - chp/2, yc, xc + chp/2, yc, cwp); % line
    CV = traceedges(BC); % line vias

    % Trace mesh
    for li=1:lnpar.nl
        mesh.layers(end + 1) = mklayer(BL, BV*( li > 1 ), li + 3, ccopper);
    end

    if 1 == lnpar.cv
        mesh.layers(end + 1) = mklayer(BC, 0*CV, 2, ccopper);
        mesh.layers(end + 1) = mklayer(BC, CV, 3, ccopper);
    elseif 2 == lnpar.cv
        mesh.layers(end + 1) = mklayer(0*BC, CV, 2, ccopper);
        mesh.layers(end + 1) = mklayer(BC, CV, 3, ccopper);
    end
    
    % Identify ports
    b1 = findbases(mesh, nx, ny, 0, 0, 0, 1);
    b2 = findbases(mesh, nx, ny, 1, 0, 1, 1);

    ports = { b1' b2'};
    portw = { b1'*0-dy b2'*0+dy };

end

function [ Y I wg mesh ports portw Z ] = simline(freq, lnpar)
    [ wg mesh ports portw ] = mymkmesh(freq, lnpar);
    [ Y I Z ]=solvey(wg, mesh, ports, portw);
end

% Obtain port discontinuity matrix via simulation of the calibration
% standards of single and double lengths
function [ Ad R Y1 Y2 Yd wg mesh ports portw ] = calcporta(freq, lnpar)

    % signle-l
    lnpar1 = lnpar;
    lnpar1.a  = lnpar.a/2;
    lnpar1.nx = lnpar.nx/2;
    lnpar1.cv = 0; % no cavity
    Y1 = simline(freq, lnpar1);

    % double-l
    lnpar2 = lnpar;
    lnpar2.cv = 0; % no cavity
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
freqs = linspace(1e7, 2e10, 100)*2*pi;
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
    
    % simulate line
    [ Ynd I wg mesh ports portw Z ] = simline(freq, lnpar);

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

    resultsName = 'microstrip_3d_temp';
    writesp([ resultsName, '.s2p' ], freqs, Yf);
    writesp([ resultsName, '_nd.s2p' ], freqs, Yndf);
    writesp([ resultsName, '_l1.s2p' ], freqs, Y1f);
    writesp([ resultsName, '_l2.s2p' ], freqs, Y2f);

end

%% [ wg mesh ] = mymkmesh(1e9, lnpar);
%% I = zeros(20000,1);

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

%% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I(:,1));
%% %graphics_toolkit('fltk')
%% trisurf(Tri, X, Y, Z, C);
%% shading interp
%% %trimesh(Tri, X, Y, Z);
%% xlim([ -(wg.a/wg.nx)  wg.a+2*(wg.a/wg.nx) ])
%% ylim([ -(wg.b/wg.ny)  wg.b+2*(wg.b/wg.ny) ])
%% %zlim([ -(wg.b/wg.ny)  wg.b+2*(wg.b/wg.ny) ])
%% %zlim([ 0 ( lnpar.h1 + lnpar.h2 + lnpar.h3 + lnpar.h4 + lnpar.h5 ) ])
%% zlim([ 0 lnpar.l/4 ])
