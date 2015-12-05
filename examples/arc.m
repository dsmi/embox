% arc

addpath(genpath([ pwd, '/../..' ]));

function [ Y I Z wg mesh ports ] = simulate(freq, h, t, w, d, len, ac, ar, doarc)

    % Simulation settings
    nl = 2;   % Number of metal layers in the microstrip model
    nx = 64;  % cells along x
    a2w = 32; % a/w, wavegiude width to trace width ratio
    ny = len/(w*a2w/nx); % to get dy=dx
    vi = 1;  % 1 - use vias, 0 - no vias

    % parametes of the enclosure to pass to mkzmat
    a = w*a2w; % x-size of the waveguide
    b = len;   % y-size of the waveguide
    epsc = eps0 - j*ccopper/freq; % permittivity for copper - incorporates conductivity
    epsd = eps0*debye(4.05, 0.02, 1e9, freq/(2*pi)); % dielectric
    weps = [ epsc  epsd  repmat(eps0,     1, nl-1)  eps0 ]; % layer eps
    wh   = [ h     h     repmat(t/(nl-1), 1, nl-1)  h    ]; % layer h

    % Mesh cell size
    dx = a/nx;
    dy = b/ny;

    wg = wgparams(freq, a, b, wh, nx, ny);
    wg.weps = weps; 
    wg.Gls0 = 0; % no bottom ground
    wg.Ggr0 = 0; % no top ground
    wg.cnx = 4;  % to make the things a little faster
    wg.cny = 4;

    % Arc intermediate values
    acn = ac/a;               % arc center in normalized coordinates
    arn = ar/a;               % arc radius in normalized coordinates
    sep = (d+w)/a;            % center-to-center separation
    ct1 = acn + arn - sep/2;  % x- and y-centers of the straight traces
    ct2 = acn + arn + sep/2;

    if doarc
        % Traces bitmap - arc
        B = zeros(nx+2, ny+2); 
        B = linefromto(B, ct1, -1.0, ct1, acn, w/a);
        B = linefromto(B, ct2, -1.0, ct2, acn, w/a);
        B = linefromto(B, -1.0, ct1, acn, ct1, w/a);
        B = linefromto(B, -1.0, ct2, acn, ct2, w/a);
        B = drawarcn(B, acn, acn, arn - sep/2, w/a, 0, pi/2);
        B = drawarcn(B, acn, acn, arn + sep/2, w/a, 0, pi/2);
    else
        % Traces bitmap - straight line
        B = zeros(nx+2, ny+2); 
        B = linefromto(B, ct1, -1.0, ct1, 2.0, w/a);
        B = linefromto(B, ct2, -1.0, ct2, 2.0, w/a);
    end

    % Make the mesh
    clear mesh
    layer = mklayer(B); % first layer - no vias
    %layer = rmfield(layer, 'conductivity'); % perfect conductor
    layer.pos = 2; % copper ground layer, then dielectric, then mstrip
    mesh.layers(1) = layer;
    layer = mklayer(B, B*vi);
    %layer = rmfield(layer, 'conductivity'); % perfect conductor
    for l = 2:nl
        layer.pos = 1 + l; % copper ground layer, then dielectric, then mstrip
        mesh.layers(l) = layer;
    end

    if doarc
        b1 = findbases(mesh, nx, ny, ct1-sep/2, 0, ct1+sep/2, 0);
        b2 = findbases(mesh, nx, ny, ct2-sep/2, 0, ct2+sep/2, 0);
        b3 = findbases(mesh, nx, ny, 0, ct1-sep/2, 0, ct1+sep/2);
        b4 = findbases(mesh, nx, ny, 0, ct2-sep/2, 0, ct2+sep/2);

        ports = { b1'       b2'       b3'       b4'       };
        portw = { b1'*0-dy  b2'*0-dy  b3'*0-dy  b4'*0-dy  };
    else
        b1 = findbases(mesh, nx, ny, ct1-sep/2, 0, ct1+sep/2, 0);
        b2 = findbases(mesh, nx, ny, ct2-sep/2, 0, ct2+sep/2, 0);
        b3 = findbases(mesh, nx, ny, ct1-sep/2, 1, ct1+sep/2, 1);
        b4 = findbases(mesh, nx, ny, ct2-sep/2, 1, ct2+sep/2, 1);

        ports = { b1'       b2'       b3'       b4'       };
        portw = { b1'*0-dy  b2'*0-dy  b3'*0+dy  b4'*0+dy  };
    end


    [ Y I Z ] = solvey(wg, mesh, ports, portw);

end

% Microstrip params
h = 1e-4   % height above ground
t = 3e-5   % thickness
w = 1e-4   % width
d = 1e-4   % separation
len = w*32 % length
% Arc geometry parameters
ac = w*8  % arc center, both x and y
ar = w*16 % arc radius, center line

% Some useful values
arc_len = ac*2 + ar*pi/2
len_a1 = (ar - d/2 - w/2)*pi/2 % inner arc segment len
len_a2 = (ar + d/2 + w/2)*pi/2 % outer arc segment len
len_c1 = ac*2 + (ar - d/2 - w/2)*pi/2 % first conductor len
len_c2 = ac*2 + (ar + d/2 + w/2)*pi/2 % second conductor len

% angular frequencies
freqs = linspace(1e6,5e10,120)*2*pi;
%% freqs = freqs(1:60);
%% freqs = [ 1e6 1e7 1e8 1e9 1e10 ]*2*pi;

% Simulation results
Yf = [];
Ydf = [];

% angular frequency
for freq = freqs

    step = find(freqs == freq);
    steps = length(freqs);
    fprintf( 'Solving for frequency %.8e, step %i of %i\n', freq, step, steps );

    wavelen = 2*pi/(freq * sqrt(eps0 * mu0))

    fprintf('Running deembedding-l simulation...\n')
    Y1 = simulate(freq, h, t, w, d, w*8, ac, ar, 0);
    fprintf('Running deembedding-2*l simulation...\n')
    Y2 = simulate(freq, h, t, w, d, w*16, ac, ar, 0);

    A1 = y2abcd(Y1);
    A2 = y2abcd(Y2);

    % Obtain the port discontinuity
    [ D dtol ] = splitd(A1*inv(A2)*A1);
    fprintf('Deembedding accuracy is %.8e\n', dtol)

    % Finally the microstrip simulation
    fprintf('Running the main simulation...\n')
    [ Y I Z wg mesh ports ] = simulate(freq, h, t, w, d, len, ac, ar, 1);

    % De-embedded Y
    invD = inv(D);
    Yd = abcd2y(invD*y2abcd(Y)*invD);

    % Add y-parameters for this frequency to the overall results
    Yf = cat(3, Yf, Y);
    Ydf = cat(3, Ydf, Yd);

end

tswrite('arc.y4p', freqs/(2*pi), Yf)
tswrite('arcd.y4p', freqs/(2*pi), Ydf)

%% % Branches:
%% %  1-4  - Y-parameters
%% %  5, 6 - voltage sources
%% %  7, 8 - shunt termination resistors
%% %  9    - series termination resistor
%% branches=[ 2 1 ; 3 1 ; 4 1 ; 5 1 ; 2 1 ; 3 1 ; 4 1 ; 5 1 ; 4 5 ];

%% nb=size(branches, 1);
%% W = [ 0 ; 0 ; 0 ; 0 ; 1 ; -1 ; 0 ; 0 ; 0 ];
%% K = 0*W;
%% Y = zeros(nb,nb);
%% Y(1:4,1:4) = Yd; % de-embedded Y-paramters of the arc
%% Y(5,5) = 1e10;
%% Y(6,6) = 1e10;
%% Y(7,7) = 1/78.5; % from the fieldsolver
%% Y(8,8) = 1/78.5;
%% Y(9,9) = 1/328.5;

%% [ F, V, I, dYd ] = solve(branches,Y,W,K);

%% % excitation voltage vector
%% XV = zeros(size(Z,1), 1);
%% XV(ports{1}) = V(1);
%% XV(ports{2}) = V(2);
%% XV(ports{3}) = -V(3);
%% XV(ports{4}) = -V(4);


%%I = (Z\XV).';
%% I = zeros(1, 30000);
%% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I);
%% trimesh(Tri, X, Y, Z);
%% %trisurf(Tri, X, Y, Z, C);
%% xlim([ -wg.a/wg.nx      wg.a+wg.a/wg.nx ])
%% ylim([ -wg.b/wg.ny      wg.b+wg.b/wg.ny ])
%% zlim([ 0                sum(wg.h) ])

