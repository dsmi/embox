% x-directed transmission line of nonzero thickness represented by
% by two or more metal layers, deembedded by simulating line of length L
% and then 2L

addpath(genpath([ pwd, '/..' ]));

function [ Y I wg mesh ]=simline(freq, nx, ny, a, b, h, weps, w, nl)

    % parametes of the enclosure to pass to mkzmat
    wg      = wgparams(freq,a,b,h,nx,ny);
    wg.weps = weps; 
    wg.Gls0 = 0; % no bottom ground
    wg.Ggr0 = 0; % no top ground

    % Mesh cell size
    dy = wg.b/ny;

    % Make mesh from the layers
    clear mesh
    for lidx=1:nl
        % x-directed trace
        B0 = zeros(nx+2,ny+2);
        B = drawline(B0, 0, (ny+2)/2, nx+2, (ny+2)/2, w/dy);
        if lidx > 1 && lidx < nl
            B = B - drawline(B0, 0, (ny+2)/2, nx+2, (ny+2)/2, w/dy - 2);
        end
        layer = mklayer(B);
        layer.pos = lidx + 1;
        mesh.layers(lidx) = layer;
    end
    
    % Identify ports
    b1 = findbases(mesh, nx, ny, 0, 0, 0, 1); 
    b2 = findbases(mesh, nx, ny, 1, 0, 1, 1);

    ports = { b1' b2'};
    portw = { b1'*0-dy b2'*0+dy };
    [ Y I ]=solvey(wg, mesh, ports, portw);
end

w=5.0e-5;   % stripline width
l=2.0e-3;   % stripline length
t=1.5e-5;   % stripline thickness
d=3e-5;     % height above ground
nl=2;       % number of the metal layers in the stripline model
nx=16;      % cells along x
ny=64;      % cells along y
a=l;        % x-size of the waveguide
b=w*8;      % y-size of the waveguide

% angular frequencies
%freqs = logspace(6,11,300)*2*pi;
freqs = linspace(1e6,1e11,300)*2*pi;
freqs = [ 1e9 1e10 1e11 ]*2*pi;

% Simulation results
Y1f = [];
Y2f = [];
Yf = [];

for freq = freqs

    step = find(freqs == freq);
    steps = length(freqs);
    fprintf( 'Solving for frequency %.8e, step %i of %i\n', freq, step, steps );

    % Layer of copper at the bottom is the ground - it is not an ideal conductor
    epsc = eps0 - j*ccopper/freq; % permittivity for copper - incorporates conductivity
    h    = [ d     d     repmat(t/(nl-1), 1, nl-1)  d    ]; 
    weps = [ epsc  eps0  repmat(eps0,     1, nl-1)  eps0 ];

    % simulate l and 2*l to deembed
    [ Y1 I1 wg mesh ] = simline(freq, nx,   ny, a,   b, h, weps, w, nl);
    [ Y2 I2 wg mesh ] = simline(freq, nx*2, ny, a*2, b, h, weps, w, nl);

    A1=y2abcd(Y1);
    A2=y2abcd(Y2);

    Add = A1*inv(A2)*A1; % double port discontinuity
    N = size(Add,1)/2;  % number of ports on each side
    Add11=Add(1:N,1:N);
    Add12=Add(1:N,N+1:end);
    Add21=Add(N+1:end,1:N);
    Add22=Add(N+1:end,N+1:end);
    Ad  = [ Add11 Add12 ; Add21*0.5 Add22 ]; % port discontinuity
    A = inv(Ad)*A1*inv(Ad); % de-embedded line of length l
    Y = abcd2y(A) % admittance from simulation

    % No deembedding for a while
    %% [ Y I wg mesh ] = simline(freq, nx, ny, a, b, h, weps, w, nl);
    %% A = y2abcd(Y);

    % Add y-parameters for this frequency to the overall results
    Yf = cat(3, Yf, Y);

    % Capture y-parameters before de-embedding as well
    %% Y1f = cat(3, Y1f, Y1);
    %% Y2f = cat(3, Y2f, Y2);

    % Characteristic impedance and delay from simulation
    N=size(A,1)/2; % number of ports on each side
    A12=A(1:N,N+1:end);
    A21=A(N+1:end,1:N);
    Z0s=sqrt(A12*inv(A21))
    tds=acos(A(1,1))./freq

    wavelen = 2*pi/(freq * sqrt(eps0 * mu0))
    skin_depth = sqrt(2./(freq*mu0*ccopper))

    % Calculated using the 2d fieldsolver
    C = 3.8918761e-011;
    L = 2.8589042e-007;
    Ro = 2.2986668e+001;
    Rs = 3.9169897e-003;

    Zo = j*freq*L + Ro + Rs*(1+j)*sqrt(freq/(2*pi)); % impedance per len
    Yo = j*freq*C; % admittance per len
    Z0 = sqrt(Zo/Yo)

    gamma = sqrt(Zo*Yo); % propagation constant
    vp = j*freq/gamma;
    td = l/vp;

    Z11 = Z0/tanh(gamma.*l);
    Z12 = Z0/sinh(gamma.*l);
    Ztl = [ Z11 Z12 ; Z12 Z11 ];
    Ytl = inv(Ztl);
end

%% freqs_hz = freqs/(2*pi);
%% writey2p('microstrip2mm_embox.y2p', freqs_hz, Yf);
%% writey2p('microstrip2mm_embox_1.y2p', freqs_hz, Y1f);
%% writey2p('microstrip2mm_embox_2.y2p', freqs_hz, Y2f);

%% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I);
%% trisurf(Tri, X, Y, Z, C);
%% xlim([ 0  wg.a+3*(wg.a/nx) ])
%% ylim([ 0  wg.b+3*(wg.b/ny) ])
%% zlim([ d  2*d+t*2          ])

% The improved colormap
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
%colorbar
