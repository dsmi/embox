%
% Stripline with two reference planes.
%


addpath(genpath([ pwd, '/..' ]));

% The geometry
inch2meter = 2.54e-2;
mil2meter = 1.0e-3*inch2meter;
t  = 0.675*mil2meter % stripline thickness
lw = 8.0*mil2meter;   % stripline width
d1 = 5*mil2meter  % line to lower plane
d2 = 5*mil2meter  % line to upper plane
d = (d1 + d2 + t) % plane-to-plane separation
w = 100.0*mil2meter;  % cavity width
l = 200.0*mil2meter;  % cavity length

% Mesh settings
cny = 50;
cnx = cny*2;
xo = 30; % edge-to-enclosure doubled
yo = 30;

% Mesh cell size
dx = l/cnx;
dy = w/cny;

% Dielectric params
lt = 0.02;
er0 = 4.3;
fr = 1e9;

% metal conductivity
sigma = 5.8e7;

% angular frequencies
%freqs = 1e9*2*pi;
freqs = linspace(1e7, 4e10, 100)*2*pi;

% Simulation results
Yf = [];

for ifr = 1:length(freqs)

    freq = freqs(ifr);
    wavelen = 2*pi/(freq * sqrt(eps0 * mu0));
    fprintf( 'Solving for frequency %.8e, step %i of %i, wlen %.8e\n', freq, ifr, length(freqs), wavelen );

    % Layers stack
    h    = [ d d1 t d2 d ];
    er   = debye(er0, lt, fr, freq/(2*pi));
    weps = [ eps0*er eps0*er eps0*er eps0*er eps0*er ];

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

    % Port vias
    B2a = drawline(B0, lxb, yc, lxb + 1.0e-5, yc, lw/dy);
    B2 = drawline(B2a, lxe, yc, lxe + 1.0e-5, yc, lw/dy);

    BL = drawline(B0, lxb, yc, lxe, yc, lw/dy);
    BV = traceedges(BL); % trace vias

    % Make mesh from the layers
    mesh.layers = struct( [] );
    mesh.layers(end + 1) = mklayer(B1, 0*B2, 1, ccopper); % Lower wall
    mesh.layers(end + 1) = mklayer(BL, B2, 2, ccopper); % Lower wall of the line and vias
    mesh.layers(end + 1) = mklayer(BL, BV, 3, ccopper); % Upper wall of the line
    mesh.layers(end + 1) = mklayer(B1, B2, 4, ccopper); % Upper wall and vias

    % Identify ports
    b1 = findbases(mesh, nx, ny, 0.0, 0, 0.5, 1, @(l) l == 2);
    b2 = findbases(mesh, nx, ny, 0.5, 0, 1.0, 1, @(l) l == 2);
    b3 = findbases(mesh, nx, ny, 0.0, 0, 0.5, 1, @(l) l == 4);
    b4 = findbases(mesh, nx, ny, 0.5, 0, 1.0, 1, @(l) l == 4);

    ports = { b1' b2' b3' b4' };
    portw = { b1'*0+dy b2'*0+dy b3'*0+dy b4'*0+dy }; % via currents are scaled by 1/dx
    [ Y I ]=solvey(wg, mesh, ports, portw);

    Yf = cat(3, Yf, Y);

end

tswrite([ 'stripline_in_cavity_3d', '.y4p' ], freqs/(2*pi), Yf, 'Y', 50);

%% I = zeros(100000, 1);
%% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I(:,1));
%% %trisurf(Tri, X, Y, Z, C);
%% trimesh(Tri, X, Y, Z);
%% xlim([ -(wg.a/wg.nx)  wg.a+2*(wg.a/wg.nx) ])
%% ylim([ -(wg.b/wg.ny)  wg.b+2*(wg.b/wg.ny) ])
%% zlim([ 0 3*d ])

