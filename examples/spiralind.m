%
% Draw and simulate square spiral inductor.
%

addpath(genpath([ pwd, '/..' ]));

% Inductor parameters
w  = 1.0e-5;  % trace width
s  = 1.0e-5;  % spacing
th = 2.0e-6;  % inductor thickness
h  = 2.0e-5;  % dielectric thickness --  between the inductor and the substrate
t2 = 2.0e-6;  % thickness of the underneath trace
h2 = 1.0e-5; % spacing between the underneath trace and the substrate
nseg = 15;    % number of segments

a = w*32; % x-size of the enclosure/waveguide
b = w*32; % y-size of the enclosure/waveguide

% Mesh options
nx = 64;    % cells along x
ny = 64;    % cells along y
nl = 1;     % number of the layers in the traces
nv = 3;     % number of the layers in the feeding vias

% Start meshing -- this is an empty bitmap, to be used for drawing
b0 = zeros(nx+2, ny+2);

% Segment lenght, increased while we are drawing
l = w*5;

% Starting point
x0 = x = a/2;
y0 = y = b/2 - l/2;

% Segment direction
sdx = 1;
sdy = 0;

% Per-step length increase
dl = (w/2+s/2);

% Inductor bitmap, empty so far
bs = b0;

% Draw the inductor segment by segment
for seg=1:nseg

    % Length of the segment to draw, first and last are adjusted
    segl = l - l/2*( seg == 1 ) - (l/2 + dl)*( seg == nseg );
    
    x1 = x + sdx*segl;
    y1 = y + sdy*segl;

    xw2 = sdx*w/2; % increase each line length by w/2
    yw2 = sdy*w/2; % to avoid the cut corners
    bs = linefromto( bs, (x-xw2)/a, (y-yw2)/b, (x1+xw2)/a, (y1+yw2)/b, w/a );
    
    [ x, y ] = deal( x1, y1 );

    % Change direction for the next segment, increase length
    [ sdx, sdy ] = deal( -sdy, sdx );
    l = l + dl;
    
end

% Line from the inductor to the wall
bs = linefromto( bs, x1/a, y1/a, x1/a, 1.0 + 1.0/ny, w/a );

% Position of the via connecting the  underneath trace outside the inductor
x3 = x0;
y3 = y0 - ceil(nseg/4)*(w+s);

% Line from the via to the other wall
bs = linefromto( bs, x3/a, (y3+w/2)/a, x3/a, -1.0/ny, w/a );

% Draw the feeding vias on the vias bitmap
bv = linefromto( b0, (x3-w/2)/a, y3/b, (x3+w/2)/a, y3/b, w/a );
bv = linefromto( bv, (x0-w/2)/a, y0/b, (x0+w/2)/a, y0/b, w/a );

% The final bitmap -- the underneath trace
bu = linefromto( b0, x3/a, (y3-w/2)/b, x0/a, (y0+w/2)/b, w/a );

undermesh = mkhull( bu, traceedges( bu ), 2, nl );
viamesh = mkhull( 0*bv, bv, nl+2, nv );
spiralmesh = mkhull( bs, traceedges( bs ), nl + nv + 2, nl );

mesh.layers = [ spiralmesh.layers viamesh.layers undermesh.layers ];

mesh = mergelayers(mesh);

% angular frequencies
freqs = linspace( 1e7, 1e10, 21 );

% Simulation results
Yf   = [];

for freq = freqs

    fprintf( 'Simulation %i of %i...\n', find(freq == freqs), length(freqs) )

    wavelen = 2*pi/(freq * sqrt(eps0 * 4.0 * mu0));
    fprintf( 'freq = %g GHz, wlen = %g\n', freq/(2*pi)/1e9, wavelen )

    % 'Permittivity' of the conducting substrate. Now it is copper
    epss = eps0 - j*ccopper/freq;

    % Dielectric between the substrate and the inductor. er = 4, lt = 0.02
    epsd = eps0*debye( 4.0, 0.02, 1e9, freq/(2*pi) );

    nlo = ones( 1, nl );
    nvo = ones( 1, nv );

    % Layers: substrate  diel.  lower trace  vias               inductor
    wh   = [  h         h2     nlo*t2/nl    nvo*(h-h2-t2)/nv   nlo*th/nl ];
    weps = [  epss      epsd   nlo*epsd     nvo*epsd           nlo*eps0  ];

    % enclosure/waveguide parameters
    wg      = wgparams(freq, a, b, wh, nx, ny);
    wg.weps = weps; 
    wg.cnx  = 4;
    wg.cny  = 4;
    wg.Ggr0 = 0; % no top ground

    b1 = findbases( mesh, nx, ny, 0, 0, 1, 0 );
    b2 = findbases( mesh, nx, ny, 0, 1, 1, 1 );

    ports = { b1' b2' };
    portw = { b1'*0-a/nx b2'*0+a/nx };

    [ Y I ] = solvey( wg, mesh, ports, portw );

    Yf = cat( 3, Yf, Y );

end

tswrite( 'inductor.y2p', freqs/(2*pi), Yf )

%% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I(:,1));
%% trisurf(Tri, X, Y, Z, C);
%% xlim([ -a/nx-1e-10  a+a/nx+1e-10 ])
%% ylim([ -b/ny-1e-10  b+b/ny+1e-10 ])
%% zlim([ 0    a+2*a/nx    ])
%% colormap(jet)
