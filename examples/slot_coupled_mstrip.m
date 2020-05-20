% Microstrip lines coupled through a slot.

addpath(genpath([ pwd, '/..' ]));

a=1e-2; % x-size of the waveguide
b=1e-2; % y-size of the waveguide
h=5e-4; % migrostrip-to-ground
w=b/16; % microstrip width
t=w/6;  % thickness, both line and ground
v=b/4;  % slot width
l=a/16; % slot length

% mesh parameters
nx = 32;  % cells along x
ny = 32;  % cells along y
nl = 1;   % Layers in the microstrip
ng = 1;   % Layers in the ground

% cell size
dx = a/nx; 
dy = b/ny;

% dielectric
er = 4.2;

% mesh bitmaps
b0 = zeros(nx+2, ny+2); % canvas
bl = linefromto(b0, -1.0, 0.5, 2.0, 0.5, w/b);        % line
bs = b0+1;%1 - linefromto(b0, 0.5-l/a, 0.5, 0.5+l/a, 0.5, v/b); % slot

% lower line
line1 = mkhull( bl, traceedges(bl), 1, nl );

% ground plane
plane = mkhull( bs, traceedges(bs), nl+2, ng );

% upper line
line2 = mkhull( bl, traceedges(bl), nl+ng+3, nl );

% Combine the meshes
mesh.layers = [ line1.layers plane.layers line2.layers ];

% angular frequencies
freqs = linspace( 1e6, 2e10, 31 )*2*pi;
%% freqs = ( 2e10 )*2*pi;

% Simulation results
Yf   = [];


for freq = freqs

    fprintf( 'Simulation %i of %i...\n', find(freq == freqs), length(freqs) )
    wavelen = 2*pi/(freq * sqrt(eps0 * er * mu0))

    % dielectric
	 epsd = eps0*debye( er, 0.02, 1e9, freq/(2*pi) ); 

	 % Line layers parameters
	 lh = repmat( t/max( nl, 1 ), 1, nl );
	 le = repmat( eps0, 1, nl );

	 % Ground layers parameters
	 gh = repmat( t/max( ng, 1 ), 1, ng );
	 ge = repmat( epsd, 1, ng );

	 % Simulation stackup
	 %           line1     ground       line2
	 hw = [ h     lh   h    gh     h     lh   h    ];
	 ew = [ eps0  le   epsd ge     epsd  le   eps0 ];

    % parametes of the enclosure to pass to mkzmat
	 wg = wgparams( freq, a, b, hw, nx, ny );
	 wg.weps = ew;
	 wg.cnx  = 4;
	 wg.cny  = 4;
	 wg.Gls0 = 0; % no bottom ground
	 wg.Ggr0 = 0; % no top ground

    % Identify ports
	 b1 = findbases(mesh, nx, ny, 0, 0, 0, 1, @( l ) l <= nl + 1);
	 b2 = findbases(mesh, nx, ny, 1, 0, 1, 1, @( l ) l <= nl + 1);
	 b3 = findbases(mesh, nx, ny, 0, 0, 0, 1, @( l ) l > nl + ng + 2);
	 b4 = findbases(mesh, nx, ny, 1, 0, 1, 1, @( l ) l > nl + ng + 2);

	 ports = { b1' b2' b3' b4' };
	 portw = { b1'*0-dy b2'*0+dy b3'*0-dy b4'*0+dy };

	 [ Y I ] = solvey( wg, mesh, ports, portw );

    % Add y-parameters for this frequency to the overall results
    Yf = cat(3, Yf, Y);

end

tswrite( 'slot_coupled_mstrip.y4p', freqs/(2*pi), Yf )

%% % For the drawing
%% Ip = zeros(50000,1);

%% Ip = I(:,1);

%% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, Ip);
%% %trimesh(Tri, X, Y, Z);
%% trisurf(Tri, X, Y, Z, C);
%% xlim([ -dx  a+dx ])
%% ylim([ -dy  b+dy ])
%% zlim([ 0    a/3    ])
