%
% NEEDS MORE WORK
% "Overclocked capacitor" - impedance of the waveguide formed by a pair of
% round plates. Compare the results against the analytical solution.
%

addpath(genpath([ pwd, '/..' ]));

% Impedance of the round cavity
function Z=rz(freq, eps1, r1, r2, d)

    % Parameters of the plane.
    Yplane = j*freq*eps1/d;
    Zplane = j*freq*mu0*d;
    k = sqrt(-Yplane*Zplane);

    % Analytical solution
    M = [ besselj(0, k*r1) bessely(0, k*r1) ; ...
	    besselj(0, k*r2) bessely(0, k*r2) ];

    dM = [ -besselj(1, k*r1)*k*r1 -bessely(1, k*r1)*k*r1 ; ...
	   -besselj(1, k*r2)*k*r2 -bessely(1, k*r2)*k*r2 ];

    b = [ 1 ; 0 ];
    x = dM\b;

    Zplane = j*freq*mu0*d;
    V = M*x;
    I = -2*pi/Zplane;
    Z = V(1)/I;
end

nx=64;  % cells along x
ny=64;  % cells along y
a=1e-2; % x-size of the waveguide
b=1e-2; % y-size of the waveguide

% Mesh cell size
dx = a/nx;
dy = b/ny;

% Impedances calculated for the freqencies below
Z1f = [];
Z0f = [];
Zcf = [];

% frequencies to solve
freqs = logspace(11, 12, 50);

for freq=freqs

    % angular frequency
    freq

    % parameters of the 'capacitor' - separation, radius
    d  = 1e-6;    % separation
    r2 = nx/3*dx; % outer radius
    r1 = dx/2;
    eps1 = eps0;

    % parametes of the enclosure to pass to mkzmat
    h=[ dx*10 d dx*10 ];
    % h=[ d dx*10 ];
    wg=wgparams(freq,a,b,h,nx,ny);
    wg.cnx = 4; % to make the things a little faster
    wg.cny = 4;
    wg.Gls0 = 0; % no bottom ground
    wg.Ggr0 = 0; % no top ground

    B=zeros(nx+2,ny+2);
    B=drawcir(B, nx/2+1+0.5, ny/2+1+0.5, r2/dx);

    % bottom layer - no via
    clear mesh
    layer = mklayer(B);
    layer = rmfield(layer, 'conductivity'); % perfect conductor
    layer.pos = 1;
    mesh.layers(1) = layer;

    % top layer - with via (which will be a port)
    layer.pos = 2;
    layer.vi=[ nx/2 ];
    layer.vj=[ ny/2 ];
    mesh.layers(2) = layer;

    % one and only layer - with via to the ground (which will be a port)
    %% clear mesh
    %% layer = mklayer(B);
    %% layer = rmfield(layer, 'conductivity'); % perfect conductor
    %% layer.pos = 1;
    %% layer.vi=[ nx/2 ];
    %% layer.vj=[ ny/2 ];
    %% mesh.layers(1) = layer;

    nxy = 2*(size(layer.xi, 1) + size(layer.yi, 1)); % number of x- and y- directed
    ports = { [ nxy + 1 ] }; % the very last b.f. is the via
    portw = { [ dy      ] };

    %% b1 = findbases(mesh, nx, ny, 1, 0, 1, 1)
    %% ports = { [ b1 ] }; % one which connects to the wall
    %% portw = { [ dy ] };

    [ Y I ]=solvey(wg, mesh, ports, portw);
    Z1 = 1/Y

    %% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I);
    %% trisurf(Tri, X, Y, Z, C);
    %% xlim([ 0 a+2*dx ])
    %% ylim([ 0 b+2*dy ])
    %% zlim([ 0 dx*25 ])

    Z0=rz(freq, eps1, r1, r2, d)

    wavelen = 2*pi/(freq * sqrt(eps0 * mu0))

    % capacitor approximation
    C = eps0*r2^2*pi/d;
    Zc = 1/(j*freq*C)

    Z1f = [ Z1f Z1 ];
    Z0f = [ Z0f Z0 ];
    Zcf = [ Zcf Zc ];

end

%semilogx(freqs,imag(Z1f),'-r',freqs,imag(Z0f),'-b',freqs,imag(Zcf),'-g')
%semilogx(freqs,imag(Z0f),'-b',freqs,imag(Zcf),'-g')
