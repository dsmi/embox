%
% x-directed differential microstrip with optional tabs of nonzero thickness
% represented by two or more metal layers, deembedded by simulating line of
% length L and then 2L
%

addpath(genpath([ pwd, '/..' ]));

mil2meter = 2.54e-5;

% Dimensions
lnpar.w  = 4*mil2meter;     % trace width
lnpar.l  = 68*mil2meter;    % trace length
lnpar.g  = 0.5*mil2meter;   % topmost layer thickness
lnpar.t  = 1.9*mil2meter;   % trace thickness
lnpar.d  = 2.7*mil2meter;   % height above ground
lnpar.s  = 12*mil2meter;    % center-to-center separation
lnpar.ts = lnpar.l/8;       % tab-to-tab separation (on different sides)
lnpar.tw = 4*mil2meter;     % tab width
lnpar.tl = 6*mil2meter;     % tab length (from the center)


% Mesh options
lnpar.nl = 4;        % number of the metal layers in the stripline mesh
lnpar.nx = (68/4)*2*2;   % cells along x
lnpar.ny = (68/4)*2*2;   % cells along y
lnpar.a  = lnpar.l;  % x-size of the waveguide
lnpar.b  = lnpar.a;  % y-size of the waveguide


% frequency sweep (frequency is angular)
nf = 41;
minf = 1e6;
maxf = 4e10;
freqs = [ minf linspace(maxf/(nf-1), maxf, nf-1) ]*2*pi;

resultsName = 'tabbed.s4p';

function [ Y I wg mesh ports portw ] = simline(freq, lnpar)

    % number of 'steps' between metal layers
    nls = lnpar.nl - 1; 

    %  Layers stack
    %     ground  separation metal                       topmost  air
    h = [ lnpar.d lnpar.d    repmat(lnpar.t/nls, 1, nls) lnpar.g  lnpar.g ];

    % Ground layer is copper (not an ideal conductor)
    epsg = eps0 - j*ccopper/freq; % permittivity for copper
    eps1 = eps0 * debye2(4.0, 0.02, 1e9*2*pi, freq); % between line and ground
    eps2 = eps0 * debye2(3.8, 0.02, 1e9*2*pi, freq); % top dielectric
    weps = [ epsg eps1 repmat(eps2, 1, nls) eps2 eps0 ];

    % parametes of the enclosure to pass to mkzmat
    nx      = lnpar.nx;
    ny      = lnpar.ny;
    wg      = wgparams(freq, lnpar.a, lnpar.b, h, nx, ny);
    wg.weps = weps; 
    wg.Ggr0 = 0; % no top ground
    wg.Gls0 = 0; % no bottom ground (we have layer of copper)
    wg.cnx  = 4; % to speed up things (change to 8 for accuracy)
    wg.cny  = 4;

    % Mesh cell size
    dx = wg.a/lnpar.nx;
    dy = wg.b/lnpar.ny;

    w  = lnpar.w;
    nl = lnpar.nl;
    l  = lnpar.l;
    s2 = lnpar.s/2; % half-separation

    % x-directed traces
    B0 = zeros(nx+2, ny+2);
    BL1 = drawline(B0, 0, (ny+2)/2-s2/dy, nx+2, (ny+2)/2-s2/dy, w/dy);
    BL2 = drawline(B0, 0, (ny+2)/2+s2/dy, nx+2, (ny+2)/2+s2/dy, w/dy);

    % tabs
    ts = lnpar.ts;
    tw = lnpar.tw;
    tl = lnpar.tl;
    nt = round( l / ts );
    for ti=1:nt
      cy1 = (ny+2)/2-s2/dy;                 % y of the first line center
      cy2 = (ny+2)/2+s2/dy;                 % y of the second line center
      tcx = (nx+2)/2-(l/2-(ti-0.5)*ts)/dx;  % tab center x  
      BL1 = drawline(BL1, tcx, cy1, tcx, cy1+(-1)^ti*tl/dy, tw/dx);
      BL2 = drawline(BL2, tcx, cy2, tcx, cy2+(-1)^ti*tl/dy, tw/dx);
    end
    
    % Make mesh from the layers
    mesh.layers = struct( [] );

    % trace mesh
    for li=1:lnpar.nl
        
        BL = BL1 + BL2;
        
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
    b1 = findbases(mesh, nx, ny, 0, 0.0, 0, 0.5);
    b2 = findbases(mesh, nx, ny, 0, 0.5, 0, 1.0);
    b3 = findbases(mesh, nx, ny, 1, 0.0, 1, 0.5);
    b4 = findbases(mesh, nx, ny, 1, 0.5, 1, 1.0);

    ports = { b1' b2' b3' b4' };
    portw = { b1'*0-dy b2'*0-dy b3'*0+dy b4'*0+dy };
    [ Y I ]=solvey(wg, mesh, ports, portw);
    
end

% Obtain port discontinuity matrix via simulation of the calibration
% standards of single and double lengths
function [ Ad R Y1 Y2 Yd wg mesh ports portw ] = calcporta(freq, lnpar)

    % signle-l
    lnpar1 = lnpar;
    Y1 = simline(freq, lnpar1);

    % double-l
    lnpar2 = lnpar;
    lnpar2.l  = lnpar.l*2;
    lnpar2.a  = lnpar.a*2;
    lnpar2.nx = lnpar.nx*2;
    [ Y2 I2 wg mesh ports portw ] = simline(freq, lnpar2);

    A1=y2abcd(Y1);
    A2=y2abcd(Y2);

    Add = A1*inv(A2)*A1; % double port discontinuity
    N = size(Add,1)/2;   % number of ports on each side
    A = Add(1:N,1:N);
    B = Add(1:N,N+1:end);
    C = Add(N+1:end,1:N);
    D = Add(N+1:end,N+1:end);

    Yd = C*0.5; % port discontinuity admittance
    Ad  = [ A B ; Yd D ]; 
    
    % Estimate accuracy of deembedding.
    % Expected values of the shunt ABCD
    % are: A = I; B = 0; C = 1/Zs; D = I
    R = norm( A - eye(N, N) ) + norm( B ) + norm( D - eye(N, N) );

end

function writesp(fileName, freq, Y)
  znorm1 = ones(1,size(Y,1));
  tswrite(fileName, freq/(2*pi), renorms(y2s(Y), znorm1, znorm1*50), 'S', 50)
end

% Simulation results
Yf = [];

for ifr = 1:length(freqs)

    freq = freqs(ifr);
    wavelen = 2*pi/(freq * sqrt(eps0 * mu0));
    fprintf( 'Solving for frequency %.4e wavelen %.4e, step %i of %i\n', ...
             freq, wavelen, ifr, length(freqs) );

    % deembedding simulations
    [ Ad R Y1 Y2 Yd wg mesh ports portw ] = calcporta(freq, lnpar);

    fprintf( 'De-embdeeing accuracy %.8e\n', R );

    And = y2abcd(Y2);

    invAd = inv(Ad);
    A = invAd * And * invAd;

    Y = abcd2y(A);

    % Add y-parameters for this frequency to the overall results
    Yf = cat(3, Yf, Y);
    
    % Characteristic impedance and delay from simulation
    N=size(A,1)/2; % number of ports on each side
    A12=A(1:N,N+1:end);
    A21=A(N+1:end,1:N);
    Z0s=sqrt(A12*inv(A21))
    tds=acos(A(1,1))./freq

    writesp( resultsName, freqs, Yf);
    
end

%% [ Y I wg mesh ports portw ] = simline(1e10, lnpar);
%% Y

%% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I(:,1), @abs);
%% hplot = trisurf(Tri, X, Y, Z, C);

%% xlim([ -(wg.a/wg.nx)  wg.a+2*(wg.a/wg.nx) ])
%% ylim([ -(wg.b/wg.ny)  wg.b+2*(wg.b/wg.ny) ])
%% zlim([ -(wg.b/wg.ny)  wg.b+2*(wg.b/wg.ny) ])

%% colormap(jet)
%% % caxis([ -7e5 7e5 ]) 
%% set(hplot,'edgecolor','none')
%% shading interp
%% colorbar
