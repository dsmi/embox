function test_thick
% x-directed transmission line of nonzero thickness represented by
% by a pair of layers, deembedded by simulating line of length L
% and then 2L

%addpath(genpath([ pwd, '/..' ]));

% angular frequency
freq=5.0e10*2*pi;

nx=8;     % cells along x
ny=32;    % cells along y
a=1.25e-3;  % x-size of the waveguide
b=1e-2;     % y-size of the waveguide
h=[ 5e-4 1e-5 5e-4 ]; % thickness of the dielectric layers - three of them
w=b/ny*4; % stripline width

[ Y1 I1 ]=simline(freq, nx, ny, a, b, h, w);
[ Y2 I2 ]=simline(freq, nx*2, ny, a*2, b, h, w);

A1=y2abcd(Y1);
A2=y2abcd(Y2);

Add = A1*inv(A2)*A1; % double port discontinuity
N = size(Add,1)/2; % number of ports on each side
Add11=Add(1:N,1:N);
Add12=Add(1:N,N+1:end);
Add21=Add(N+1:end,1:N);
Add22=Add(N+1:end,N+1:end);
Ad  = [ Add11 Add12 ; Add21*0.5 Add22 ]; % port discontinuity
A = inv(Ad)*A1*inv(Ad); % de-embedded line of length l
Y = abcd2y(A); % admittance from simulation

% Characteristic impedance and delay from simulation
A12=A(1:N,N+1:end);
A21=A(N+1:end,1:N);
Z0s=sqrt(A12*inv(A21));
tds=acos(A(1,1))./freq;

wavelen = 2*pi/(freq * sqrt(eps0 * mu0));

l=a;   % length

C = 6.2415337e-011;
L = 1.7826551e-007;
Ro = 2.7583998e-001;
Rs = 2.1258788e-004;
%% conductivity = 5.8e7; % Copper
%% skin_depth = sqrt(2./(freq*mu0*conductivity));
%% Rs = 1./(w*skin_depth*conductivity);
Zo = j*freq*L + Ro + Rs*(1+j)*sqrt(freq/(2*pi)); % impedance per len
Yo = j*freq*C; % admittance per len
Z0 = sqrt(Zo/Yo);

gamma=sqrt(Zo*Yo); % propagation constant
vp=j*freq/gamma;
td=l/vp;

Z11=Z0/tanh(gamma.*l);
Z12=Z0/sinh(gamma.*l);
Ztl=[ Z11 Z12 ; Z12 Z11 ];
Ytl=inv(Ztl);

% slightly differs from the tline model result, but we use quite coarse mesh
test_Y = ...
   [ (5.79897500294733e-006 -i*4.66359233478252e-003) (-3.98480013709619e-006+i*1.81062642911279e-002) ;...
     (-3.98480013709619e-006+i*1.81062642911279e-002) (5.79897500294733e-006-i*4.66359233478252e-003) ];

assertEquals(test_Y, Y, 1.0e-14)

end

function [ Y I ]=simline(freq, nx, ny, a, b, h, w)

    % parametes of the enclosure to pass to mkzmat
    wg=wgparams(freq,a,b,h,nx,ny);

    % x-directed trace
    ypos = 0.5+mod(w*ny/b,2)*0.5*b/ny;
    B=zeros(nx+2,ny+2);
    B=linefromto(B, -1.0, ypos, 2.0, ypos, w*ny/(nx*b));
    layer = mklayer(B, B*0, 1, ccopper);
    mesh = mkmesh(layer, 1, layer, 2); % two identical layers
    %mesh = mkmesh(layer, 1); % one layer

    % Identify ports
    b1 = findbases(mesh, nx, ny, 0, 0, 0, 1); 
    b2 = findbases(mesh, nx, ny, 1, 0, 1, 1);

    % Mesh cell size
    dy = wg.b/ny;

    ports = { b1' b2'};
    portw = { b1'*0-dy b2'*0+dy };
    [ Y I ]=solvey(wg, mesh, ports, portw);
end
