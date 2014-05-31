function test_tline
% This is an accelerated version (coarser mesh) of examples/tline. Simulate
% (and deembed) a stripline and compare with values from the transmission
% line model.
% 
% x-directed transmission line, deembedded by simulating line of length L
% and then 2L
%

% angular frequency
freq=5.0e10;

nx=16;  % cells along x
ny=16;  % cells along y
a=1e-2; % x-size of the waveguide
b=1e-2; % y-size of the waveguide
h=5e-4; % height of the metallization above ground
c=1e-3; % height of the upper ground
w=b/ny*4; % stripline width

Y1=simline(freq, nx, ny, a, b, h, c, w);
Y2=simline(freq, nx*2, ny, a*2, b, h, c, w);

A1=y2a(Y1);
A2=y2a(Y2);

Add = A1*inv(A2)*A1; % double port discontinuity
Ad  = [ Add(1,1) Add(1,2) ; Add(2,1)/2 Add(2,2) ]; % port discontinuity
A = inv(Ad)*A1*inv(Ad); % de-embedded line of length l

% Characteristic impedance and delay from simulation
Z0s=sqrt(A(1,2)./A(2,1));
tds=acos(A(1,1))./freq;

wavelen = 2*pi/(freq * sqrt(eps0 * mu0));

l=a;   % length
b=c;   % distance between planes

%Z0=30*pi*b/w/(1+0.441*b/w)
% Calculated using mom: see mom/examples/cstrip
C = 1.0425e-010;
L = 1.0673e-007;
%% C = 5.9988e-011;
%% L = 1.8548e-007;
Z0 = sqrt(L/C);

beta=freq*sqrt(eps0*mu0); % wavenumber
gamma=j*beta; % propagation constant
vp=freq/beta;
td=l/vp;

Z11=Z0/tanh(gamma.*l);
Z12=Z0/sinh(gamma.*l);
Ztl=[ Z11 Z12 ; Z12 Z11 ];
Ytl=inv(Ztl);

% Values from the tline model
A0=y2a(Ytl);
assertEquals(Z0, Z0s, 2);
assertEquals(td, tds, 5e-15);

end

function Y=simline(freq, nx, ny, a, b, h, c, w)

    % parametes of the enclosure to pass to mkzmat
    wg=wgparams(freq,a,b,h,c,nx,ny);

    % x-directed trace
    B=zeros(nx+2,ny+2);
    B=linefromto(B, -1.0, 0.5, 2.0, 0.5, w/a);
    mesh=mkmesh(B);

    % Identify ports
    b1 = findbases(mesh, nx, ny, 0, 0, 0, 1);
    b2 = findbases(mesh, nx, ny, 1, 0, 1, 1);

    % Mesh cell size
    dy = wg.b/ny;

    Y=solvey(wg, mesh, { b1' b2' }, { b1'*0-dy b2'*0+dy });
end
