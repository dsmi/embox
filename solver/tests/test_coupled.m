function test_coupled
% x-directed two-conductor (plus ground) transmission line, deembedded by
% simulating line of length L and then 2L

addpath(genpath([ pwd, '/..' ]));

% angular frequency
freq=1.0e9*2*pi;

nx=8;    % cells along x
ny=32;   % cells along y
a=5e-3;  % x-size of the waveguide
b=1e-2;  % y-size of the waveguide
h=[ 5e-4 5e-4 5e-4 ]; % thickness of the dielectric layers - three of them
w=b/ny*4; % stripline width

Y1=simline(freq, nx, ny, a, b, h, w);
Y2=simline(freq, nx*2, ny, a*2, b, h, w);

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

Y=abcd2y(A);

test_Y = ...
  i*[ -1.62219231950047e-1  6.49812343517737e-2   1.63114021528466e-1 -6.53396653406959e-2 ; ...
       6.49812343517729e-2 -1.62219231949220e-1  -6.53396653405713e-2  1.63114021527704e-1 ; ...
       1.63114021528330e-1 -6.53396653406110e-2  -1.62219231950178e-1  6.49812343518555e-2 ; ...
      -6.53396653406106e-2  1.63114021527502e-1   6.49812343517310e-2 -1.62219231949416e-1 ];

assertEquals(test_Y, Y, 1.0e-12)

end

function Y=simline(freq, nx, ny, a, b, h, w)

    % parametes of the enclosure to pass to mkzmat
    wg=wgparams(freq,a,b,h,nx,ny);

    % x-directed trace
    ypos = 0.5+mod(w*ny/b,2)*0.5*b/ny;
    B=zeros(nx+2,ny+2);
    B=linefromto(B, -1.0, ypos, 2.0, ypos, w*ny/(nx*b));
    layer = mklayer(B);
    layer = rmfield(layer, 'conductivity'); % perfect conductor
    mesh = mkmesh(layer, 1, layer, 2); % two identical layers

    % Identify ports : first two ports on one side
    b1 = findbases(mesh, nx, ny, 0, 0, 0, 1, @(l) l == 1);
    b2 = findbases(mesh, nx, ny, 0, 0, 0, 1, @(l) l == 2);
    % and the second two ports on the other
    b3 = findbases(mesh, nx, ny, 1, 0, 1, 1, @(l) l == 1);
    b4 = findbases(mesh, nx, ny, 1, 0, 1, 1, @(l) l == 2);

    % Mesh cell size
    dy = wg.b/ny;

    ports = { b1' b2' b3' b4' };
    portw = { b1'*0-dy b2'*0-dy b3'*0+dy b4'*0+dy };
    Y=solvey(wg, mesh, ports, portw);
end

