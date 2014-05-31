function test_tlinexy
% Simulate x-directed and y-directed transmission lines of the same
% dimensions. Obviously, the result should be the same.
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


% parametes of the enclosure to pass to mkzmat
wg=wgparams(freq,a,b,h,c,nx,ny);

% Mesh cell size
dx = wg.a/nx;
dy = wg.b/ny;

% x-directed trace
B=zeros(nx+2,ny+2);
B=linefromto(B, -1.0, 0.5, 2.0, 0.5, w/a);
mesh=mkmesh(B);

% Identify ports
b1 = findbases(mesh, nx, ny, 0, 0, 0, 1);
b2 = findbases(mesh, nx, ny, 1, 0, 1, 1);

Y1=solvey(wg, mesh, { b1' b2' }, { b1'*0-dy b2'*0+dy });

% y-directed trace
B=zeros(nx+2,ny+2);
B=linefromto(B, 0.5, -1.0, 0.5, 2.0, w/b);
mesh=mkmesh(B);

% Identify ports
b1 = findbases(mesh, nx, ny, 0, 0, 1, 0);
b2 = findbases(mesh, nx, ny, 0, 1, 1, 1);

Y2=solvey(wg, mesh, { b1' b2' }, { b1'*0-dx b2'*0+dx });

assertEquals(Y1, Y2, 1e-14)
