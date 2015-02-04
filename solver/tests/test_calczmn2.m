function test_calczmn2
% test_calczmn2
%
% Compare results of calczmn2 (sligtly reorganized version of calczmn)
% against the calczmn
% 


% angular frequency
freq=1.0e11;

a=0.01; % x-size of the waveguide
b=0.01; % y-size of the waveguide
h=5e-4; % height of the metallization above ground
c=1e-3; % height of the upper ground
nx=16;  % cells along x
ny=16;  % cells along y

% parametes of the enclosure to pass to calczmn
wg=wgparams(freq,a,b,[h,c-h],nx,ny);

% testing
mi=8;
mj=9;

ni=10;
nj=5;

% x directed testing, x directed source
Zxx=calczmn(wg, mi, mj, 1, ni, nj, 1);
Zxx2=calczmn2(wg, mi, mj, 1, ni, nj, 1);
assertEquals(Zxx, Zxx2, 1e-18);

% y directed testing, x directed source
Zyx=calczmn(wg, mi, mj, 0, ni, nj, 1);
Zyx2=calczmn2(wg, mi, mj, 0, ni, nj, 1);
assertEquals(Zyx, Zyx2, 1e-18);

% x directed testing, y directed source
Zxy=calczmn(wg, mi, mj, 1, ni, nj, 0);
Zxy2=calczmn2(wg, mi, mj, 1, ni, nj, 0);
assertEquals(Zxy, Zxy2, 1e-18);

% y directed testing, y directed source
Zyy=calczmn(wg, mi, mj, 0, ni, nj, 0);
Zyy2=calczmn2(wg, mi, mj, 0, ni, nj, 0);
assertEquals(Zyy, Zyy2, 1e-18);
