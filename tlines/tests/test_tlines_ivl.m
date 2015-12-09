function test_tlines_ivl
% Test of calc_ivl function

% Tline parameters
R   = [ 10    5      4     8     3     ]; 
L   = [ 1e-6  1.5e-6 8e-7  5e-7  1e-7  ];
G   = [ 5e2   4e2    5e1   1e2   2e2   ];
C   = [ 1e-11 3e-11  5e-11 5e-12 4e-11 ];
len = [ 1e-3  4e-3   4e-4  2e-3  6e-3  ];

% Frequency and angular frequency
freq=1e6;
afreq=2*pi*freq;

% Characteristic impedances
Z0=sqrt( (R+j*afreq*L)./(G+j*afreq*C) );

% Propagation constants
tl_k=sqrt( (R+j*afreq*L).*(G+j*afreq*C) );

% Setup the endpoint coordinates to be used by calculator.
tl_z = [ 0 cumsum(len) ];

% Matched termination at both ends
Gls1=0;
GgrN=0;

% Run the coefficients precomputation
tl = calc_tlines(tl_z, Z0, tl_k, Gls1, GgrN);

% Linear source coefficient
l = 2/tl.d(3);

n=100000;
dl=tl.d(3)/n;
zi=linspace(tl.z(3)+dl/2, tl.z(4)-dl/2, n);
lm=(zi-tl.z(3))*l;

test_v3 = sum(calc_iv(tl, tl_z(3)+3e-4, 3, zi, 3).*lm)*dl;
v3 = calc_ivl(tl, tl_z(3)+3e-4, 3, 3, l);
assertEquals(test_v3,v3,1e-12);

test_v2 = sum(calc_iv(tl, tl_z(2)+1e-3, 2, zi, 3).*lm)*dl;
v2 = calc_ivl(tl, tl_z(2)+1e-3, 2, 3, l);
assertEquals(test_v2,v2,1e-14);

test_v1 = sum(calc_iv(tl, tl_z(1)+3e-4, 1, zi, 3).*lm)*dl;
v1 = calc_ivl(tl, tl_z(1)+3e-4, 1, 3, l);
assertEquals(test_v1, v1, 1e-14);

test_v4 = sum(calc_iv(tl, tl_z(4)+5e-4, 4, zi, 3).*lm)*dl;
v4 = calc_ivl(tl, tl_z(4)+5e-4, 4, 3, l);
assertEquals(test_v4,v4,1e-14);

test_v5 = sum(calc_iv(tl, tl_z(5)+2e-3, 5, zi, 3).*lm)*dl;
v5 = calc_ivl(tl, tl_z(5)+2e-3, 5, 3, l);
assertEquals(test_v5,v5,1e-14);

% try to place the source into the first tline
n=1000000;
dl=tl.d(1)/n;
zi=linspace(tl.z(1)+dl/2, tl.z(2)-dl/2, n);
lm=(zi-tl.z(1))*l;

test_v5v1 = sum(calc_iv(tl, tl_z(5)+2e-3, 5, zi, 1).*lm)*dl;
v5v1 = calc_ivl(tl, tl_z(5)+2e-3, 5, 1, l);
assertEquals(test_v5v1,v5v1,1e-12);
