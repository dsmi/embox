function test_tlines_iivl
% Test of calc_iivl function - use calc_ivl and do the numerical integration
% of the current

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

% integral multiplier coefficients a*z + b
a = 4.4;
b = -2;

% Linear source multiplier
lsrc = -0.7;

for iobs=1:5
    for jsrc=1:5
	n=100000;
	dl=tl.d(iobs)/n;
	zi=linspace(tl.z(iobs)+dl/2, tl.z(iobs+1)-dl/2, n);
        lm = (zi-tl.z(iobs))*a + b;
	test_ii = sum(calc_ivl(tl, zi, iobs, jsrc, lsrc).*lm)*dl;
	ii      = calc_iivl(tl, iobs, a, b, 1, jsrc, lsrc);
	tol = 1e-100+abs(ii)*1e-9;
	assertEquals(test_ii, ii, tol);
    end
end
