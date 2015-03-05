function test_tlines_ivd
% Test of calc_ivd function - use calc_iv and simulate the distributed
% source by numerical integration

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

for iobs=1:5
    for jsrc=1:5
	n=100000;
	dl=tl.d(jsrc)/n;
	zi=linspace(tl.z(jsrc), tl.z(jsrc+1), n);
	test_i = sum(calc_iv(tl, tl_z(iobs)+len(iobs)*0.3, iobs, zi, jsrc))*dl;
	i = calc_ivd(tl, tl_z(iobs)+len(iobs)*0.3, iobs, jsrc);
	tol = 1e-12+abs(i)*1e-6;
	assertEquals(test_i,i,tol);
    end
end
