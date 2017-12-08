function test_tlines_ivid
% Test of calc_ivid function - use calc_vid and do the numerical integration
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

iobs=2;
jsrc=3;

for iobs=1:5
    for jsrc=1:5
        n=10001;
        dl=tl.d(iobs)/n;
        zi=linspace(tl.z(iobs) + dl/2, tl.z(iobs+1) - dl/2, n);

        test_iv = sum(calc_vid(tl, zi, iobs, jsrc))*dl;
        iv      = calc_ivid(tl, iobs, jsrc);
        tol = 1e-15+abs(iv)*1e-10;
        assertEquals(test_iv, iv, tol);
    end
end


