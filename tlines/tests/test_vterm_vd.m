function test_vterm_vd
% Test of the calc_vterm_vd function

% Tline parameters
R   = [ 10    5      4     8     3     ]; 
L   = [ 1e-6  1.5e-6 8e-7  5e-7  1e-7  ];
G   = [ 5e2   4e2    5e1   1e2   2e2   ];
C   = [ 1e-11 3e-11  5e-11 5e-12 4e-11 ];
len = [ 1e-3  4e-3   4e-5  2e-3  6e-3  ];

% Frequency and angular frequency
freq=1e6;
afreq=2*pi*freq;

% Characteristic impedances
Z0=sqrt( (R+j*afreq*L)./(G+j*afreq*C) );

% Propagation constants
tl_k=sqrt( (R+j*afreq*L).*(G+j*afreq*C) );

% Output the tline parameters
if 0,
	for l=1:length(R),
		printf("Tline %i\n",l);
		printf("    R=%e, L=%e, G=%e, c=%e\n", R(l), L(l), G(l), C(l));
		printf("    len=%e\n", len(l));
		printf("    k=%e + %ei\n", real(tl_k(l)), imag(tl_k(l)));
		printf("    Z0=%e + %ei\n", real(Z0(l)), imag(Z0(l)));
	end
	printf("\n");
end

% Output the matched terminations
if 0,
	% Matched terminations
	Rterm1=real(Z0(1));
	Lterm1=imag(Z0(1))/afreq;
	Rterm5=real(Z0(5));
	Lterm5=imag(Z0(5))/afreq;

	printf("Matched termination for line 1: R=%e, L=%e\n", Rterm1, Lterm1);
	printf("Matched termination for line 5: R=%e, L=%e\n", Rterm5, Lterm5);
	printf("\n");
end

% Setup the endpoint coordinates to be used by calculator.
tl_z = [ 0 cumsum(len) ];

% Matched termination at both ends
Gls1=0;
GgrN=0;

% Run the coefficients precomputation
tl = calc_tlines(tl_z, Z0, tl_k, Gls1, GgrN);

vt3 = calc_vterm_vd(tl,3);

% This is how the test value has been obtained - this is quite slow,
% so the result is just recorded below
%% n=1000000;
%% test_vt3 = 0;
%% dl=tl.d(3)/n;
%% for i=1:n
%%     test_vt3 = test_vt3 + calc_vterm_v(tl,tl.z(3)+i*dl,3)*dl;
%% end

test_vt3 = -2.09323310480321e-005 - i*2.56585068783313e-006;

assertEquals(test_vt3,vt3,1e-13);

