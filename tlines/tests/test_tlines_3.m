%function test_tlines_3
% One more test based on the exponent overflow problem
  
z  = [ 0.00000000000000e+000  5.00000000000000e-004  5.01000000000000e-004 1.00100000000000e-003 ];
Z0 = [ 0.00000000000000e+000 - 2.55207916846424e+009i  0.00000000000000e+000 - 2.55207916846424e+009i 0.00000000000000e+000 - 2.55207916846424e+009i ];
k  = [ 7.09892756048521e+005  7.09892756048521e+005  7.09892756048521e+005 ];

% Matched termination at both ends
Gls1=0;
GgrN=0;

% Run the coefficients precomputation
tl = calc_tlines(z, Z0, k, Gls1, GgrN);

ztest = 5.01000000000000e-004;
ltest = 2;
vl    = 2;             % layer
ivd   = calc_ivd(tl, ztest, ltest, vl);
vvd   = calc_vvd(tl, ztest, ltest, vl);
iii   = calc_iii(tl, vl, ztest, ltest);

assertEquals(vvd,-iii,1e-25); % by reciprocity

n=1000000;
dl=tl.d(vl)/n;
zi=linspace(tl.z(vl), tl.z(vl+1), n);
test_ivd = sum(calc_iv(tl, ztest, ltest, zi, vl))*dl;
test_vvd = sum(calc_vv(tl, ztest, ltest, zi, vl))*dl;

assertEquals(test_ivd,ivd,1e-22);
assertEquals(test_vvd,vvd,1e-12);
