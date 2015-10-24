function test_linsteps
%

ends = [ 0 1 5 ];
N = [ 3 5 ];
x = linsteps(ends, N);
xtest = [ 0 0.5 1 2 3 4 5 ];

assertEquals(xtest, x, 1e-12);
