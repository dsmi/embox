function test_abcd2y()

N=3;

A=rand(N*2,N*2) + j*rand(N*2,N*2); % bad idea to use rand?

v2=rand(N,1);
i2=rand(N,1);
v1i1=A*[ v2 ; -i2 ];
v1 = v1i1(1:N);
i1 = v1i1(N+1:end);

Y = abcd2y(A);

i = Y*[ v1 ; v2 ];
test_i = [ i1 ; i2 ];

assertEquals(test_i, i, 1e-14);

% This tests abcd2y against y2abcd
A2 = y2abcd(Y);
assertEquals(A, A2, 1e-14);
