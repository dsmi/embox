function test_y2abcd()

N=3;

Y=rand(N*2,N*2) + j*rand(N*2,N*2); % bad idea to use rand?

A=y2abcd(Y);
v=rand(N*2,1);
i=Y*v;

v1 = v(1:N);
v2 = v(N+1:end);
i1 = i(1:N);
i2 = i(N+1:end);

v1i1 = A*[ v2 ; -i2 ]; % inverted sign of the current
test_v1i1 = [ v1 ; i1 ];

assertEquals(test_v1i1, v1i1, 1e-14);


