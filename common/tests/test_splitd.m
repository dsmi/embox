function test_splitd()

% test for 2-by-2 first
C11=10e-10;
w = 1.0e3*2*pi;
Y11 = j*w*C11;

D2 = [ 1 0 ; Y11 1 ];

[ D, tol ] = splitd(D2);

assertEquals(D2, D*D, 1e-15)


% now for 6-by-6
C11=10e-10;
C22=20e-10;
C33=30e-10;
C12=1e-10;
C13=2e-10;
C23=4e-10;

w = 1.0e3*2*pi;
Y11 = j*w*C11;
Y22 = j*w*C22;
Y33 = j*w*C33;
Y12 = j*w*C12;
Y13 = j*w*C13;
Y23 = j*w*C23;

N = 3;

Y2 = [ (Y11+Y12+Y13) -Y12 -Y13 ; -Y12 (Y22+Y12+Y23) -Y23 ; -Y13 -Y23 (Y33+Y13+Y23) ];
D2 = [ eye(N) zeros(N) ; Y2 eye(N) ];

[ D, tol ] = splitd(D2);

assertEquals(D2, D*D, 1e-15)

% check how tol works

D2t = D2;
D2t(1:N,1:N) = D2t(1:N,1:N) + ones(N,N)*1e-10;
[ D, tol ] = splitd(D2t);
assertEquals(3e-10, tol, 1e-15)

D2t = D2;
D2t(N+1:end,N+1:end) = D2t(N+1:end,N+1:end) + ones(N,N)*1e-10;
[ D, tol ] = splitd(D2t);
assertEquals(3e-10, tol, 1e-15)

D2t = D2;
D2t(1:N,N+1:end) = D2t(1:N,N+1:end) + ones(N,N)*1e-10/norm(Y2);
[ D, tol ] = splitd(D2t);
assertEquals(1.5e-10, tol, 1e-15)
