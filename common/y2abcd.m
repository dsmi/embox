function A = y2abcd(Y);
% A = y2abcd(Y)
%
% Admittance to ABCD transformation
% for 2xN-by-2xN matrices

N = size(Y,1)/2;

Y11=Y(1:N,1:N);
Y12=Y(1:N,N+1:end);
Y21=Y(N+1:end,1:N);
Y22=Y(N+1:end,N+1:end);

invY21 = inv(Y21);

A=zeros(size(Y,1), size(Y,2));
A(1:N,1:N)         = -invY21*Y22;
A(1:N,N+1:end)     = -inv(Y21);
A(N+1:end,1:N)     = Y12 - Y11*inv(Y21)*Y22;
A(N+1:end,N+1:end) = -Y11*invY21;
