function Y = abcd2y(A);
% Y = abcd2y(A)
%
% ABCD to admittance transformation
% for 2xN-by-2xN matrices

N = size(A,1)/2;

A11=A(1:N,1:N);
A12=A(1:N,N+1:end);
A21=A(N+1:end,1:N);
A22=A(N+1:end,N+1:end);

invA12=inv(A12);

Y=zeros(size(A,1), size(A,2));
Y(1:N,1:N)         = A22*invA12;
Y(1:N,N+1:end)     = A21 - A22*invA12*A11;
Y(N+1:end,1:N)     = -invA12;
Y(N+1:end,N+1:end) = invA12*A11;
