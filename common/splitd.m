function [ D, tol ] = splitd(D2)
% [ D, tol ] = splitd(D2)
%
% Performs one of the stages of the deembedding process - obtains the port
% discontinuity matrix D from the double discontinuity matrix which is the 
% discontinuity matrix cascaded with itself or D*D since we work with ABCD
% parameters. This only works for certain kind of discontinuities which can
% be represented as a shunt admittance per port plus mutual admittances
% between the ports. In this case the ABCD matrix has the following form
%   I  0
%   Y  I
% Where I is the identity matrix and Y is
%  Y11 + Y12 + .. + Y1N         -Y12                  -Y13
%         -Y21           Y21 + Y22 + .. + Y2N         -Y23
%         -Y31                  -Y32            Y21 + Y22 + .. + Y2N
% and Ymn = Ynm.
% The double discontinuity matrix can be obtained from simulations of
% transmission lines of lengths l and 2*l as
%   A1*inv(A2)*A1
% Where A1 and A2 are the ABCD matrices from l and 2*l simulations
% correspondingly.
% The double discontinuity matrix is expected to have some form because
% of the assumptions about the discontinuity we made, and the tol returned
% value shows how much the D2 matrix matches the expected form.
%

N = size(D2,1)/2; % number of ports

if N*2 ~= size(D2,1)
    error('Odd number of ports')
end

Y2 = D2(N+1:end,1:N);

% Off-diagonal admittances
Y2o = -tril(Y2,-1);

% Diagonal admittances can be obtained as: 
% Y2d = diag(Y2) - sum(Y2o, 2);

% Divide by 2 and symmetrize
Y = diag(diag(Y2))*0.5 - Y2o*0.5 - transpose(Y2o)*0.5;

% The resulting port discontinuity matrix
D = [ eye(N) zeros(N) ; Y eye(N) ];

% estimate the error
tol = norm(D2(1:N,1:N) - eye(N)) + norm(D2(N+1:end,N+1:end) - eye(N));

% difference between the given D2 and the resulting D doubled (ideally should
% be zero)
Ddiff = D2 - D*D;

% The norm of the Y block is used to 'normalize' the errors in the Y and Z
% blocks of the ABCD matrix
normY = norm(Y);

% Finally, include the Y and Z blocks errors normalized using the Y block
% norm in the error estimate
tol = tol + norm(Ddiff(N+1:end,1:N))/normY + norm(D2(1:N,N+1:end))*normY;
