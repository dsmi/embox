function S1 = renorms(S0, Z0, Z1)
% S1 = renorms(S0, Z0, Z1)
%
% Renormalize S-parameters
%

T = diag(sqrt(Z1./Z0));
Q = diag(sqrt(Z0./Z1));
D = diag(Z0./Z1);

for i=1:size(S0,3)
    S = S0(:,:,i);
    I = eye(size(S));
    S1(:,:,i) = Q*inv( (D+I) + S*(D-I) )*( (D-I) + S*(D+I) )*T;
end
