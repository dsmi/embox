function g=gtri(dx,kx)
% g=gtri(dx,kx)
% 
% Calculates constant resulting from integration of the triangular basis
% function.
% The target integrals are:
%  F=int_(x0-dx)^(x0+dx)(f(x)*sin(kx*x))dx
%  F=int_(x0-dx)^(x0+dx)(f(x)*cos(kx*x))dx
% Where f(x) is the basis function defined as
%  f(x) = (x-x0)/dx+1   if x0-dx <= x < x0
%         (x0-x)/dx+1   if x0    <= x < x0+dx
%         0             otherwise
% The integration result is
%   F=g(dx,kx)*sin(kx*x0)
%   F=g(dx,kx)*cos(kx*x0)
% Where g(dx,kx) is the constant returned by this function.
%
zidx=find(abs(kx)<1e-300); % indices of the kx=0
kx(zidx)=1e-30; % to avoid division by zero

g=2./(dx*kx.*kx).*(1-cos(dx*kx));
g(zidx)=dx; % special case kx=0

