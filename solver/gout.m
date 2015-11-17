function g=gout(dx,kx)
% g=gout(dx,kx)
% 
% Calculates constant resulting from integration of the 'out' basis function
% function - one which represent current leaving/entering a cell and used as
% the top and bottom of the via.
% The target integrals are:
%  F=int_(x0-dx/2)^(x0+dx/2)(f(x)*sin(kx*x))dx
%  F=int_(x0-dx/2)^(x0+dx/2)(f(x)*cos(kx*x))dx
% Where f(x) is the basis function defined as
%   f(x) = (x0-x)/dx   if (x0 - dx/2) < x < (x0 + dx/2)
%          0           otherwise
% The integration result is (notice sin and cos swapped and minus sign)
%   F=g(dx,kx)*cos(kx*x0)
%   F=-g(dx,kx)*sin(kx*x0)
% Where g(dx,kx) is the constant returned by this function.
%
zidx=find(abs(kx)<1e-300); % indices of the kx=0
kx(zidx)=1e-30; % to avoid division by zero

dx2 = dx/2;
g = -2.*(dx2.*kx.*cos(dx2.*kx)-sin(dx2.*kx))./(dx2*kx.*kx);
g(zidx)=0; % special case kx=0
