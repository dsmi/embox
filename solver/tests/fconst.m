function y = fconst(x, x0, dx)
% y = fconst(x, x0, dx)
%
% Piecewise-constant base function given by
% y = 1    if (x0 - dx/2) < x < (x0 + dx/2)
%     0    otherwise
%

dx2 = dx/2;
y = ((x > (x0-dx2)) & (x < (x0+dx2)))*1.0;
