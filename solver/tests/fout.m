function y = fout(x, x0, dx)
% y = fout(x, x0, dx)
%
% 'Out' basis function representing the current leaving the cell in +x and -x
% directions, used as the bottom or top (with negative sign) of the via.
% Is given by
% y = (x0-x)/(dx/2)  if (x0 - dx/2) < x < (x0 + dx/2)
%     0              otherwise
%

dx2 = dx/2;
y = ((x > (x0-dx2)) & (x < (x0+dx2))).*((x-x0)/dx2);
