function y = ftri(x, x0, dx)
% y = ftri(x, x0, dx)
%
% Triangular basis function give by
%     (x-x0)/dx+1   if (x0-dx) < x < x0
% y = (x0-x)/dx+1   if x0 < x < (x0 + dx)
%     0  otherwise
%

y1 = ((x > (x0-dx)) & (x < x0)).*((x-x0)/dx+1);
y2 = ((x > x0) & (x < (x0+dx))).*((x0-x)/dx+1);
y = y1 + y2;
