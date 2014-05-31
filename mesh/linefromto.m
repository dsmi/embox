function B=linefromto(B0, x0, y0, x1, y1, th)
% B=linefromto(B0, x0, y0, x1, y1, th)
%
% Draws a line on a monochrome bitmap. Unlike drawline the linefromto
% uses normalized coordinates 0..1, but - the first and last pixels
% are not included in this unity range, i.e. the coordinates are
% recalculated as: x=xpx*(npx-2)+1 (where xpx is the coordinate
% in pixels and npx is the bitmap size in the corresponding dimension).
% This is because the first and last pixels are reserved for basis functions
% crossing the waveguide boundary.
%
% Inputs:
%    B0        - source bitmap
%    x0, y0    - starting point
%    x1, y1    - ending point
%    th        - thickness of the line
% Outputs:
%    B         - bitmap with the line drawn
%

% bitmap size
nx = size(B0, 1);
ny = size(B0, 2);

% coordinates and thickness in pixels
xpx0 = x0*(nx-2)+1;
ypx0 = y0*(ny-2)+1;
xpx1 = x1*(nx-2)+1;
ypx1 = y1*(ny-2)+1;
thpx = th*(nx-2);

B=drawline(B0, xpx0, ypx0, xpx1, ypx1, thpx);
