function B=drawcir(B0, xc, yc, r)
% B=drawline(B0, xc, yc, r)
%
% Draws a filled circle on a monochrome bitmap. The coordinates and thickness
% are in pixels, but are not neccessarily integers.
%
% Inputs:
%    B0        - source bitmap
%    xc, yc    - center of the circle
%    r         - radius
% Outputs:
%    B         - bitmap with the circle drawn
%

% beginning and end points
c = [ xc yc ];

% bitmap size
nx = size(B0, 1);
ny = size(B0, 2);

% coordinates of the bitmap pixel centers
[ px, py ] = ndgrid((1:nx)-0.5, (1:ny)-0.5);
p = cat(3, px, py);

% pixel to center distance
d = sqrt(sum((p - repmat(shiftdim(c, -1), nx, ny)).^2, 3));

B = B0;

% Find the pixels closer than r
tol = 1e-10;
lpx = find(abs(d) <= r+tol);
B(lpx)=1;
