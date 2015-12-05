function B=drawarc(B0, xc, yc, r, th, a0, a1)
% B=drawarc(B0, xc, yc, r, th, a0, a1)
%
% Draws an arc on a monochrome bitmap. The coordinates and thickness
% are in pixels, but are not neccessarily integers.
%
% Inputs:
%    B0        - source bitmap
%    xc, yc    - center of the arc
%    r         - radius of the arc
%    th        - thickness of the arc
%    a0, a1    - starting and ending angles in the range [-pi,pi]
% Outputs:
%    B         - bitmap with the arc drawn
%

% beginning and end points
c = [ xc yc ];

% bitmap size
nx = size(B0, 1);
ny = size(B0, 2);

% coordinates of the bitmap pixel centers
[ px, py ] = ndgrid((1:nx)-0.5, (1:ny)-0.5);
p = cat(3, px, py);

% pixel to center vectors
rv = p - repmat(shiftdim(c, -1), nx, ny);

% pixel to center distance
d = sqrt(sum(rv.^2, 3));

% pixel angles
a = atan2(rv(:,:,2), rv(:,:,1));

B = B0;

% Find the conforming pixels
tol = 1e-10;
lpx = find(d <= r+th/2+tol & d >= r-th/2-tol & a <= a1 + tol & a >= a0 - tol);
B(lpx)=1;
