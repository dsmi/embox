function B=drawline(B0, x0, y0, x1, y1, th)
% B=drawline(B0, x0, y0, x1, y1, th)
%
% Draws a line on a monochrome bitmap. The coordinates and thickness are in
% pixels, but are not neccessarily integers.
%
% Inputs:
%    B0        - source bitmap
%    y0, y0    - starting point
%    y0, y0    - ending point
%    th        - thickness of the line
% Outputs:
%    B         - bitmap with the line drawn
%

% beginning and end points
p0 = [ x0 y0 ];
p1 = [ x1 y1 ];

% tangent vector, non-normalized
tg = p1 - p0;

% squared lenght
l2 = sum(tg.*tg);

% bitmap size
nx = size(B0, 1);
ny = size(B0, 2);

% coordinates of the bitmap pixel centers
[ px, py ] = ndgrid((1:nx)-0.5, (1:ny)-0.5);
p = cat(3, px, py);

% The line the segment lies on is normalized as p0+t*tg
% This finds t for all the pixels
t = dot(p - repmat(shiftdim(p0, -1), nx, ny), repmat(shiftdim(tg, -1), nx, ny), 3) ./ l2;

% Projection on the line
proj = repmat(shiftdim(p0, -1), nx, ny) + repmat(t, [ 1 1 2 ]).*repmat(shiftdim(tg, -1), nx, ny);

% Distance from the line
r = sqrt(sum((p-proj).^2, 3));

B = B0;

% Find the conforming pixels
tol = 1e-10;
lpx = find(abs(r) <= th/2+tol & t >= -tol & t <= 1+tol);
B(lpx)=1;
