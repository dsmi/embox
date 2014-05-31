function b = findbases(mesh, nx, ny, x0, y0, x1, y1)
% b = findbases(mesh, nx, ny, x0, y0, x1, y1)
%
% Finds basis functions, whether x-directed or y-directed, crossing the
% given line. Uses normalized coordinates (0..1). Follows the common agreement
% that x-directed bases come first and y-directed follow them.
% 

% mesh cell sizes
dx=1/nx;
dy=1/ny;

% tolerance used
tol = 1e-15;

% x-directed
if x0 == x1
    % coordinates of the basis centers
    xc=mesh.xi*dx;
    yc=mesh.xj*dy+dy/2;
    b = find(abs(xc-x0) < tol & yc >= y0-tol & yc <= y1+tol);
else % y-directed
    % coordinates of the basis centers
    xc=mesh.yi*dx+dx/2;
    yc=mesh.yj*dy;
    b = find(xc >= x0-tol & xc <= x1+tol & abs(yc-y0) < tol)+length(mesh.xi);
end
