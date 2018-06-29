function b = findbases(mesh, nx, ny, x0, y0, x1, y1, flay)
% b = findbases(mesh, nx, ny, x0, y0, x1, y1, flay)
%
% Finds basis functions:
%   a) x-directed if or y-directed, crossing the given line, 
%         if x0 = x1 or y0 = y1.
%   b) z-directed (vias) in x0,y0 -- x1,y1 rectangle otherwise.
% Uses normalized coordinates (0..1). Follows the common agreement
% that x-directed bases come first and y-directed follow them, and
% the layers are processed in the order as they follow in the mesh
% structure.
% Finds the basis functions on all layers unless lay is specified.
% 

% mesh cell sizes
dx=1/nx;
dy=1/ny;

% tolerance used
tol = 1e-15;

% total number of basis functions on the already processed layers
onprev = 0;

% start with empty, process layers one by one
b = [];

for lidx = 1:length(mesh.layers)

    layer = mesh.layers(lidx);

    % see if user specified a particular layer
    if ~exist('flay') || flay(layer.pos)

	% x-directed
	if x0 == x1
	    % coordinates of the basis centers
	    xc=layer.xi*dx;
	    yc=layer.xj*dy+dy/2;
	    fb = find(abs(xc-x0) < tol & yc >= y0-tol & yc <= y1+tol);
	elseif y0 == y1 % y-directed
	    % coordinates of the basis centers
	    xc=layer.yi*dx+dx/2;
	    yc=layer.yj*dy;
	    fb = find(xc >= x0-tol & xc <= x1+tol & abs(yc-y0) < tol);
	    fb = fb + length(layer.xi); % x-directed come first
	else            
	    % coordinates of the basis centers
	    xc = layer.vi*dx + dx/2;
	    yc = layer.vj*dy + dy/2;
	    fb = find(xc >= x0-tol & xc <= x1+tol & yc >= y0-tol & yc <= y1+tol);
	    fb = fb + length(layer.xi) + length(layer.yi); % x- and y-directed come first
	end

	b = [ b ; (fb + onprev) ];

    end
	
    onprev = onprev + length(layer.xi) + length(layer.yi) + length(layer.vi);

end
