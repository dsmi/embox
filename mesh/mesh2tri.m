function [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I, fi2c)
% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I, fi2c)
% 'Triangulates' mesh so it can be drawn with trisurf/trimesh
%
% wg     - shiedling parameters, see wgparams
% mesh   - meshed metal, see mkmesh
% I      - currents
% fi2c   - function to use to get color from current phasor, real/abs etc
%

% mesh cell sizes
dx = wg.a/wg.nx;
dy = wg.b/wg.ny;

% z coordinate of the layers
lz = [ 0 cumsum(wg.h) ];

% Cumulative sums of numbers of the basis functions in the layers up to
% the given one which is then used to find currents in the vector
cumx = cumsum(cellfun(@(v) length(v), { mesh.layers(:).xi }));
cumy = cumsum(cellfun(@(v) length(v), { mesh.layers(:).yi }));
cumv = cumsum(cellfun(@(v) length(v), { mesh.layers(:).vi }));
cumbf = [ 0 (cumx + cumy + cumv) ];

xyz = [];
Tri = [];
C = [];

if ~exist('fi2c')
   fi2c = @real;
end


for lidx = 1:length(mesh.layers)

    lay = mesh.layers(lidx);
    B = zeros(wg.nx+2,wg.ny+2);

    % Number of x-, y- and z-directed basis functions for this layer
    Nx = length(lay.xi);
    Ny = length(lay.yi);
    Nz = length(lay.vi);

    % currents of this layer
    Ix = I( cumbf(lidx)           + 1 : cumbf(lidx) + Nx );
    Iy = I( cumbf(lidx) + Nx      + 1 : cumbf(lidx) + Nx + Ny );
    Iz = I( cumbf(lidx) + Nx + Ny + 1 : cumbf(lidx) + Nx + Ny + Nz );

    % indices of the basis functions in B
    xind1 = sub2ind(size(B), lay.xi+1, lay.xj+2); % source
    xind2 = sub2ind(size(B), lay.xi+2, lay.xj+2); % sink
    yind1 = sub2ind(size(B), lay.yi+2, lay.yj+1); % source
    yind2 = sub2ind(size(B), lay.yi+2, lay.yj+2); % sink

    % populate bitmap - ones correspond to basis functions
    B(xind1) = 1;
    B(xind2) = 1;
    B(yind1) = 1;
    B(yind2) = 1;

    % find nonzero pixels in B, and generate verices
    [ ii, jj ] = find(B);
    nr = length(ii); % number of rectangles
    iijj = [ kron(ii, ones(4,1)) , kron(jj, ones(4,1)) ];
    dxdy = repmat([ dx dy ], nr*4, 1);
    % 2 subtracted to get coordinates in ranges [ -dx a+dx ] [ -dy a+dy ] 
    xy = (iijj + repmat([ 0 0 ; 0 1 ; 1 1 ; 1 0 ] - 2, nr, 1)).*dxdy;

    % Color for each vertex of the rect is computed indivitually. Color of
    % vertex1 is sum of incoming x and y currents since this vertex has
    % 'incoming' edges adjacent to it, vertex2 has incoming x and outgoing
    % y edges and so on.
    color1 = color2 = color3 = color4 = B*0;
    color1(xind2) = color1(xind2) + fi2c(Ix);
    color1(yind2) = color1(yind2) + fi2c(Iy);
    color2(xind2) = color2(xind2) + fi2c(Ix);
    color2(yind1) = color2(yind1) + fi2c(Iy);
    color3(xind1) = color3(xind1) + fi2c(Ix);
    color3(yind1) = color3(yind1) + fi2c(Iy);
    color4(xind1) = color4(xind1) + fi2c(Ix);
    color4(yind2) = color4(yind2) + fi2c(Iy);

    % Interleave the vertex colors
    colors = xy(:,1)*0;
    colors(1:4:end) = color1(sub2ind(size(B), ii, jj));
    colors(2:4:end) = color2(sub2ind(size(B), ii, jj));
    colors(3:4:end) = color3(sub2ind(size(B), ii, jj));
    colors(4:4:end) = color4(sub2ind(size(B), ii, jj));

    % and then triangles
    voffs = repmat(kron((0:4:(nr*4-4))', ones(2,1)), 1, 3); % vertex index offsets
    T = repmat([ 1 2 3 ; 1 3 4 ], nr, 1) + voffs;

    % now add triangulation of this layer to the overall
    if numel(xy)
	Tri = [ Tri ; T + size(xyz,1) ];
	xyz = [ xyz ; [ xy xy(:,1)*0 + lz(lay.pos + 1) ] ];
	C   = [ C   ; colors ];
    end

    % Via meshing. One via per cell so no bitmap needed
    ii = lay.vi(:) + 2; % to make it similar to the above
    jj = lay.vj(:) + 2;
    nr = length(ii); % number of rectangles
    iijj = [ kron(ii, ones(8,1)) , kron(jj, ones(8,1)) ];
    % 2 subtracted to get coordinates in ranges [ -dx a+dx ] [ -dy a+dy ] 
    vivj = iijj + repmat([ 0 0 ; 0 1 ; 1 1 ; 1 0 ] - 2, nr*2, 1);
    dxdy = repmat([ dx dy ], nr*8, 1);
    xy = vivj.*dxdy;
    z1 = lz(lay.pos);
    z2 = lz(lay.pos + 1);
    z = repmat([ repmat(z1, 4, 1) ; repmat(z2, 4, 1) ], nr, 1);
    
    % triangles for the via
    edget = [ 1 2 6 ; 1 6 5 ];
    viat = [ edget ; edget + 1 ; edget + 2 ; [ 4 1 5 ; 4 5 8 ] ];
    voffs = repmat(kron((0:8:(nr*8-8))', ones(8,1)), 1, 3); % vertex index offsets
    T = repmat(viat, nr, 1) + voffs;

    % color from current
    color = fi2c(Iz);

    % add triangulation of this layer vias to the overall
    if numel(xy)
    	Tri = [ Tri ; T + size(xyz,1) ];
    	xyz = [ xyz ; [ xy z ] ];
	C   = [ C   ; kron(color, ones(8,1)) ];
    end

end
 
X = xyz(:,1);
Y = xyz(:,2);
Z = xyz(:,3);
