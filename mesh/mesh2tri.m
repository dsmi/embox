function [ Tri, X, Y, Z ] = mesh2tri(wg, mesh)
% [ Tri, X, Y, Z ] = mesh2tri(wg, mesh)
% 'Triangulates' mesh so it can be drawn with trisurf/trimesh
%
% wg     - shiedling parameters, see wgparams
% mesh   - meshed metal, see mkmesh
% 

% mesh cell sizes
dx = wg.a/wg.nx;
dy = wg.b/wg.ny;

% z coordinate of the layers
lz = [ 0 cumsum(wg.h) ];

xyz = [];
Tri = [];

for lidx = 1:length(mesh.layers)

    lay = mesh.layers(lidx);
    B=zeros(wg.nx+2,wg.ny+2);

    % populate bitmap - ones corresponds to basis functions
    B(sub2ind(size(B), lay.xi+1,lay.xj+2)) = 1;
    B(sub2ind(size(B), lay.xi+2,lay.xj+2)) = 1;
    B(sub2ind(size(B), lay.yi+2,lay.yj+1)) = 1;
    B(sub2ind(size(B), lay.yi+2,lay.yj+2)) = 1;

    % find nonzero pixels in B, and generate verices
    [ ii, jj ] = find(B);
    nr = length(ii); % number of rectangles
    iijj = [ kron(ii, ones(4,1)) , kron(jj, ones(4,1)) ];
    dxdy = repmat([ dx dy ], nr*4, 1);
    xy = (iijj + repmat([ 0 0 ; 0 1 ; 1 1 ; 1 0 ], nr, 1)).*dxdy;

    % and then triangles
    voffs = repmat(kron((0:4:(nr*4-4))', ones(2,1)), 1, 3); % vertex index offsets
    T = repmat([ 1 2 3 ; 1 3 4 ], nr, 1) + voffs;

    % now add triangulation of this layer to the overall
    if numel(xy)
	Tri = [ Tri ; T + size(xyz,1) ];
	xyz = [ xyz ; [ xy xy(:,1)*0 + lz(lay.pos + 1) ] ];
    end

    % now the vias in this layer
    B=zeros(wg.nx+2,wg.ny+2);
    B(sub2ind(size(B), lay.vi+2,lay.vj+2)) = 1;

    % similar to the above - generate vertices for nonzero pixels of B,
    % but now we generate two quads one on top of another
    [ ii, jj ] = find(B);
    nr = length(ii); % number of rectangles
    iijj = [ kron(ii, ones(8,1)) , kron(jj, ones(8,1)) ];
    vivj = iijj + repmat([ 0 0 ; 0 1 ; 1 1 ; 1 0 ], nr*2, 1);
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

    % add triangulation of this layer vias to the overall
    if numel(xy)
    	Tri = [ Tri ; T + size(xyz,1) ];
    	xyz = [ xyz ; [ xy z ] ];
    end

end
 
X = xyz(:,1);
Y = xyz(:,2);
Z = xyz(:,3);
