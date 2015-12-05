function [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I)
% [ Tri, X, Y, Z, C ] = mesh2tri(wg, mesh, I)
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

% Cumulative sums of numbers of the basis functions in the layers up to
% the given one which is then used to find currents in the vector
cumx = cumsum(cellfun(@(v) length(v), { mesh.layers(:).xi }));
cumy = cumsum(cellfun(@(v) length(v), { mesh.layers(:).yi }));
cumv = cumsum(cellfun(@(v) length(v), { mesh.layers(:).vi }));
cumbf = [ 0 (cumx + cumy + cumv) ];

xyz = [];
Tri = [];
C = [];

for lidx = 1:length(mesh.layers)

    lay = mesh.layers(lidx);
    B = zeros(wg.nx+2,wg.ny+2);

    % currents of this layer
    Il = I(cumbf(lidx)+1:cumbf(lidx+1));

    % horizontal currents
    ix = zeros(wg.nx+2,wg.ny+2);
    iy = zeros(wg.nx+2,wg.ny+2);

    % indices of the basis functions in B, ix and iy
    xind1 = sub2ind(size(B), lay.xi+1,lay.xj+2);
    xind2 = sub2ind(size(B), lay.xi+2,lay.xj+2);
    yind1 = sub2ind(size(B), lay.yi+2,lay.yj+1);
    yind2 = sub2ind(size(B), lay.yi+2,lay.yj+2);

    % populate bitmap - ones correspond to basis functions
    B(xind1) = 1;
    B(xind2) = 1;
    B(yind1) = 1;
    B(yind2) = 1;

    % calculate currents
    ix(xind1) = ix(xind1) + Il(1:length(xind1)).' ./ dx;
    ix(xind2) = ix(xind2) + Il(1:length(xind1)).' ./ dx;
    iy(yind1) = iy(yind1) + Il(length(xind1)+1:length(xind1)+length(yind1)).' ./ dy;
    iy(yind2) = iy(yind2) + Il(length(xind1)+1:length(xind1)+length(yind1)).' ./ dy;

    % find nonzero pixels in B, and generate verices
    [ ii, jj ] = find(B);
    nr = length(ii); % number of rectangles
    iijj = [ kron(ii, ones(4,1)) , kron(jj, ones(4,1)) ];
    dxdy = repmat([ dx dy ], nr*4, 1);
    % 2 subtracted to get coordinates in ranges [ -dx a+dx ] [ -dy a+dy ] 
    xy = (iijj + repmat([ 0 0 ; 0 1 ; 1 1 ; 1 0 ] - 2, nr, 1)).*dxdy;

    % magnitute of the current at t=0
    im = sqrt(real(ix).^2 + real(iy).^2);
    color = im(sub2ind(size(B), ii, jj));

    % and then triangles
    voffs = repmat(kron((0:4:(nr*4-4))', ones(2,1)), 1, 3); % vertex index offsets
    T = repmat([ 1 2 3 ; 1 3 4 ], nr, 1) + voffs;

    % now add triangulation of this layer to the overall
    if numel(xy)
	Tri = [ Tri ; T + size(xyz,1) ];
	xyz = [ xyz ; [ xy xy(:,1)*0 + lz(lay.pos + 1) ] ];
	C   = [ C   ; kron(color, ones(4,1)) ];
    end

    % now the vias in this layer
    B=zeros(wg.nx+2,wg.ny+2);
    B(sub2ind(size(B), lay.vi+2,lay.vj+2)) = 1;

    % similar to the above - generate vertices for nonzero pixels of B,
    % but now we generate two quads one on top of another
    [ ii, jj ] = find(B);
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

    % add triangulation of this layer vias to the overall
    if numel(xy)
    	Tri = [ Tri ; T + size(xyz,1) ];
    	xyz = [ xyz ; [ xy z ] ];
	C   = [ C   ; xyz(:,1)*0 ]; % ignored now
    end

end
 
X = xyz(:,1);
Y = xyz(:,2);
Z = xyz(:,3);
