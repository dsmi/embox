function mesh = mkhull( tbp, edg, l0, n )
% mesh = mkhull( tbp, edg, l0, n )
%
% Makes a mesh of a hull or body, consisting of the top, bottom
% and edges.
%   tbp  top and bottom planes bitmap
%   edg  edges bitmap
%   l0   starting layer
%   n    number of layers

% Start from an empty structure array
mesh.layers = struct( [] );

% And add layers
for lidx = 1:(n+1)
    bh = max( tbp*( lidx == 1 || lidx == (n+1) ), edg ); % horizontal metal
    bv = edg*( lidx > 1 );                               % vias
    mesh.layers( end + 1 ) = mklayer( bh, bv, l0 - 1 + lidx, ccopper );
end
