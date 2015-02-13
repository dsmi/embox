function mesh = mkmesh(varargin)
% mesh = mkmesh(layer1,pos1,...,layerN,pos1)
%
% Composes a mesh structure.
% In this solver the mesh consists of a number of vertical/flat metal layers,
% and this function adds the position field to each layer (which indicates
% the layer position in the stackup) and composes the mesh structure which
% carries the layers.
% The mesh structure currently has only one field - mesh.layers - which is
% array of structures represening the layer.
% See mklayer for a possible way of creating a layer.
% 

for i=1:nargin/2
    layer = varargin{2*i-1};
    layer.pos = varargin{2*i};
    mesh.layers(i) = layer;
end
