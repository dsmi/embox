function mesh = mergelayers(srcmesh)
% mesh = mergelayers(srcmesh)
%
% Merge layers with identical positions (pos). Layers are merged using
% merge2layers which removes duplicated basis functions.
%

% Source layer positions. Duplicated needs to be merged.
pos_cell = { srcmesh.layers(:).pos };
pos = [ pos_cell{:} ]; 

% this determines unique layers
[ upos, isrc, iunq ] = unique(pos);

% Fill unique layers
mesh.layers = srcmesh.layers(isrc);

% And merge
for li = 1:length(iunq)
    iun = iunq(li); % index in unique layers
    mesh.layers( iun ) = merge2layers( mesh.layers( iun ), srcmesh.layers( li ) );
end
