function layer = mklayer(B, BV, conduct)
% layer = mklayer(B, BV, conduct)
%
% Creates a vertical/flat mesh representing a layer of metallization.
% Given a bitmap representing the metal shape (recall that we use the
% rectilinear mesh - one consisting of axis-parallel rectangles/squares)
% where 0 means no metal and 1 means metal returns the set of basis functions
% which can be used to approximate the current distribution on this
% metallization layer. There are x- and y-directed basis functions,
% and the resulting structure has the following fields:
%   xi, xj - in-grid coordinates of the x-directed basis functions
%   yi, yj - in-grid coordinates of the y-directed basis functions
% In addition the layer may contain vias, which are the connections between
% this layer and the next one. If this function is given the second bitmap
% BV it creates vias from it. Although the vias can not penetrate the boundary
% and therefore the vias bitmap should be two pixels less in both X- and Y-
% dimensions, for the caller convenience this function expects both B and BV
% bitmaps to be of the same size and crops the BV by cuttint one pixel wide
% strips from each side.
% The corresponding fields of the resulting layer struct are:
%   vi, vj - in-grid coordinates of the via basis functions
% Each basis function is represented by coordinates of its supporting point,
% see drawing below where the supporting point is marked by 'o'
% 
%         +--+
% y       |  |   
% ^       o--+   +--+--+   +--+
% |       |  |   |  |  |   |  | <- this is via
% +-> x   +--+   +--o--+   o--+
%

% Identify x-directed basis functions
Bx = B+[ B(2:end,:) ; zeros(1, size(B, 2)) ];
[ xi, xj ] = find(Bx(:,2:end-1) > 1.5);

% Identify y-directed basis functions
By = B+[ B(:,2:end)  zeros(size(B, 1), 1) ];
[ yi, yj ] = find(By(2:end-1,:) > 1.5);

layer=struct('xi', xi-1, 'xj', xj-1, 'yi', yi-1, 'yj', yj-1);

if exist('BV')
    % strip edges - not needed for the vias
    [ vi, vj ] = find(BV(2:end-1, 2:end-1) > 0.5);
    layer.vi = vi - 1;
    layer.vj = vj - 1;
else
    % empty fields for the via
    layer.vi = ones(0,1);
    layer.vj = ones(0,1);
end

% Add conductivity if defined
if exist('conduct')
    layer.conductivity = conduct;
end

