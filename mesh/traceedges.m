function B = traceedges(B0)
% B=drawline(B0, x0, y0, x1, y1, th)
%
% Traces image edges on a monochrome bitmap and returns another bitmap
% with edges only. Edges are inside the shape. Used to place vias when
% they are to be placed at the edges of the conductor shapes.
%

% A row and a column of zeros
row0 = zeros(1, size(B0, 2));
col0 = zeros(size(B0, 1), 1);

% Summed with offset images
bs = B0 + [ B0(2:end,:) ; row0 ] + [ row0 ; B0(1:(end-1),:) ] ...
        + [ B0(:,2:end) , col0 ] + [ col0, B0(:,1:(end-1)) ];

% Find the conforming pixels
B = B0*0;
B(find(B0 > 0 & bs < 5)) = 1;
