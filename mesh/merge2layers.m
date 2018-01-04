function layer = merge2layers(lay1, lay2)
% layer = merge2layers(lay1, lay2)
%
% Merge two layers together, duplicated basis functions in the merged
% layer are removed. For how to create layer please refer to mklayer.
%

layer = lay1;

layer.xi = [ lay1.xi ; lay2.xi ];
layer.xj = [ lay1.xj ; lay2.xj ];

layer.yi = [ lay1.yi ; lay2.yi ];
layer.yj = [ lay1.yj ; lay2.yj ];

layer.vi = [ lay1.vi ; lay2.vi ];
layer.vj = [ lay1.vj ; lay2.vj ];

% Remove duplicated basis functions. Attempt to preserve ordering
% as it comes from mklayer, so xjxi not xixj
xjxi  = [ layer.xj layer.xi ];
uxjxi = unique(xjxi, "rows");

if size(xjxi, 1) ~= size(uxjxi, 1)
    layer.xi = uxjxi(:, 2);
    layer.xj = uxjxi(:, 1);
end

yjyi  = [ layer.yj layer.yi ];
uyjyi = unique(yjyi, "rows");

if size(yjyi, 1) ~= size(uyjyi, 1)
    layer.yi = uyjyi(:, 2);
    layer.yj = uyjyi(:, 1);
end

vjvi  = [ layer.vj layer.vi ];
uvjvi = unique(vjvi, "rows");

if size(vjvi, 1) ~= size(uvjvi, 1)
    layer.vi = uvjvi(:, 2);
    layer.vj = uvjvi(:, 1);
end
