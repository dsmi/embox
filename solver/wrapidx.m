function widx = wrapidx(idx, size)
% widx = wrapidx(idx, size)
%
% Wraps index in a cyclic fasion. Handles negative indices correctly,
% 0 becomes the last element, -2 second to last etc.
% 
%  wrapidx([ -2 -1 0 1 5 ], 4) => [ 2 3 4 1 1 ]
%
widx = idx - floor((idx-1)./size)*size;
