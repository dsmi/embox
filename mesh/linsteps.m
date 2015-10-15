function x = linsteps(ends, N)
% x = linsteps(ends, N)
%
% Subdivides adjacent intervals with linearly spaced points.
% i'th interval is given by [ ends(i) ends(i+1) ], number of points for
% this interval is N(i). Endpoints are not duplicated.
%

%x=ends;

num = 1:length(ends)-1;
lns = @(i) linspace(ends(i), ends(i+1), N(i));
nolast = @(a) a(1:end-1);
x = cell2mat(arrayfun( @(i) nolast(lns(i)), num, 'UniformOutput', false));
x = [ x ends(end) ];
