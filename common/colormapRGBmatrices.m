function mymap = colormapRGBmatrices( N, rm, gm, bm)
% mymap = colormapRGBmatrices( N, rm, gm, bm)
%
% It expects 4 input parameters: N is the number of intermediate points that
% your colormap should have. The other three are matrices that contain
% the transitions for each channel. Such a matrix should have the following
% form:
% M = [ 0, b1;
%       x1,b2;
%       ...
%       x_n, b_(n+1);
%       1, b_(n+1);
%     ];
%  the first column give the fractions, where a brightness value should be
%  defined. The second column should contain the brightness levels. Make
%  sure to start the first column at 0 and end it with 1!
%  A simple, linear grayscale map can be created with
%   M = [0,0;1,1;];
%   simplegray = colormapRGBmatrices( 256, M, M, M);
%
% Taken from http://cresspahl.blogspot.com/2012/03/expanded-control-of-octaves-colormap.html
%

  x = linspace(0,1, N);
  rv = interp1( rm(:,1), rm(:,2), x);
  gv = interp1( gm(:,1), gm(:,2), x);
  mv = interp1( bm(:,1), bm(:,2), x);
  mymap = [ rv', gv', mv'];
  %exclude invalid values that could appear
  mymap( isnan(mymap) ) = 0;
  mymap( (mymap>1) ) = 1;
  mymap( (mymap<0) ) = 0;
end
