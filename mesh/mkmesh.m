function mesh = mkmesh(B)
% mesh = mkmesh(B)
%
% Given a bitmap representing the metal shape (recall that we use the
% rectilinear mesh - one consisting of axis-parallel rectangles/squares)
% where 0 means no metal and 1 means metal returns the set of basis functions
% which can be used to approximate the current distribution. There are x- and
% y-directed basis functions, and the resulting structure has the following
% fields:
%   xi, xj - in-grid coordinates of the x-directed basis functions
%   yi, yj - in-grid coordinates of the y-directed basis functions
% Each basis function is represented by coordinates of its supporting point,
% see drawing below where the supporting point is marked by 'o'
% 
%         +--+
% y       |  |   
% ^       o--+   +--+--+
% |       |  |   |  |  |
% +-> x   +--+   +--o--+
%

% Identify x-directed basis functions
Bx = B+[ B(2:end,:) ; zeros(1, size(B, 2)) ];
[ xi, xj ] = find(Bx(:,2:end-1) > 1.5);

% Identify y-directed basis functions
By = B+[ B(:,2:end)  zeros(size(B, 1), 1) ];
[ yi, yj ] = find(By(2:end-1,:) > 1.5);

mesh=struct('xi', xi-1, 'xj', xj-1, 'yi', yi-1, 'yj', yj-1);
