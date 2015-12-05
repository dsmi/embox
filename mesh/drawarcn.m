function B=drawarcn(B0, xc, yc, r, th, a0, a1)
% B=drawarcn(B0, xc, yc, r, th, a0, a1)
%
% Wrapper for the drawarc - draws an arc on the monochrome bitmap.
% Unlike drawarc uses the normalized coordinates 0..1, but - the first
% and last pixels are not included in this unity range, i.e. the coordinates
% are recalculated as: x=xpx*(npx-2)+1 (where xpx is the coordinate
% in pixels and npx is the bitmap size in the corresponding dimension).
% This is because the first and last pixels are reserved for basis functions
% crossing the waveguide boundary.
%
% Inputs:
%    B0        - source bitmap
%    xc, yc    - center of the arc
%    r         - radius of the arc
%    th        - thickness of the arc
%    a0, a1    - starting and ending angles in the range [-pi,pi]
% Outputs:
%    B         - bitmap with the arc drawn
%

% bitmap size
nx = size(B0, 1);
ny = size(B0, 2);

% coordinates and thickness in pixels
xcp = xc*(nx-2)+1;
ycp = yc*(ny-2)+1;
rp  = r*(nx-2);
thp = th*(nx-2);

B=drawarc(B0, xcp, ycp, rp, thp, a0, a1);
