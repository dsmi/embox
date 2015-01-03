function Z=mkzsmat(wg, mesh)
% Z=mkzsmat(wg, mesh)
%
% Populates the surface impedance matrix, which is to be added to the
% PEC impedance matrix to take into account the finite conductance.
% Elements of the matrix are defined as:
%  Zmn = <Fn, Wm>
% where the Fn and Wm are the expansion and testing functions correspondingly
% and <,> is the scalar product.
% The matrix Zmn is then to be multiplied by the surface impedance.
%
% wg     - shiedling parameters, see wgparams
% mesh   - meshed metal, see mkmesh
% 

[ mxi, nxi ] = ndgrid(mesh.xi, mesh.xi);
[ mxj, nxj ] = ndgrid(mesh.xj, mesh.xj);

% area of one cell - used in the integrals
acell = (wg.a/wg.nx)*(wg.b/wg.ny);

% x-directed basis functions
ho = (abs(mxi - nxi) == 1)*(1/6); % half overlap case
fo = (mxi == nxi).*(1/3 + 1/3*((mxi > 0) & (mxi < wg.nx))); % full overlap
Zxx = acell*(mxj==nxj).*(ho+fo);

[ myi, nyi ] = ndgrid(mesh.yi, mesh.yi);
[ myj, nyj ] = ndgrid(mesh.yj, mesh.yj);

% y-directed basis functions
ho = (abs(myj - nyj) == 1)*(1/6); % half overlap case
fo = (myj == nyj).*(1/3 + 1/3*((myj > 0) & (myj < wg.ny))); % full overlap
Zyy = acell*(myi==nyi).*(ho+fo);

% compose the entire matrix
Zxy = zeros(size(Zxx,1), size(Zyy,2));
Zyx = zeros(size(Zyy,1), size(Zxx,2));
Z = [ Zxx Zxy ; Zyx Zyy ];
