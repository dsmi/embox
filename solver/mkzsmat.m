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

% We want to pre-allocate the Z matrix, for that we need to know its size.
% To calculate the size we need to count the basis functions on all layers.
numx = sum(cellfun(@(v) length(v), { mesh.layers(:).xi }));
numy = sum(cellfun(@(v) length(v), { mesh.layers(:).yi }));
numv = sum(cellfun(@(v) length(v), { mesh.layers(:).vi }));
numbf = numx + numy + numv;

% We also compute the cumulative sums of numbers of the basis functions in
% the layers up to the given one which is then used to place the blocks
% of the matrix (see below)
cumx = cumsum(cellfun(@(v) length(v), { mesh.layers(:).xi }));
cumy = cumsum(cellfun(@(v) length(v), { mesh.layers(:).yi }));
cumv = cumsum(cellfun(@(v) length(v), { mesh.layers(:).vi }));
cumbf = [ 0 (cumx + cumy + cumv) ];

% Pre-allocate it!
Z = sparse(numbf,numbf);

% Here we start popolating the impedance matrix. The geometry consists of a
% number of layers of metallization (vias are on their way) and the Z matrix
% consists of N-by-N blocks where N is the number of layers. The block M(m,n)
% corresponds to the m-th observation and n-th source layer - and in case of
% the surface impedance it only has nonzeros if m == n.
for lidx = 1:length(mesh.layers)

    layer = mesh.layers(lidx);

    [ mxi, nxi ] = ndgrid(layer.xi, layer.xi);
    [ mxj, nxj ] = ndgrid(layer.xj, layer.xj);

    % area of one cell - used in the integrals
    acell = (wg.a/wg.nx)*(wg.b/wg.ny);

    % x-directed basis functions
    ho = (abs(mxi - nxi) == 1)*(1/6); % half overlap case
    fo = (mxi == nxi).*(1/3 + 1/3*((mxi > 0) & (mxi < wg.nx))); % full overlap
    Zxx = acell*sparse( (mxj==nxj).*(ho+fo) );

    [ myi, nyi ] = ndgrid(layer.yi, layer.yi);
    [ myj, nyj ] = ndgrid(layer.yj, layer.yj);

    % y-directed basis functions
    ho = (abs(myj - nyj) == 1)*(1/6); % half overlap case
    fo = (myj == nyj).*(1/3 + 1/3*((myj > 0) & (myj < wg.ny))); % full overlap
    Zyy = acell*sparse( (myi==nyi).*(ho+fo) );

    % clean up memory
    clear mxi nxi mxj nxj myi nyi myj nyj ho fo

    % compose the entire matrix block for this pair of layers
    Zxy = sparse(size(Zxx,1), size(Zyy,2));
    Zyx = sparse(size(Zyy,1), size(Zxx,2));
    Zxv = sparse(length(layer.xi), length(layer.vi));
    Zyv = sparse(length(layer.yi), length(layer.vi));
    Zvx = sparse(length(layer.vi), length(layer.xi));
    Zvy = sparse(length(layer.vi), length(layer.yi));
    Zvv = sparse(length(layer.vi), length(layer.vi)); % vias are ideal conductors
    Zl = [ Zxx Zxy Zxv ; Zyx Zyy Zyv ; Zvx Zvy Zvv ];

    if isfield(layer, 'conductivity')
	% If the conductivity is defined
	conductivity = layer.conductivity;
	% compute the surface impedance (as the impedance of a half
	% space) for this layer
	Zs = (1+j)*sqrt(mu0*wg.freq/(2*conductivity));
    else
	Zs = 0; % perfect conductor otherwise
    end

    % Finally scale it with this layer's surfece impedance
    Zls = Zs*Zl;

    % clean up memory
    clear Zl Zxx Zxy Zxv Zyx Zyy Zyv Zvx Zvy Zvv

    % And, finally, put this block into the overall matrix
    Z(cumbf(lidx)+1:cumbf(lidx+1), cumbf(lidx)+1:cumbf(lidx+1)) = Zls;

end
