function Z=mkzmat2(wg, mesh)
% Z=mkzmat2(wg, mesh)
% Populates the impedance/reactions matrix in a straightforward and slow
% manner, used for tests.
%
% wg     - shiedling parameters, see wgparams
% mesh   - meshed metal, see mkmesh
%

% mesh cell sizes
dx=wg.a/wg.nx;
dy=wg.b/wg.ny;

% We want to pre-allocate the Z matrix, for that we need to know its size.
% To calculate the size we need to count the basis functions on all layers.
numx = sum(cellfun(@(v) length(v), { mesh.layers(:).xi }));
numy = sum(cellfun(@(v) length(v), { mesh.layers(:).yi }));
numv = sum(cellfun(@(v) length(v), { mesh.layers(:).vi }));
numbf = numx + numy;

% We also compute the cumulative sums of numbers of the basis functions in
% the layers up to the given one which is then used to place the blocks
% of the matrix (see below)
cumx = cumsum(cellfun(@(v) length(v), { mesh.layers(:).xi }));
cumy = cumsum(cellfun(@(v) length(v), { mesh.layers(:).yi }));
cumv = cumsum(cellfun(@(v) length(v), { mesh.layers(:).vi }));
cumbf = [ 0 (cumx + cumy + cumv) ];

Z=zeros(numbf,numbf);

% via current scaling
viac = 1/dx;

for mli = 1:length(mesh.layers)

    mlay = mesh.layers(mli);

    % Position of the m-th layer in the stackup
    mpos = mlay.pos;

    for nli = 1:length(mesh.layers)

	nlay = mesh.layers(nli);

	% Position of the n-th layer in the stackup
	npos = nlay.pos;

	% Number of basis functions in m and n layers
	nmx = length(mlay.xi);
	nmy = length(mlay.yi);
	nmv = length(mlay.vi);
	nnx = length(nlay.xi);
	nny = length(nlay.yi);
	nnv = length(nlay.vi);

	% Zxx
	for m=1:nmx,
	    for n=1:nnx,
		xim = mlay.xi(m);
		xjm = mlay.xj(m);
		xin = nlay.xi(n);
		xjn = nlay.xj(n);
		z = calczmn(wg, xim, xjm, mpos, 1, xin, xjn, npos, 1);
		Z(cumbf(mli)+m,cumbf(nli)+n) = z;
	    end
	end

	% Zxy
	for m=1:nmx,
	    for n=1:nny,
		xim = mlay.xi(m);
		xjm = mlay.xj(m);
		yin = nlay.yi(n);
		yjn = nlay.yj(n);
		z = calczmn(wg, xim, xjm, mpos, 1, yin, yjn, npos, 0);
		Z(cumbf(mli)+m,cumbf(nli)+nnx+n) = z;
	    end
	end

	% Zxv
	for m=1:nmx,
	    for n=1:nnv,
		xim = mlay.xi(m);
		xjm = mlay.xj(m);
		vin = nlay.vi(n);
		vjn = nlay.vj(n);
		z = calczmn(wg, xim, xjm, mpos, 1, vin, vjn, npos, 2);
		Z(cumbf(mli)+m,cumbf(nli)+nnx+nny+n) = z*viac;
	    end
	end

	% Zyx
	for m=1:nmy,
	    for n=1:nnx,
		yim = mlay.yi(m);
		yjm = mlay.yj(m);
		xin = nlay.xi(n);
		xjn = nlay.xj(n);
		z = calczmn(wg, yim, yjm, mpos, 0, xin, xjn, npos, 1);
		Z(cumbf(mli)+nmx+m,cumbf(nli)+n) = z;
	    end
	end

	% Zyy
	for m=1:nmy,
	    for n=1:nny,
		yim = mlay.yi(m);
		yjm = mlay.yj(m);
		yin = nlay.yi(n);
		yjn = nlay.yj(n);
		z = calczmn(wg, yim, yjm, mpos, 0, yin, yjn, npos, 0);
		Z(cumbf(mli)+nmx+m,cumbf(nli)+nnx+n) = z;
	    end
	end
	
	% Zyv
	for m=1:nmy,
	    for n=1:nnv,
		yim = mlay.yi(m);
		yjm = mlay.yj(m);
		vin = nlay.vi(n);
		vjn = nlay.vj(n);
		z = calczmn(wg, yim, yjm, mpos, 0, vin, vjn, npos, 2);
		Z(cumbf(mli)+nmx+m,cumbf(nli)+nnx+nny+n) = z*viac;
	    end
	end
	
	% Zvx
	for m=1:nmv
	    for n=1:nnx
		vim = mlay.vi(m);
		vjm = mlay.vj(m);
		xin = nlay.xi(n);
		xjn = nlay.xj(n);
		z = calczmn(wg, vim, vjm, mpos, 2, xin, xjn, npos, 1);
		Z(cumbf(mli)+nmx+nmy+m,cumbf(nli)+n) = z*viac;
	    end
	end

	% Zvy
	for m=1:nmv
	    for n=1:nny
		vim = mlay.vi(m);
		vjm = mlay.vj(m);
		yin = nlay.yi(n);
		yjn = nlay.yj(n);
		z = calczmn(wg, vim, vjm, mpos, 2, yin, yjn, npos, 0);
		Z(cumbf(mli)+nmx+nmy+m,cumbf(nli)+nnx+n) = z*viac;
	    end
	end

	% Zvv
	for m=1:nmv
	    for n=1:nnv
		vim = mlay.vi(m);
		vjm = mlay.vj(m);
		vin = nlay.vi(n);
		vjn = nlay.vj(n);
		z = calczmn(wg, vim, vjm, mpos, 2, vin, vjn, npos, 2);
		Z(cumbf(mli)+nmx+nmy+m,cumbf(nli)+nnx+nny+n) = z*viac*viac;
	    end
	end

	% Identify segments which cross the waveguide boundary - the
	% corresponding elements of the Z matrix need to be multiplied by 0.5
	% The source segmetns which cross:
	mxmul=-(~mlay.xi | ~(mlay.xi-wg.nx))*0.5+1.0;
	mymul=-(~mlay.yj | ~(mlay.yj-wg.ny))*0.5+1.0;
	mvmul = mlay.vi*0 + 1.0; % vias can not cross the boundary
	Mm=diag([ mxmul(:) ; mymul(:) ; mvmul(:) ]);

	% And the observation segmetns which cross:
	nxmul=-(~nlay.xi | ~(nlay.xi-wg.nx))*0.5+1.0;
	nymul=-(~nlay.yj | ~(nlay.yj-wg.ny))*0.5+1.0;
	nvmul = nlay.vi*0 + 1.0; % vias can not cross the boundary
	Mn=diag([ nxmul(:) ; nymul(:) ; nvmul(:) ]);

	% scale the boundary crossing segments in the block we just populated
	Zl = Z(cumbf(mli)+1:cumbf(mli+1), cumbf(nli)+1:cumbf(nli+1));
	Zl = Mm*Zl*Mn;
	Z(cumbf(mli)+1:cumbf(mli+1), cumbf(nli)+1:cumbf(nli+1)) = Zl;

    end
end
