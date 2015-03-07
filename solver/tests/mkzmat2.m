function Z=mkzmat2(wg, mesh)
% Z=mkzmat(wg, mesh)
% Populates the impedance/reactions matrix in a straightforward and slow
% manner, used for tests. Only works correctly if:
%  The waveguide has two layers
%  The mesh has one layer of metallization
%
% wg     - shiedling parameters, see wgparams
%

% mesh cell sizes
dx=wg.a/wg.nx;
dy=wg.b/wg.ny;

% this is supposed to be the only layer in the mesh
layer = mesh.layers(1);

nx=length(layer.xi);
ny=length(layer.yi);
n=nx+ny;

Z=zeros(n,n);

% Zxx
for m=1:nx,
    for n=1:nx,
	xim = layer.xi(m);
	xjm = layer.xj(m);
	posm = layer.pos;
	xin = layer.xi(n);
	xjn = layer.xj(n);
	posn = layer.pos;
	Z(m,n)=calczmn(wg, xim, xjm, posm, 1, xin, xjn, posn, 1);
    end
end

% Zxy
for m=1:nx,
    for n=1:ny,
	xim = layer.xi(m);
	xjm = layer.xj(m);
	posm = layer.pos;
	yin = layer.yi(n);
	yjn = layer.yj(n);
	posn = layer.pos;
	Z(m,nx+n)=calczmn(wg, xim, xjm, posm, 1, yin, yjn, posn, 0);
    end
end

% Zyx
for m=1:ny,
    for n=1:nx,
	yim = layer.yi(m);
	yjm = layer.yj(m);
	posm = layer.pos;
	xin = layer.xi(n);
	xjn = layer.xj(n);
	posn = layer.pos;
	Z(nx+m,n)=calczmn(wg, yim, yjm, posm, 0, xin, xjn, posn, 1);
    end
end

% Zyy
for m=1:ny,
    for n=1:ny,
	yim = layer.yi(m);
	yjm = layer.yj(m);
	posm = layer.pos;
	yin = layer.yi(n);
	yjn = layer.yj(n);
	posn = layer.pos;
	Z(nx+m,nx+n)=calczmn(wg, yim, yjm, posm, 0, yin, yjn, posn, 0);
    end
end

% Identify segments which cross the waveguide boundary - the corresponding
% elements of the Z matrix need to be multiplied by 0.5
xmul=-(~layer.xi | ~(layer.xi-wg.nx))*0.5+1.0;
ymul=-(~layer.yj | ~(layer.yj-wg.ny))*0.5+1.0;
M=diag([ xmul(:) ; ymul(:) ]);

% scale the impedance matrix elements for the boundary-crossing segments
Z=M*Z*M;
