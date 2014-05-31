function Z=mkzmat2(wg, mesh)
% Z=mkzmat(wg, mesh)
% Populates the impedance/reactions matrix
%
% wg     - shiedling parameters, see wgparams
%

% mesh cell sizes
dx=wg.a/wg.nx;
dy=wg.b/wg.ny;

nx=length(mesh.xi);
ny=length(mesh.yi);
n=nx+ny;

Z=zeros(n,n);

% Zxx
for m=1:nx,
    for n=1:nx,
	Z(m,n)=calczmn(wg, mesh.xi(m), mesh.xj(m), 1, mesh.xi(n), mesh.xj(n), 1);
    end
end

% Zxy
for m=1:nx,
    for n=1:ny,
	Z(m,nx+n)=calczmn(wg, mesh.xi(m), mesh.xj(m), 1, mesh.yi(n), mesh.yj(n), 0);
    end
end

% Zyx
for m=1:ny,
    for n=1:nx,
	Z(nx+m,n)=calczmn(wg, mesh.yi(m), mesh.yj(m), 0, mesh.xi(n), mesh.xj(n), 1);
    end
end

% Zyy
for m=1:ny,
    for n=1:ny,
	Z(nx+m,nx+n)=calczmn(wg, mesh.yi(m), mesh.yj(m), 0, mesh.yi(n), mesh.yj(n), 0);
    end
end

% Identify segments which cross the waveguide boundary - the corresponding
% elements of the Z matrix need to be multiplied by 0.5
xmul=-(~mesh.xi | ~(mesh.xi-wg.nx))*0.5+1.0;
ymul=-(~mesh.yj | ~(mesh.yj-wg.ny))*0.5+1.0;
M=diag([ xmul(:) ; ymul(:) ]);

% scale the impedance matrix elements for the boundary-crossing segments
Z=M*Z*M;
