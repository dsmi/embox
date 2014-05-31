function [ Y Z ]=solvey(wg, mesh, ports, portw)
% [ Y Z ]=solvey(wg, mesh, ports, portw)
%
% Given a system of conductors with N ports, this function calculates
% the admittance matrix Y of size N-by-N.
% Inputs:
%    wg     - shiedling parameters, see wgparams
%    mesh   - meshed conductors, see mkmesh
%    ports  - cell array of vectors, indices of the basis functions forming
%             the ports.
%    portw  - cell array of vectors of the same size as ports, widths of the
%             corresponding basis functions.
% Outputs:
%    Y      - the admittance matrix.
%    Z      - the generalized impedance matrix (for all the basis functions)
%

% The impedance matrix!
Z=mkzmat(wg, mesh);

% Number of the basis functions
N=size(Z,1);

% Number of the ports.
np = length(ports);

% Width of the basis functions forming the ports
portwi = cell2mat(portw);

% Build the voltage excitation matrix - the number of columns is the number
% of ports, each column applies unity voltage to the basis function(s)
% destignated for the corresponding port and zero voltage for all the others.
[ faceidx fportidx ] = ports2subs(ports);
V = zeros(N,np);
V(sub2ind(size(V), faceidx, fportidx)) = portwi;

Y=-V.'*(Z\V);
