function [ Y I Z ]=solvey(wg, mesh, ports, portw)
% [ Y I Z ]=solvey(wg, mesh, ports, portw)
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
%    I      - mesh_size-by-num_of_ports matrix, I(:,n) gives the currents due
%             to the unit voltage applied to the port n
%    Z      - the generalized impedance matrix (for all the basis functions)
%

% The impedance matrix assuming perfect conductivity
Zpc = mkzmat(wg, mesh);

% surface impedance
Zs = mkzsmat(wg, mesh);

% Impedance matrix with surface impedace
Z = Zpc - Zs;

% The linear system solution will need memory
clear Zpc Zs

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

% Currents due to the voltages applied to the ports
I = Z\V;

% And the admittance matrix
Y = -V.'*I;
