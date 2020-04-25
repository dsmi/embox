function [ D, Y1, Y2 ] = deembsims( fsim1, fsim2 )
% [ D, Y1, Y2 ] = deembsims( fsim1, fsim2 )
%
%  Runs two simulations of the deembedding standards of length l and 2*l,
% by calling fsim1 and fsim2 correspondingly, and calculates the port(s)
% discontinuity characterization.
%

fprintf('Running deembedding-l simulation...\n')
Y1 = fsim1( );

fprintf('Running deembedding-l*2 simulation...\n')
Y2 = fsim2( );

A1 = y2abcd(Y1);
A2 = y2abcd(Y2);

% Obtain the port discontinuity
[ D dtol ] = splitd( A1*inv(A2)*A1 );
fprintf( 'Deembedding accuracy is %.8e\n', dtol );


