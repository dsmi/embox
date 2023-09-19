
addpath(genpath([ pwd, '/..' ]));

[ fhz, Svf ] = SXPParse( 'via_96_via.y2p', stdout );

Yvf = s2y( renorms( Svf, [ 50 50 ], [ 1 1 ] ) );
freqs = 2*pi*fhz;

Z1 = Z2 = Z3 = 0*freqs;
Lp = Cp1 = Cp2 = 0*freqs;

for fidx = 1:length(freqs)

    freq = freqs(fidx);
    fhz = freq/(2*pi);
    Y = Yvf(:,:,fidx);
    Z = inv(Y);

    % T-network
    Z1(fidx) = Z(1,1) - Z(2,1); % half-series
    Z2(fidx) = Z(2,1);          % shunt
    Z3(fidx) = Z(2,2) - Z(2,1); % half-series 2


    % pi-network
    Y1(fidx) = Y(1,1) + Y(1,2); % half-shunt
    Y2(fidx) = -Y(1,2);         % series
    Y3(fidx) = Y(2,2) + Y(1,2); % half-shunt 2
    
    
end

plot( freqs/(2*pi), real(Z1), '-*r', freqs/(2*pi), imag(Z1), '-*b' )
legend( 'real', 'imag' )
xlabel('freq, Hz')
ylabel('Z_1, Ohm')
title('Series admittance')

%% plot( freqs/(2*pi), 0*real(Z2), '-*r', freqs/(2*pi), imag(1./Z2), '-*b' )
%% legend( 'real', 'imag' )
%% xlabel('freq, Hz')
%% ylabel('1/Z_3, S')
%% title('Shunt conductance')

%% Ct = j*imag(1./Y2)./(j*freqs)
%% Lt1 = j*imag(Z1)/(j*freq);
%% Lt2 = j*imag(Z3)/(j*freq);

%% Lp = j*imag(1/Y2)/(j*freq);
%% Cp1 = j*imag(Y1)/(j*freq);
%% Cp2 = j*imag(Y3)/(j*freq);
