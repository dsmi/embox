function [ ne, nm ] = wnorm(a, b, maxm, maxn)
% function [ ne, nm ] = wnorm(a, b, maxm, maxn)
%
% Retruns the normalization coefficients for the modes of the rectangular
% waveguide. The modes are normalized such that integral of magnitude of
% the mode vector dotted with itself over the guide crossection is unity.
%  Params:
%    a, b       - dimensions of the waveguide
%    maxm, maxn - upper limits of the mode indices
%  Outputs:
%    ne, nm     - normalization coefficients for the te and tm modes
%                 correspondingly
%

% wm(m,n)=m-1, wn(m,n)=n-1
[ wm, wn ] = ndgrid(0:maxm-1, 0:maxn-1);

% Normalization coefficient for TE modes
ne=1/pi*sqrt(2*2*a*b./((wm.*b).^2+(wn.*a).^2));
ne(1,:)=1/pi*sqrt(2*a*b./((wm(1,:).*b).^2+(wn(1,:).*a).^2)); % 2 in numerator
ne(:,1)=1/pi*sqrt(2*a*b./((wm(:,1).*b).^2+(wn(:,1).*a).^2)); % 2 in numerator
ne(1,1)=0; % no TE00 mode

% Normalization coefficient for TM modes
nm=1/pi*sqrt(2*2*a*b./((wm.*b).^2+(wn.*a).^2));
nm(1,1)=0; % no TM00 mode
