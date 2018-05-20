
addpath(genpath([ pwd, '/..' ]))

mil2meter = 2.54e-5;
ha = pi/4; % hatch angle
ho = 0;    % hatch offset
hp = 16/sin(ha)*mil2meter;  % hatch pitch
hw = 3.5/sin(ha)*mil2meter;   % hatch width
l = 12*16*mil2meter


N = 48*4;
B0 = zeros(N+2, N+2);

BH = drawhatch(B0, hp/l, hw/l, ha, ho/l);
imagesc(BH)
lay = mklayer(BH, 0*BH);

% number of basis functions and Z matrix memory
nbf = numel(lay.xi) + numel(lay.yi)
memZgb = ((nbf^2)*16)/(1024^3)

% total memory with line mesh
lbf = 12320
tbf = lbf + nbf
totZgb = ((tbf^2)*16)/(1024^3)




