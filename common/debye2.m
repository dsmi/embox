function epsr = debye2(Dk0, Df0, w0, w1)
% Another debye to compare
% given Dk, Df at w0, finds Dk and Df at w1
% 

wmin = 1.0e4;
wmax = 1.0e14;
k12 = log(wmin / wmax);

enorm = @(wf) 0.5*log((wmin.^2 + wf.^2) ./ (wmax.^2 + wf.^2)) ./ k12 ...
                         + j * (atan(wf ./ wmin) - atan(wf ./ wmax)) ./ k12;

z = enorm(w0);
K = -Dk0 .* sin(Df0) ./ imag(z);
einf = Dk0 .* cos(Df0) - K .* real(z);

z = enorm(w1);
epsr = einf + K * z;

%% Dk1 = abs(epsr);
%% Df1 = -atan(imag(epsr) / real(epsr));
