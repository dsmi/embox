
freq = linspace(0, 1e14, 500);

e1 = debye(4.3, 0.02, 1e9*2*pi, freq*2*pi);
e2 = debye2(4.3, 0.02, 1e9*2*pi, freq*2*pi);

plot(freq, imag(e1), '-r', freq, imag(e2), '-b')
% plot(freq, real(e1), '-*r', freq, real(e2), '-*b')
