function writeTouch(fileName, freq, y)
% writeTouch(fileName, freq, y)
%
% Writes multiport network parameters in Touchstone format
% Now only handles 2-port y-parameters
%

comment = 'no comment supplied';

f = fopen(fileName, 'wt');

fprintf(f, '%s\n', ['! ' comment]);
fprintf(f, '%s\n', '# Hz Y RI R 50');

N = size(y, 3);
for i = 1:N
    fprintf(f, '%.11e %.11e %.11e %.11e %.11e %.11e %.11e %.11e %.11e\n',  ...
            freq(i), ...
            real(y(1,1,i)), imag(y(1,1,i)), real(y(2,1,i)), imag(y(2,1,i)), ...
            real(y(1,2,i)), imag(y(1,2,i)), real(y(2,2,i)), imag(y(2,2,i)) );
end

fclose(f);
