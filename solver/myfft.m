function [ cc ss cs sc ] = myfft(f)
% [ cc ss cs sc ] = myfft(f)
%
% Computes the following sums
%  cc(k,l) = sum(sum(f(m,n)*cos(2*pi*(k-1)*(m-1)/M)*cos(2*pi*(l-1)*(n-1)/N)))
%  ss(k,l) = sum(sum(f(m,n)*sin(2*pi*(k-1)*(m-1)/M)*sin(2*pi*(l-1)*(n-1)/N)))
%  cs(k,l) = sum(sum(f(m,n)*cos(2*pi*(k-1)*(m-1)/M)*sin(2*pi*(l-1)*(n-1)/N)))
%  sc(k,l) = sum(sum(f(m,n)*sin(2*pi*(k-1)*(m-1)/M)*cos(2*pi*(l-1)*(n-1)/N)))
% where the summation ranges are m=1..M n=1..N using the FFT.
%
% Here is an explanation of how this works
% 2d fft computes the following sum:
%  F(k,l)=sum(sum(f(m,n)*exp(-2*pi*j*(k-1)*(m-1)/M)*exp(-2*pi*j*(l-1)*(n-1)/N)))
%       =sum(sum(f(m,n)*exp(-2*pi*j*((k-1)*(m-1)/M+(l-1)*(n-1)/N))))
%
% Flipping the fft matrix along first or second dimension (around the zero
% frequency point) is equivalent to changing the sign of the exponent of the
% first or second terms (ones with k-m and l-n correspondingly)
%
% The sums we are interested in can be represented as:
%   cos(a)*cos(b)=(exp(-j*a)+exp(j*a))*(exp(-j*b)+exp(j*b))/4=
%            =(exp(j*(-a-b))+exp(j*(-a+b))+exp(j*(a-b))+exp(j*(a+b)))/4
%   sin(a)*sin(b)=-(exp(-j*a)-exp(j*a))*(exp(-j*b)-exp(j*b))/4=
%            =(-exp(j*(-a-b))+exp(j*(a-b))+exp(j*(-a+b))-exp(j*(a+b)))/4
%   cos(a)*sin(b)=(exp(-j*a)+exp(j*a))*(exp(-j*b)-exp(j*b))/(-j*4)=
%            =(exp(j*(-a-b))-exp(j*(-a+b))+exp(j*(a-b))-exp(j*(a+b)))/(-j*4)
%   sin(a)*cos(b)=(exp(-j*a)-exp(j*a))*(exp(-j*b)+exp(j*b))/(-j*4)=
%            =(exp(j*(-a-b))+exp(j*(-a+b))-exp(j*(a-b))-exp(j*(a+b)))/(-j*4)
%
% The rest is obvious.
%

F=fft2(f);

F10=fftshift(F);
F01=fftshift(F);
F11=fftshift(F);

if mod(size(F, 1), 2)
    F10=flip(F10,1); % if M is odd
    F11=flip(F11,1);
else 
    F10(2:end,:)=flip(F10(2:end,:),1); % if even
    F11(2:end,:)=flip(F11(2:end,:),1);
end 

if mod(size(F, 2), 2)
    F01=flip(F01,2); % if N is odd
    F11=flip(F11,2);
else 
    F01(:,2:end)=flip(F01(:,2:end),2); % if even
    F11(:,2:end)=flip(F11(:,2:end),2);
end 

F10=ifftshift(F10);
F01=ifftshift(F01);
F11=ifftshift(F11);

cc=(F+F10+F01+F11)./4;
ss=(-F+F10+F01-F11)./4;
cs=(F+F10-F01-F11)./(-j*4);
sc=(F-F10+F01-F11)./(-j*4);
