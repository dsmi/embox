function test_myfft
% test_myfft
%
% Test myfft by evaluating the coresponding sums directly.
% 

M=8;
N=8;
f=rand(M,N)+j*rand(M,N);

F=f*0;
cct=f*0;
sst=f*0;
cst=f*0;
sct=f*0;
for k=1:M
    for l=1:N
	for m=1:M
	    for n=1:N
		% This is fft2
		%F(k,l)=F(k,l)+f(m,n)*exp(-2*pi*j*(k-1)*(m-1)/M)*exp(-2*pi*j*(l-1)*(n-1)/N);
		cct(k,l)=cct(k,l)+f(m,n)*cos(2*pi*(k-1)*(m-1)/M)*cos(2*pi*(l-1)*(n-1)/N);
		sst(k,l)=sst(k,l)+f(m,n)*sin(2*pi*(k-1)*(m-1)/M)*sin(2*pi*(l-1)*(n-1)/N);
		cst(k,l)=cst(k,l)+f(m,n)*cos(2*pi*(k-1)*(m-1)/M)*sin(2*pi*(l-1)*(n-1)/N);
		sct(k,l)=sct(k,l)+f(m,n)*sin(2*pi*(k-1)*(m-1)/M)*cos(2*pi*(l-1)*(n-1)/N);
	    end
	end
    end
end

[ cc, ss, cs, sc ] = myfft(f);

assertEquals(cct, cc, 1e-12);
assertEquals(sst, ss, 1e-12);
assertEquals(cst, cs, 1e-13);
assertEquals(sct, sc, 1e-13);
