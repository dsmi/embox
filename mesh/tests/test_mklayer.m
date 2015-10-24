function test_mklayer
%

% X-directed two-cells wide line
B=zeros(8,8);

B(:,4)=1;
B(:,5)=1;

layer=mklayer(B);

xixj=layer.xi*10+layer.xj;
xixjt = [ 2 ; 12 ; 22 ; 32 ; 42 ; 52 ; 62 ; 3 ; 13 ; 23 ; 33 ; 43 ; 53 ; 63 ];
assertEquals(sort(xixjt), sort(xixj));

yiyj=layer.yi*10+layer.yj;
yiyjt = [ 3 ; 13 ; 23 ; 33 ; 43 ; 53 ];
assertEquals(sort(yiyjt), sort(yiyj));

% Y-directed two-cells wide line
B=zeros(8,8);

B(4,:)=1;
B(5,:)=1;

layer=mklayer(B);

xixj=layer.xi*10+layer.xj;
xixjt = [ 30 ; 31 ; 32 ; 33 ; 34 ; 35 ];
assertEquals(sort(xixjt), sort(xixj));

yiyj=layer.yi*10+layer.yj;
yiyjt = [ 20 ; 30 ; 21 ; 31 ; 22 ; 32 ; 23 ; 33 ; 24 ; 34 ; 25 ; 35 ; 26 ; 36 ];
assertEquals(sort(yiyjt), sort(yiyj));

% Fully-populated layer
% Y-directed two-cells wide line
B=ones(6,6);

layer=mklayer(B);
assertEquals(20, length(layer.xi));
assertEquals(20, length(layer.xj));
assertEquals(20, length(layer.yi));
assertEquals(20, length(layer.yj));

% one via in the left-bottom (0, 0) corner
B0=zeros(8,8);
B=zeros(8,8);

B(2,2)=1;

layer=mklayer(B0, B);

assertEquals(0, layer.vi);
assertEquals(0, layer.vj);

% some vias
B0=zeros(8,8);
B=zeros(8,8);

B(2:3,2:4)=1;

layer=mklayer(B0, B);

vivj=layer.vi*10+layer.vj;
vivjt = [ 0 ; 1 ; 2 ; 10 ; 11 ; 12 ];
assertEquals(sort(vivjt), sort(vivj));
