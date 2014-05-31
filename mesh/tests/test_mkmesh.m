function test_mkmesh
%

% X-directed two-cells wide line
B=zeros(8,8);

B(:,4)=1;
B(:,5)=1;

mesh=mkmesh(B);

xixj=mesh.xi*10+mesh.xj;
xixjt = [ 2 ; 12 ; 22 ; 32 ; 42 ; 52 ; 62 ; 3 ; 13 ; 23 ; 33 ; 43 ; 53 ; 63 ];
assertEquals(sort(xixjt), sort(xixj));

yiyj=mesh.yi*10+mesh.yj;
yiyjt = [ 3 ; 13 ; 23 ; 33 ; 43 ; 53 ];
assertEquals(sort(yiyjt), sort(yiyj));

% Y-directed two-cells wide line
B=zeros(8,8);

B(4,:)=1;
B(5,:)=1;

mesh=mkmesh(B);

xixj=mesh.xi*10+mesh.xj;
xixjt = [ 30 ; 31 ; 32 ; 33 ; 34 ; 35 ];
assertEquals(sort(xixjt), sort(xixj));

yiyj=mesh.yi*10+mesh.yj;
yiyjt = [ 20 ; 30 ; 21 ; 31 ; 22 ; 32 ; 23 ; 33 ; 24 ; 34 ; 25 ; 35 ; 26 ; 36 ];
assertEquals(sort(yiyjt), sort(yiyj));

% Fully-populated mesh
% Y-directed two-cells wide line
B=ones(6,6);

mesh=mkmesh(B);
assertEquals(20, length(mesh.xi));
assertEquals(20, length(mesh.xj));
assertEquals(20, length(mesh.yi));
assertEquals(20, length(mesh.yj));
