function test_findbases
%

% X-directed two-cells wide line
B=zeros(8,8);

B(:,4)=1;
B(:,5)=1;

mesh=mkmesh(B);

b = findbases(mesh, 6, 6, 0, 0, 0, 1);
assertEquals(0, mesh.xi(b));

b = findbases(mesh, 6, 6, 1, 0, 1, 1);
assertEquals(6, mesh.xi(b));

% Y-directed two-cells wide line
B=zeros(8,8);

B(4,:)=1;
B(5,:)=1;

mesh=mkmesh(B);

b = findbases(mesh, 6, 6, 0, 0, 1, 0);
assertEquals(0, mesh.yj(b-length(mesh.xi)));

b = findbases(mesh, 6, 6, 0, 1, 1, 1);
assertEquals(6, mesh.yj(b-length(mesh.xi)));
