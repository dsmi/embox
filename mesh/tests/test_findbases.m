function test_findbases
%

% X-directed two-cells wide line
B=zeros(8,8);

B(:,4)=1;
B(:,5)=1;

layer = mklayer(B);
layer.pos = 2;

mesh.layers(1) = layer;

b = findbases(mesh, 6, 6, 0, 0, 0, 1);
assertEquals(0, mesh.layers(1).xi(b));

b = findbases(mesh, 6, 6, 1, 0, 1, 1);
assertEquals(6, mesh.layers(1).xi(b));

% Y-directed two-cells wide line
B=zeros(8,8);

B(4,:)=1;
B(5,:)=1;

layer=mklayer(B);
layer.pos = 2;

mesh.layers(1) = layer;

b = findbases(mesh, 6, 6, 0, 0, 1, 0);
assertEquals(0, mesh.layers(1).yj(b-length(mesh.layers(1).xi)));

b = findbases(mesh, 6, 6, 0, 1, 1, 1);
assertEquals(6, mesh.layers(1).yj(b-length(mesh.layers(1).xi)));
