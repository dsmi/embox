function test_mergelayers
%

% Full fill
B1 = ones(10,10);

% Line
B2 = zeros(10,10);
B2(:,6) = 1;
B2(:,7) = 1;

% Another Line
B3 = zeros(10,10);
B3(5,:) = 1;
B3(6,:) = 1;

mesh.layers = struct( [] );
mesh.layers(end+1) = mklayer(B1, B1, 10, 1e5);
mesh.layers(end+1) = mklayer(B2, B2, 20, 1e5);
mesh.layers(end+1) = mklayer(B1, B1, 30, 1e5);
mesh.layers(end+1) = mklayer(B3, B2, 20, 1e5); % duplicate 20
mesh.layers(end+1) = mklayer(B1, B1, 40, 1e5);
mesh.layers(end+1) = mklayer(B2, B2, 10, 1e5); % duplicate 10
mesh.layers(end+1) = mklayer(B1, B1, 50, 1e5);

merged_mesh = mergelayers( mesh );

assertTrue( isequal( merge2layers( mesh.layers(1), mesh.layers(6) ), merged_mesh.layers(1) ) )
assertTrue( isequal( merge2layers( mesh.layers(2), mesh.layers(4) ), merged_mesh.layers(2) ) )
assertTrue( isequal( mesh.layers(3), merged_mesh.layers(3) ) )
assertTrue( isequal( mesh.layers(5), merged_mesh.layers(4) ) )
assertTrue( isequal( mesh.layers(7), merged_mesh.layers(5) ) )
