function test_findbases_vias
%

% parametes of the enclosure to pass to mkzmat
cny = 9;
cnx = cny*2;
xo  = 2; % edge-to-enclosure doubled
yo  = 2;
nx  = cnx + xo;
ny  = cny + yo;

% Canvas
B0 = zeros(nx+2, ny+2);

% Cavity upper/lower wall
B1 = drawline(B0, 1.5 + xo/2, (ny + 2)/2, nx + 0.5 - xo/2, (ny + 2)/2, cny);

% Port vias
B2a = drawline(B0, 1.5 + xo/2, (ny + 2)/2, 1.6 + xo/2, (ny + 2)/2, 1);
B2 = drawline(B2a, nx + 0.5 - xo/2, (ny + 2)/2, nx + 0.6 - xo/2, (ny + 2)/2, 1);

% Make mesh from the layers
mesh.layers = struct( [] );
mesh.layers(end + 1) = mklayer(B1, 0*B2, 1, ccopper); % Just a wall
mesh.layers(end + 1) = mklayer(B1, B2, 2, ccopper); % wall and vias
mesh.layers(end + 1) = mklayer(B1, B2, 3, ccopper); % wall and vias
mesh.layers(end + 1) = mklayer(B1, B2, 4, ccopper); % wall and vias

% Identify ports
b1 = findbases(mesh, nx, ny, 0.0, 0, 0.5, 1, @(l) l == 2);
b2 = findbases(mesh, nx, ny, 0.5, 0, 1.0, 1, @(l) l == 2);
b3 = findbases(mesh, nx, ny, 0.0, 0, 0.5, 1, @(l) l == 3);
b4 = findbases(mesh, nx, ny, 0.5, 0, 1.0, 1, @(l) l == 3);
b5 = findbases(mesh, nx, ny, 0.0, 0, 0.5, 1, @(l) l == 4);
b6 = findbases(mesh, nx, ny, 0.5, 0, 1.0, 1, @(l) l == 4);

assertEquals(b1, 2 * (length(mesh.layers(1).xi) + length(mesh.layers(1).yi)) + 1);
assertEquals(b2, 2 * (length(mesh.layers(1).xi) + length(mesh.layers(1).yi)) + 2);
assertEquals(b3, 3 * (length(mesh.layers(1).xi) + length(mesh.layers(1).yi)) + 2 + 1);
assertEquals(b4, 3 * (length(mesh.layers(1).xi) + length(mesh.layers(1).yi)) + 2 + 2);
assertEquals(b5, 4 * (length(mesh.layers(1).xi) + length(mesh.layers(1).yi)) + 4 + 1);
assertEquals(b6, 4 * (length(mesh.layers(1).xi) + length(mesh.layers(1).yi)) + 4 + 2);
