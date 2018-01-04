function test_merge2layers
%

% X-directed two-cells wide line
B1 = zeros(5,5);
B1(:,2) = 1;
B1(:,3) = 1;

% X-directed two-cells wide line with offset
B2 = zeros(5,5);
B2(:,3) = 1;
B2(:,4) = 1;

layer1 = mklayer(B1, 0*B1, 3, 1e5);

layer2 = mklayer(B2, 0*B2, 3, 1e5);

merged_layer = merge2layers(layer1, layer2);

xixjt = [ 0 0 ; ...
          0 1 ; ...
          0 2 ; ...
          1 0 ; ...
          1 1 ; ... 
          1 2 ; ...
          2 0 ; ...
          2 1 ; ...
          2 2 ; ...
          3 0 ; ...
          3 1 ; ...
          3 2 ];
assertEquals( sortrows( xixjt ), sortrows( [ merged_layer.xi merged_layer.xj ] ) );

yiyjt = [ 0 1 ; ...
          0 2 ; ... 
          1 1 ; ...
          1 2 ; ...
          2 1 ; ...
          2 2 ];
assertEquals( sortrows( yiyjt ), sortrows( [ merged_layer.yi merged_layer.yj ] ) );

assertEquals( 3, merged_layer.pos );
assertEquals( 1e5, merged_layer.conductivity );

% Full fill
B1 = ones(10, 10);

% A line
B2 = zeros(10, 10);
B2(:,3) = 1;
B2(:,4) = 1;

layer1 = mklayer(B1, B1, 2, 2e5);

layer2 = mklayer(B2, B2, 2, 2e5);

merged_layer = merge2layers(layer1, layer2);


assertEquals( sortrows( [ layer1.xi layer1.xj ] ), ...
              sortrows( [ merged_layer.xi merged_layer.xj ] ) );

assertEquals( sortrows( [ layer1.yi layer1.yj ] ), ...
              sortrows( [ merged_layer.yi merged_layer.yj ] ) );

assertEquals( sortrows( [ layer1.vi layer1.vj ] ), ...
              sortrows( [ merged_layer.vi merged_layer.vj ] ) );

assertEquals( 2, merged_layer.pos );
assertEquals( 2e5, merged_layer.conductivity );

% Try to merge layer with itself. merge2layers follows the same basis function
% ordering as mklayer, so the layer merged with itself should be identical to the
% source one
layer2 = mklayer(B2, B2, 2, 2e5);
merged2 = merge2layers( layer2, layer2 );
assertTrue( isequal( layer2, merged2 ) );
