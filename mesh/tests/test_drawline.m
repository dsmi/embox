function test_drawline
%

C=eye(50);
B0=C*0;
B=drawline(B0,0,0,50,50,0.5);

assertEquals(C, B);
