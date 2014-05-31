function test_gtri
% test_gtri
%
% Evaluate the target integral numerically and compare the result
% 

x0=3;
dx=0.3;

for kx=[ 0 0.1 0.3 2 ]
    gx=gtri(dx,kx);

    f1 = @(x) (x-x0)/dx+1;
    f2 = @(x) (x0-x)/dx+1;

    f1c = @(x) f1(x)*cos(kx*x);
    f2c = @(x) f2(x)*cos(kx*x);

    ifc = gx*cos(kx*x0);
    ifc_test = quad(f1c, x0-dx, x0)+quad(f2c, x0, x0+dx);
    assertEquals(ifc_test, ifc, 1e-12);

    f1s = @(x) f1(x)*sin(kx*x);
    f2s = @(x) f2(x)*sin(kx*x);

    ifs = gx*sin(kx*x0);
    ifs_test = quad(f1s, x0-dx, x0)+quad(f2s, x0, x0+dx);
    assertEquals(ifs_test, ifs, 1e-12);
end
