function test_gtri
% test_gtri
%
% Evaluate the target integral numerically and compare the result
% 

x0=3;
dx=0.3;

for kx=[ 0 0.1 0.3 2 ]
    gx=gtri(dx,kx);

    fc = @(x) ftri(x, x0, dx)*cos(kx*x);

    ifc = gx*cos(kx*x0);
    ifc_test = quad(fc, x0-dx, x0+dx, 1e-20);
    assertEquals(ifc_test, ifc, 1e-12);

    fs = @(x) ftri(x, x0, dx)*sin(kx*x);

    ifs = gx*sin(kx*x0);
    ifs_test = quad(fs, x0-dx, x0+dx);
    assertEquals(ifs_test, ifs, 1e-12);
end
