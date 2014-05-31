function test_gflat
% test_gflat
%
% Evaluate the target integral numerically and compare the result
% 

x0=3;
dx=0.3;

for kx=[ 0 0.1 0.3 2 ]
    gx=gflat(dx,kx);

    fc = @(x) cos(kx*x);

    ifc = gx*cos(kx*x0);
    ifc_test = quad(fc, x0-dx/2, x0+dx/2);
    assertEquals(ifc_test, ifc, 1e-12);

    fs = @(x) sin(kx*x);

    ifs = gx*sin(kx*x0);
    ifs_test = quad(fs, x0-dx/2, x0+dx/2);
    assertEquals(ifs_test, ifs, 1e-12);
end
