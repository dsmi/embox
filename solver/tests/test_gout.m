function test_gout
% test_gout
%
% Evaluate the target integral numerically and compare the result
% 

x0=3;
dx=0.3;

for kx=[ 0 0.1 0.3 2 ]

    gx=gout(dx,kx);

    fc = @(x) fout(x, x0, dx)*cos(kx*x);

    ifc = -gx*sin(kx*x0);
    if kx ~= 0
        ifc_test = quad(fc, x0-dx/2, x0+dx/2);
    else
        ifc_test = 0;
    end
    assertEquals(ifc_test, ifc, 1e-12);

    fs = @(x) fout(x, x0, dx)*sin(kx*x);

    ifs = gx*cos(kx*x0);
    ifs_test = quad(fs, x0-dx/2, x0+dx/2);
    assertEquals(ifs_test, ifs, 1e-12);
end
