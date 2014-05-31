function wg = wgparams(freq,a,b,h,c,nx,ny)
% wg = wgparams(freq,a,b,h,c,nx,ny)
%
% Returns the structure which describes the enclosure (which is treated as
% the waveguide when calculating the fields) parameters initialized by default.
%
% The strucure contains the following fields:
%   freq      Angular frequency
%   a, b      x and y dimensions of the waveguide
%   h         height of the metallization above ground
%   c         height of the upper ground
%   nx, ny    number of cells along x and y axes correspondingly.
%   cnx, cny  number of divisions of each cell (along x and y) used when
%             summing the waveguide modes. nx*cnx-1 and ny*cny-2 give the
%             indices of the higher mode processed when evaluating the reaction
%             integrals. The fft2 dimensions are nx*cnx*2-by-ny*cny*2, from
%             the performace considerations it is better when the dimensions
%             are powers of two. Both cnx and cny must be divisible by 4.
%   weps, wmu Layers stackup, from bottom to top along z
%   Gls0      left-looking reflection coefficient at the bottom (z=0)
%             Gls0=-1 for perfect conductor, Gls0=0 for matched
%   Ggr0      right-looking reflection coefficient at the top (z=h)
%             Ggr0=-1 for perfect conductor, Ggr0=0 for matched
%

wg.freq = freq;
wg.a    = a;
wg.b    = b;
wg.h    = h;
wg.c    = c;
if exist('nx')
    wg.nx = nx;
else
    wg.nx = 64;
end
if exist('ny')
    wg.ny   = ny;
else
    wg.ny   = 64;
end
wg.cnx  = 16;
wg.cny  = 16;
wg.weps = [ eps0 eps0 ];
wg.wmu  = [ mu0 mu0 ];
wg.Gls0 = -1;
wg.Ggr0 = -1;
