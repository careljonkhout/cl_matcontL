function [ix, nodex] = block2d4(x1, y1, x2, y2, nx, ny)

nix = nx*ny;
nnodex = (nx+1)*(ny+1);
mex_id_ = 'block2d4(i double, i double, i double, i double, i int, i int, o int[xx], o double[xx])';
[ix, nodex] = femex(mex_id_, x1, y1, x2, y2, nx, ny, 4, nix, 2, nnodex);


