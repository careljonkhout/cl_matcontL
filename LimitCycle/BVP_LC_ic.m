% -------------------------------------------------------------
% functions defining the BVP for limitcycles
% -------------------------------------------------------------

function ic = BVP_LC_ic(ups)
global lds
% At each fine grid point compute dot product of point coordinates and the
% derivative at the corresponding point in the cycle of the previous
% continuation step. We have:
% 
% size(lds.upoldp) == size(ups) == [lds.nphase lds.ntst*lds.ncol+1]
% == [lds.nphase lds.tps]
% size(ficd) == [1  lds.tps]
ficd = dot(ups,lds.upoldp);
% In the next step the values of ficd associated to each mesh interval are
% collected in the columns of ficdmat
% the matrix lds.idxmat looks like this:
% [ 1        ncol+1       ....  tps - ncol    ;
%   2        ncol+2       ....  tps - ncol + 1;
%   ...                         ....
%   ncol+1   2*ncol+1     ....  tps]
% and we have size(lds.idxmat) == size(ficdmat) == [ncol+1 ntst]
ficdmat = ficd(lds.idxmat);
% finally the integral is computed. We have size(lds.wi) = [ncol+1 1]. lds.wi
% contains the Gaussian quadrature coeffici\"ents.
ic = sum(lds.dt.*(lds.wi*ficdmat));
