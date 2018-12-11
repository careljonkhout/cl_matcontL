function DirDer = hjacDirDerSquare(x, p, Q1, vec)

% Compute hessian matrices of F(x) numerically
%
% contCISjacDirDer is a matrix:
%   hess(i,j,k) = d^2 F_i / dx_j dx_k ????

global contopts

Incr = contopts.Increment;

Qvec = Q1 * vec;
x1 = x;
x1 = x1 - (Incr/norm(Qvec)) * Qvec;
J1 = hjac(x1,p);
x2 = x;
x2 = x2 + (Incr/norm(Qvec)) * Qvec;
J2 = hjac(x2,p);
% DirDer = Q1'* (J2 - J1) * ( norm(Qvec)/(2*Incr) );
DirDer = (J2 - J1) * ( norm(Qvec)/(2*Incr) ); % DV