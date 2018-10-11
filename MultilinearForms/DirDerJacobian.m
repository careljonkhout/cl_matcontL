function DirDer = DirDerJacobian(x, vec, JacMatrix)

% (Numerically) computes hessian matrices of F(X) acting on the vector vec
% with n-1 components, i.e.
%
% DirDer = d/dz ( F_X(x + z*vec) )
%
% if JacMatrix is specified, it is assumed to be F_X(x) without the last
% column

global cds contopts

ndim = cds.ndim;
Incr = contopts.Increment;

if nargin == 3
    J1 = JacMatrix;
else
    J1 = contjac(x);
    J1 = J1(:,1:end-1);
end

x2 = x;
x2(1:ndim-1) = x2(1:ndim-1) + (Incr/norm(vec)) * vec;
J2 = contjac(x2);
J2 = J2(:,1:end-1);

DirDer = (J2 - J1) * norm(vec) / Incr;