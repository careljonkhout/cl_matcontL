function vec2 = multilinear2(q1,q2,x0,p,increment)
%----------------------------------------------------------
%This file computes the multilinear function B(q1,q2) where
%B = D^2(F(x0)), the second derivative of the vectorfield
%----------------------------------------------------
global cds contopts

if contopts.SymDerivative(2) && isfield(cds, 'Hessian') && ~isempty(cds.Hessian)
    Hess = feval(cds.Hessian, 0, x0, p{:});
    vec2 = tensor2op(Hess, q1, q2, cds.ncoo);
else
    if (q1==q2)
        vec2 = Bvv(q1,x0,p,increment);
    else
        part1 = Bvv(q1+q2,x0,p,increment);
        part2 = Bvv(q1-q2,x0,p,increment);
        vec2 = (part1-part2)/4.0;
    end
end


%----------------------------------------------------
function tempvec = Bvv(vq,x0,p,increment)

global cds

f0 = x0;
f1 = x0 + increment*(vq);
f2 = x0 - increment*(vq);
f0 = feval(cds.func, 0, f0, p{:});
f1 = feval(cds.func, 0, f1, p{:});
f2 = feval(cds.func, 0, f2, p{:});
tempvec = (f1+f2-2.0*f0)/(increment^2);
