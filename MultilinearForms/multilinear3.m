function vec3 = multilinear3(q1,q2,q3,x0,p,increment)
%--------------------------------------------------------------
% This file computes the multilinear function C(q1,q2,q3) where
% C = D^3(F(x0)), the 3rd derivative of the vectorfield wrt to phase
% variables only.
%--------------------------------------------------------------
global cds contopts

if contopts.SymDerivative(3) && isfield(cds, 'der3') && ~isempty(cds.der3)
    der3 = feval(cds.der3, 0, x0, p{:});  %MP
    vec3 = tensor3op(der3, q1, q2, q3, cds.ncoo);
else
    if (q1==q2)
        if (q1==q3)
            vec3 = Cvvv(q1,x0,p,increment);
        else
            part1 = Cvvv(q1+q3,x0,p,increment);
            part2 = Cvvv(q1-q3,x0,p,increment);
            part3 = Cvvv(q3,x0,p,increment);
            vec3 = (part1 - part2 - 2.0*part3)/6.0;
        end
    else
        part1 = Cvvv(q1+q2+q3,x0,p,increment);
        part2 = Cvvv(q1+q2-q3,x0,p,increment);
        part3 = Cvvv(q1-q2+q3,x0,p,increment);
        part4 = Cvvv(q1-q2-q3,x0,p,increment);
        vec3 = (part1 - part2 - part3 + part4)/24.0;
    end
end

%----------------------------------------------------
function tempvec = Cvvv(vq,x0,p,increment)

global cds
f1 = x0 + 3.0*increment*vq;
f2 = x0 +     increment*vq;
f3 = x0 -     increment*vq;
f4 = x0 - 3.0*increment*vq;
f1 = feval(cds.func, 0, f1, p{:});
f2 = feval(cds.func, 0, f2, p{:});
f3 = feval(cds.func, 0, f3, p{:});
f4 = feval(cds.func, 0, f4, p{:});
tempvec = (f1 - 3.0*f2 + 3.0*f3 - f4)/(8.0*increment^3);
