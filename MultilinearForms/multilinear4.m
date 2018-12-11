function vec4 = multilinear4(q1,q2,q3,q4,x0,p,increment)
%--------------------------------------------------------------
%This file computes the multilinear function D(q1,q2,q3,q4) where
%D = D^4(F(x0)), the fourth derivative of the map wrt to phase
%variables only. First we decide whether q1=q2, then the rest.
%--------------------------------------------------------------
global cds contopts

if contopts.SymDerivative(4) && isfield(cds, 'der4') && ~isempty(cds.der4)
    der4 = feval(cds.Hessian, 0, x0, p{:});
    vec4 = tensor4op(der4, q1, q2, q3, q4, cds.ncoo);
else
    if isequal(q1,q2)
        if isequal(q1,q3)
            if isequal(q1,q4)
                vec4 = Dvvvv(q1,x0,p,increment);
            else
                part1 = Dvvvv(3.0*q1+q4,x0,p,increment);
                part2 = Dvvvv(3.0*q1-q4,x0,p,increment);
                part3 = Dvvvv(q1+q4,x0,p,increment);
                part4 = Dvvvv(q1-q4,x0,p,increment);
                vec4 = (part1 - part2 - 3.0*part3 + 3.0*part4)/192.0;
            end
        elseif (q3==q4)
            part1 = Dvvvv(q1+q3,x0,p,increment);
            part2 = Dvvvv(q1-q3,x0,p,increment);
            part3 = Dvvvv(q1   ,x0,p,increment);
            part4 = Dvvvv(   q3,x0,p,increment);
            vec4 = (part1 + part2 - 2.0*part3 - 2.0*part4)/12.0;
        else
            part1 = Dvvvv(2.0*q1+q3+q4,x0,p,increment);
            part2 = Dvvvv(2.0*q1+q3-q4,x0,p,increment);
            part3 = Dvvvv(2.0*q1-q3+q4,x0,p,increment);
            part4 = Dvvvv(2.0*q1-q3-q4,x0,p,increment);
            part5 = Dvvvv(       q3+q4,x0,p,increment);
            part6 = Dvvvv(      -q3+q4,x0,p,increment);
            vec4 = (part1 - part2 - part3 + part4 - 2.0*part5 + 2.0*part6)/192.0;
        end
    else
        part1 = Dvvvv(q1+q2+q3+q4,x0,p,increment);
        part2 = Dvvvv(q1+q2+q3-q4,x0,p,increment);
        part3 = Dvvvv(q1+q2-q3+q4,x0,p,increment);
        part4 = Dvvvv(q1+q2-q3-q4,x0,p,increment);
        part5 = Dvvvv(q1-q2+q3+q4,x0,p,increment);
        part6 = Dvvvv(q1-q2+q3-q4,x0,p,increment);
        part7 = Dvvvv(q1-q2-q3+q4,x0,p,increment);
        part8 = Dvvvv(q1-q2-q3-q4,x0,p,increment);
        vec4 = (part1 - part2 - part3 + part4 - part5 + part6 + part7 - part8)/192.0;
    end
end

%----------------------------------------------------
function tempvec = Dvvvv(vq,x0,p,increment)

global cds
f0 = x0;
f1 = x0 + 4.0*increment*vq;
f2 = x0 + 2.0*increment*vq;
f3 = x0 - 2.0*increment*vq;
f4 = x0 - 4.0*increment*vq;
f0 = feval(cds.func, 0, f0, p{:});
f1 = feval(cds.func, 0, f1, p{:});
f2 = feval(cds.func, 0, f2, p{:});
f3 = feval(cds.func, 0, f3, p{:});
f4 = feval(cds.func, 0, f4, p{:});
tempvec = (f1 - 4.0*f2 + 6.0*f0 - 4.0*f3 + f4)/(16.0*increment^4);
