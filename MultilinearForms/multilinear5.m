function vec5 = multilinear5(q1,q2,q3,q4,q5,x0,p,increment)
%--------------------------------------------------------------
%This file computes the multilinear function E(q1,q2,q3,q4,q5) where
%E = D^5(F(x0)), the fifth derivative of the vectorfield wrt to phase
%variables only. We use this for normal form computations in which we
%will have q1=q2=q3=q4=q5 or q1=q2=q3\neq q4=q5. We decide on these
%cases only. Otherwise we just compute the thing directly without
%optimization.
%--------------------------------------------------------------
global cds contopts

if contopts.SymDerivative(5) && isfield(cds, 'der5') && ~isempty(cds.der5)
    der5 = feval(cds.der5, 0, x0, p{:});
    vec5 = tensor5op(der5, q1, q2, q3, q4, q5, cds.ncoo);
else
    if (isequal(q1,q2) && isequal(q1,q3))
        if (isequal(q1,q4) && isequal(q1,q5))
            vec5 = Evvvvv(q1,x0,p,increment);
        else
            part1 = Evvvvv(3.0*q1+2.0*q4,x0,p,increment);
            part2 = Evvvvv(3.0*q1-2.0*q4,x0,p,increment);
            part3 = Evvvvv(3.0*q1       ,x0,p,increment);
            part4 = Evvvvv(    q1       ,x0,p,increment);
            part5 = Evvvvv(    q1+2.0*q4,x0,p,increment);
            part6 = Evvvvv(    q1-2.0*q4,x0,p,increment);
            vec5 = (part1 + part2 - 2.0*part3 + 6.0*part4 - 3.0*part5 - 3.0*part6)/1920.0;
        end
    else
        part1 = Evvvvv(q1+q2+q3+q4+q5,x0,p,increment);
        part2 = Evvvvv(q1+q2+q3+q4-q5,x0,p,increment);
        part3 = Evvvvv(q1+q2+q3-q4-q5,x0,p,increment);
        part4 = Evvvvv(q1+q2+q3-q4+q5,x0,p,increment);
        part5 = Evvvvv(q1+q2-q3+q4+q5,x0,p,increment);
        part6 = Evvvvv(q1+q2-q3+q4-q5,x0,p,increment);
        part7 = Evvvvv(q1+q2-q3-q4-q5,x0,p,increment);
        part8 = Evvvvv(q1+q2-q3-q4+q5,x0,p,increment);
        vec5 = (part1 - part2 + part3 - part4 - part5 + part6 - part7 + part8)/1920.0;
        part1 = Evvvvv(q1-q2+q3+q4+q5,x0,p,increment);
        part2 = Evvvvv(q1-q2+q3+q4-q5,x0,p,increment);
        part3 = Evvvvv(q1-q2+q3-q4-q5,x0,p,increment);
        part4 = Evvvvv(q1-q2+q3-q4+q5,x0,p,increment);
        part5 = Evvvvv(q1-q2-q3+q4+q5,x0,p,increment);
        part6 = Evvvvv(q1-q2-q3+q4-q5,x0,p,increment);
        part7 = Evvvvv(q1-q2-q3-q4-q5,x0,p,increment);
        part8 = Evvvvv(q1-q2-q3-q4+q5,x0,p,increment);
        vec5 = vec5 - (part1 - part2 + part3 - part4 - part5 + part6 - part7 + part8)/1920.0;
    end
end

%----------------------------------------------------
function tempvec = Evvvvv(vq,x0,p,increment)

global cds
f1 = x0 + 5.0*increment*vq;
f2 = x0 + 3.0*increment*vq;
f3 = x0 + 1.0*increment*vq;
f4 = x0 - 1.0*increment*vq;
f5 = x0 - 3.0*increment*vq;
f6 = x0 - 5.0*increment*vq;
f1 = feval(cds.func, 0, f1, p{:});
f2 = feval(cds.func, 0, f2, p{:});
f3 = feval(cds.func, 0, f3, p{:});
f4 = feval(cds.func, 0, f4, p{:});
f5 = feval(cds.func, 0, f5, p{:});
f6 = feval(cds.func, 0, f6, p{:});
tempvec =  (f1 - 5.0*f2 + 10.0*f3 - 10.0*f4 + 5.0*f5 - f6)/(32.0*increment^5);
