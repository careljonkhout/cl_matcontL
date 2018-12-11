function j = ejac(x,p)
global cds contopts ejac_prev

%Prevents recalculation if last point is the same as this point
if isfield(ejac_prev, 'point_x') && ~isempty(ejac_prev.point_x)
    if isequal(x, ejac_prev.point_x) && isequal(p, ejac_prev.point_p)
        j = ejac_prev.jac;
        return;
    end
end
   
% compute jacobian
if contopts.SymDerivative(1) && isfield(cds, 'Jacobian') && ~isempty(cds.Jacobian)
    j = feval(cds.Jacobian, 0, x, p{:});
else
    Incr = contopts.Increment;
    func = cds.func;
    j = zeros(cds.ncoo);
    if contopts.contL_ParallelComputing
        parfor i=1:cds.ncoo
            x1 = x; x1(i) = x1(i)-Incr;
            x2 = x; x2(i) = x2(i)+Incr;
            j(:,i) = feval(func, 0, x2, p{:})-feval(func, 0, x1, p{:});
        end
    else
        for i=1:cds.ncoo
            x1 = x; x1(i) = x1(i)-Incr;
            x2 = x; x2(i) = x2(i)+Incr;
            j(:,i) = feval(func, 0, x2, p{:})-feval(func, 0, x1, p{:});
        end
    end
    j = j/(2*Incr);
    j = sparse(j);
end

% store to prevent recalculation
ejac_prev.jac = j;
ejac_prev.point_x = x;
ejac_prev.point_p = p;