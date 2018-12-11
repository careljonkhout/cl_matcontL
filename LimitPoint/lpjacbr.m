function j = lpjacbr(x,p)
global cds contopts

if contopts.SymDerivativeP(1) >= 1
  j = feval(cds.JacobianP, 0, x, p{:});
  j = j(:,cds.BranchParams);
else
    Incr = contopts.Cont_IncrFinDiff;
  for i = cds.BranchParams
    p1 = p; p1{i} = p1{i}-Incr;
    p2 = p; p2{i} = p2{i}+Incr;
    j(:,i) = feval(cds.func, 0, x, p2{:})-feval(cds.func, 0, x, p1{:});
  end
  j = j(:,cds.BranchParams)/(2*Incr);

end
