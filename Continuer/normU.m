function val = normU(arg)
% Allows for calculation of user defined norm for discretized PDE Problems
% etc...  Called from newtcorrL

global cds

if isempty(cds.usernorm)
    val = norm(arg);
else
    u = arg(1:cds.ncoo);
    p = arg(cds.ncoo+1:end);
    val = sqrt(feval(cds.usernorm,u)^2 + norm(p)^2);
end