function j=ejacp(x,p)
global cds contopts
if isfield(cds, 'JacobianP') && ~isempty(cds.JacobianP)
    j = feval(cds.JacobianP, 0, x, p{:});
    j = j(:,cds.ActiveParams);
else
    opts = contopts;
    if isempty(opts)
      opts = contset();
    end
    Incr = opts.Increment;
    j = zeros(cds.ncoo, cds.nap);
    for ii = 1:length(cds.ActiveParams)
        ind = cds.ActiveParams(ii);
        p1 = p; p1{ind} = p1{ind}-Incr;
        p2 = p; p2{ind} = p2{ind}+Incr;
        j(:,ii) = feval(cds.func, 0, x, p2{:})-feval(cds.func, 0, x, p1{:});
    end
    j = j/(2*Incr);
end
