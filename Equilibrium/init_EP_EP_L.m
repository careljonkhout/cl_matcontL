function [x0, v0] = init_EP_EP_L(probfile, x0, p, ap)
%print_diag(5,'In init_EP_EP_L\n');

  global cds
  load_odefile(probfile)

  if isempty(x0)
    if ~isempty(init_prob)
      p2 = num2cell(p);
      x0 = feval(init_prob, p2{:});
    else
      error('No initial point found');
    end
  end
  x0 = [x0; p(ap)];
  v0 = [];

  % initialize cds
  cds.P0           = p;
  cds.ActiveParams = ap;
  cds.ndim         = length(x0);
  cds.ncoo         = length(x0) - 1;
  cds.nap          = length(ap);
end