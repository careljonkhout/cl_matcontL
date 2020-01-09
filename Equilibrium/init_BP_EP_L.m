function [x0, v0]= init_BP_EP_L(odefile, x, p, ap, data)
  %
  % [x0,v0]= init_BP_EP_L(odefile, x, p, s, h)
  %
  % Branch switching at an Branching Point (BP) detected on an equilibrium
  % curve
  %
  % data from located branching point must be available
  % If x, p or ap are empty the values from the data struct are used
  %
  
  if isempty(data)
    error('Data structure must be given')
  end

  if isempty(x)
    x = data.x(1:end-1);
  end

  if isempty(p)
    p = data.P0;
  end

  if isempty(ap)
    ap = data.ap;
  end
  if length(ap) ~= 1
    error('Requires only one active parameter')
  end
  if ap ~= data.ap
    error('Active parameter must be the same as on the equilibrium curve')
  end
  x0 = [x; p(ap)];

  global cds
  cds = [];
  
  load_odefile(odefile); % adds fields to cds

  cds.P0           = p;
  cds.ActiveParams = ap;
  cds.ndim         = length(x0);
  cds.ncoo         = length(x0) - 1;
  cds.nap          = length(ap);
  cds.curve        = @equilibriumL;

  if isfield(data, 's') && isfield(data.s,'xnext')
    x0 = data.s.xnext;
    v0 = data.s.vnext;
  elseif isfield(data.s, 'v2') || ~isempty(data.s.v2)
    v0 = data.s.v2;
  else
    error('No tangent vector found in data. BP might be degenerate');
  end
end