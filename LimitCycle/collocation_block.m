function block = collocation_block(x,i)
  global cds lds

  ntst   = lds.ntst;
  ncol   = lds.ncol;
  nphase = lds.nphase;
  
  if ~ isfield(lds,'wt') || isempty(lds.wt)
    calc_weights
  end
  
  if ~ isfield(lds,'dt') || isempty(lds.dt)
    set_mesh()
  end
  
  
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = x(end);
  parameters                   = num2cell(parameters);

  period   = x(end-1);
  ups      = reshape(x(1:end-2), nphase, ntst * ncol + 1);
  ups_cols = (1:(ncol+1)) + (i-1) * ncol;
  sysjac   = zeros(nphase * ncol, nphase * (ncol+1));
  
 
  block_rows = 1:nphase;
  % value of polynomial in each mesh point
  xval  = ups(:,ups_cols) * lds.wt;
  wploc = lds.wp/lds.dt(i);
  % evaluate part of Jacobian matrix
  for j = 1:(ncol)
    xtmp = xval(:,j);
    jac                  = cjac(cds.dydt_ode,cds.jacobian_ode,xtmp, ...
                                                 parameters,cds.ActiveParams);
    sysjac(block_rows,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,j)',jac);
    block_rows           = block_rows + lds.nphase;
  end

  block = wploc - period * sysjac;

end

function calc_weights()
% calculate weights

  global lds

  lds.wi = nc_weight(lds.ncol)';
  lds.wt = [];
  lds.wpvec = [];
  zm = gl_pos(lds.ncol)/2 + 0.5;
  xm = (0:lds.ncol)/lds.ncol;

  ncp1 = lds.ncol+1;
  for j=1:ncp1
    xmi = xm(setdiff(1:ncp1,j));
    denom = prod(xm(j)-xmi);
    lds.wt(j,:) = prod(repmat(zm',lds.ncol,1)-repmat(xmi',1,lds.ncol))/denom;
    s = zeros(1,lds.ncol);
    for l=setdiff(1:ncp1,j)
      xmil = xm(setdiff(1:ncp1,[j l]));
      rep = repmat(zm',lds.ncol-1,1)-repmat(xmil',1,lds.ncol);
      if lds.ncol > 2
        s = s+prod(rep);
      else
        s = s+rep;
      end
    end
    lds.wpvec(j,:) = s/denom;
  end
  lds.wp = kron(lds.wpvec',eye(lds.nphase));
end

function set_mesh()
  global lds
  lds.msh = linspace(0,1,lds.ntst+1);
  lds.finemsh = get_fine_mesh(lds.msh, lds.ntst, lds.ncol);
  lds.dt = lds.msh(2:lds.ntst+1)-lds.msh(1:lds.ntst);
end
