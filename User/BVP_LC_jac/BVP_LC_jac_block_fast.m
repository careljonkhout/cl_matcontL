function blocks = BVP_LC_jac_block_fast(x)
  global cds lds

  ntst   = lds.ntst;
  ncol   = lds.ncol;
  nphase = lds.nphase;

  
  
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = x(end);
  parameters                   = num2cell(parameters);
 
  blocks = zeros(nphase*ncol,nphase*(ncol+1),ntst);

  period   = x(end - 1);
  ups      = reshape(x(1:end-2), nphase, ncol*ntst+1);
  ups_cols = 1:(ncol+1);
  sysjac   = zeros(nphase * ncol, nphase * (ncol+1));
  
  for i = 1:(ntst)
    block_rows = 1:lds.nphase;
    % value of polynomial in each mesh point
    xval  = ups(:,ups_cols) * lds.wt;
    wploc = lds.wp/lds.dt(i);
    % evaluate part of Jacobian matrix
    for j = 1:(ncol)
      xtmp = xval(:,j);
      jac                  = cjac(lds.func,lds.Jacobian,xtmp, ...
                                                   parameters,lds.ActiveParams);
      sysjac(block_rows,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,j)',jac);
      block_rows           = block_rows + lds.nphase;
    end
    % Store result
   
    blocks(:,:,i) = wploc - period * sysjac;
   
    % shift to following intervals
    ups_cols = ups_cols + lds.ncol;
  end
end