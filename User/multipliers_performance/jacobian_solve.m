function jacobian_solve(J,v)
  global lds
  
  J = [J;v'];
  

  ntst   = lds.ntst;
  ncol   = lds.ncol;
  nphase = lds.nphase;

  J_blocks = zeros(nphase*ncol,nphase*(ncol+1),ntst);
  j_row_indices = 1:ncol*nphase;
  j_col_indices = 1:((ncol+1)*nphase);
  for i=1:ntst
    J_blocks(:,:,i) = J(j_row_indices,j_col_indices);
    j_row_indices = j_row_indices + ncol*nphase;
    j_col_indices = j_col_indices + ncol*nphase;
  end
  

  J_bottom = full(J(end-nphase-2+1:end,:));

  
  L        = zeros(nphase*ncol+nphase+2,     nphase*ncol,         ntst);
  U        = zeros(nphase*ncol         ,     nphase*ncol,         ntst);
  

  bottom_indices  = (1:nphase*ncol);
  block_cols_1    = (1:nphase*ncol);
  block_cols_tail = (nphase*ncol+1:nphase*(ncol+1));
  jac_l           = spalloc(size(J,1),size(J,2),nphase^2*ncol^2*ntst);
  jac_u           = spalloc(size(J,1),size(J,2),nphase^2*ncol^2*ntst);
  jac_l_row_indices = (1:nphase*ncol);
  jac_l_bottom_row_indices = ((size(J,1)-nphase-2+1):size(J,1));
  jac_l_col_indices = (1:nphase*ncol);
  
  jac_u_col_indices = (1:nphase*ncol);
  jac_u_row_indices = (1:nphase*ncol);
  
  for i=1:ntst
    block = J_blocks(:,:,i);
    bordered_block = [block(:,block_cols_1);J_bottom(:,bottom_indices)];
    [L(:,:,i),U(:,:,i)] = lu(bordered_block);
%    jac_l(jac_l_row_indices, jac_l_col_indices) = L(1:ncol*nphase,:);
%    jac_l(jac_l_bottom_row_indices, jac_l_col_indices) =
%    L(ncol*nphase_1:end,:);%
%    jac_u(jac_u_col_indices,jac_u_row_indices) = U;
    
%    tail_block   = J_blocks(:,block_cols_tail);
    %u_tail_block = L(:,:,i) \ tail_block;

    
  end
    
   1; 
  
  
end