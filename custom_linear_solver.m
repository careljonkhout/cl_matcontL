function full_solution = custom_linear_solver(J,b)
  global lds

  ntst   = lds.ntst;
  ncol   = lds.ncol;
  nphase = lds.nphase;

  J_blocks      = zeros(nphase*ncol, nphase*(ncol+1),         ntst  );
  J_last_cols   = zeros(nphase*ncol, 2              ,         ntst  );
  J_last_rows   = zeros(2          , nphase*(ncol+1),         ntst+1);

  cs_size          = (ntst+1)*nphase+2;
  nnz              = ntst * nphase^2 + 4 * ( ntst + 1 ) * nphase  - 4;
  condensed_system = spalloc(cs_size, cs_size, nnz);
  condensed_b      = zeros((ntst+1)* nphase+2, 1);
  

  % Retrieve blocks of Jacobian first and store them as dense matrices to prevent
  % repeating expensive lookups of blocks in sparse matririces.
  j_row_indices = 1 :  ncol    * nphase;
  j_col_indices = 1 : (ncol+1) * nphase;
  for i=1:ntst
    J_blocks(:,:,i)    = J(j_row_indices, j_col_indices);
    J_last_cols(:,:,i) = J(j_row_indices, end-1:end);
    J_last_rows(:,:,i) = J(end-1:end    , j_col_indices);
    j_row_indices = j_row_indices + ncol * nphase;
    j_col_indices = j_col_indices + ncol * nphase;
  end

  
  last_rows_col_indices_A = (1:nphase);
  last_rows_col_indices_B = (1:nphase) + nphase *  ncol;
    
  p_row_indices    = (1:nphase) + (ncol-1)*nphase;
  j_col_indices_p  = (1:(ncol*nphase)) + nphase;
  j_col_indices_A  = (1:nphase);
  j_col_indices_B  = (1:nphase) + ncol*nphase;
  cs_col_indices_A = (1:nphase);
  cs_col_indices_B = (1:nphase) + nphase;
  cs_row_indices   = (1:nphase);
  
  b_row_indices    = (1:(ncol* nphase));
  
  for i=1:ntst
    sJ = J_blocks(:,j_col_indices_p,i);
    [sl,su] = lu(sJ);
    p = inv(sl);
    
    first_block_column = p * J_blocks(:,j_col_indices_A,i);
    su = [first_block_column, su];
    
    p = p(p_row_indices,:);
    
    

    condensed_system(cs_row_indices,cs_col_indices_A) = ...
            p * J_blocks(:,j_col_indices_A,i);
    condensed_system(cs_row_indices,cs_col_indices_B) = ...
            p * J_blocks(:,j_col_indices_B,i);
	  condensed_system(cs_row_indices, end-1:end) = p * J_last_cols(:,:,i);
    condensed_b(cs_row_indices,:) = p * b(b_row_indices,:);
    last_rows = J_last_rows(:,:,i);
    for j=1:nphase*(ncol-1)
      for k = 1:2
        if last_rows(k,nphase + j) ~= 0
          scaling_factor = last_rows(k,nphase + j) / su(j,nphase + j);
          last_rows(k,:) = last_rows(k,:) - su(j,:) * scaling_factor;
        end
      end
    end
    condensed_system(end-1:end, cs_col_indices_A) = ...
                              last_rows(:, last_rows_col_indices_A);
   
    J_last_rows(:,1:nphase,i+1) = last_rows(:, last_rows_col_indices_B);
    
    cs_col_indices_A = cs_col_indices_A + nphase;
    cs_col_indices_B = cs_col_indices_B + nphase;
    cs_row_indices   = cs_row_indices   + nphase;
    b_row_indices    = b_row_indices    + ncol * nphase;
  end
 
  
  condensed_system(end-1:end, end-nphase-1:end-2) = ...
                                      J_last_rows(:, 1:nphase, ntst+1);
	condensed_system(end-nphase-1:end-2, 1:nphase          ) =   eye(nphase);
  condensed_system(end-nphase-1:end-2, end-nphase-1:end-2) = - eye(nphase);
  condensed_b(end-nphase-2:end,:)                 = b(end-nphase-2,:);
%  spy(condensed_system)
  condensed_solution = condensed_system \ condensed_b;
  
  full_solution = zeros((ntst * ncol + 1) * nphase + 2, 1);
  full_solution_indices_A      = (1:nphase);
  indices_middle               = (1:(ncol-1) * nphase) + nphase;
  
  J_block_row_indices     = (1:(ncol-1) * nphase);
  J_block_col_indices_A   = 1:nphase;
  J_block_col_indices_mid = (1:(ncol-1) * nphase) + nphase;
  J_block_col_indices_B   = (1:nphase)  + ncol * nphase;   
  
  
  cond_sol_indices_A = (1:nphase);
  cond_sol_indices_B = (1:nphase) + nphase;
  
  
  for i = 1:ntst
    full_solution(full_solution_indices_A) = ...
      condensed_solution(cond_sol_indices_A);
    
   
    middle_rhs = J_blocks(J_block_row_indices, J_block_col_indices_mid);
    middle_lhs = b(indices_middle,1) - ...
      J_blocks(J_block_row_indices, J_block_col_indices_A) * ...
      condensed_solution(cond_sol_indices_A,1) - ...
      J_blocks(J_block_row_indices, J_block_col_indices_B) * ...
      condensed_solution(cond_sol_indices_B,1);
      
    full_solution(indices_middle) = middle_rhs \ middle_lhs;
      
    
    full_solution_indices_A = full_solution_indices_A + ncol * nphase;
    indices_middle          = indices_middle          + ncol * nphase;
  
    cond_sol_indices_A      = cond_sol_indices_A      + nphase;
    cond_sol_indices_B      = cond_sol_indices_B      + nphase;

  end
  
  full_solution(end-nphase-1:end) = condensed_solution(end-nphase-1:end);
  
  
  
end