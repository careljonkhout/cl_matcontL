% Linear solver that replaces the Matlab backslash operator, for matrices used
% in the Newton corrections in continuation of cycles using orhtogonal
% collocation. This specialized solver will, for large enough matrices, save a
% lot of time and memory.

% The solver reduces the linear problem to a smaller "condensed" version, which
% is then solved using the Matlab backslash operator. After that, the values of
% the variables that were eliminated, are computed.

% The methods in this function are similar to the condensation of parameters
% method used in AUTO

% written by Carel Jonkhout, note: this method is NOT documented in my master
% thesis
function full_solution = linear_solver_collocation(J,b)
  global lds contopts

  ntst   = lds.ntst;
  ncol   = lds.ncol;
  nphase = lds.nphase;

  J_blocks      = zeros(nphase*ncol, nphase*(ncol+1),         ntst  );
  J_last_cols   = zeros(nphase*ncol, 2              ,         ntst  );
  J_last_rows   = zeros(2          , nphase*(ncol+1),         ntst+1);

  cs_size          = (ntst+1)*nphase+2;
  condensed_system = zeros(cs_size, cs_size);
  condensed_b      = zeros((ntst+1)* nphase+2, 1);
  transformed_b    = b;
  
  
  
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


     
  p_row_indices    = (1:nphase) + (ncol-1)*nphase;

  cs_col_indices_A = (1:nphase);
  cs_col_indices_B = (1:nphase) + nphase;
  cs_row_indices   = (1:nphase);
  
  
  b_row_indices    = (1:(ncol* nphase));
  
  testing = false || false; % added "|| false" to prevent code analyzer warnings
  
  if testing
    J_original = J;
  end
  
  condensed_b(end-nphase-1:end,:)                 = b(end-nphase-1:end,:);
   
  bottom_rect = J(end-1:end, end-1:end);
  
   j_row_indices = 1 :  ncol    * nphase;
  j_col_indices = 1 : (ncol+1) * nphase;
  
  for i=1:ntst
    sJ = J_blocks(:, (1:(ncol*nphase)) + nphase, i);
    [lower,upper] = lu(sJ);
    p = inv(lower);
    
    first_block_column = p * J_blocks(:,1:nphase,i); %#ok<MINV>
    
    % The inverse is already available, that is why the MINV warning is ignored
    % sl \  J_blocks(:,j_col_indices_A,i), may be more accurate, but it is
    % slower
    upper  = [first_block_column, upper]; %#ok<AGROW> warning is a false positive
    transformed_b(b_row_indices) = p * b(b_row_indices); %#ok<MINV>
    if testing
      J(j_row_indices,j_col_indices) = upper; 
      J(j_row_indices,end-1:end)   = p * J_last_cols(:,:,i); %#ok<MINV>
      
      disp(sum(J\transformed_b))
      disp(sum(J_original\b))
    end
    transformed_last_cols        = p * J_last_cols(:,:,i); %#ok<MINV>
    
    
    %J(j_row_indices,:)           = sl * J(j_row_indices,:);
    
    
    %disp(det(J))
    %disp(det(J_original))
    
    p = p(p_row_indices,:);
    
    

    condensed_system(cs_row_indices, cs_col_indices_A) = p * J_blocks(:,  1:nphase               , i);
    condensed_system(cs_row_indices, cs_col_indices_B) = p * J_blocks(:, (1:nphase) + ncol*nphase, i);
	  condensed_system(cs_row_indices, end-1:end)        = p * J_last_cols(:,:,i);
    condensed_b(cs_row_indices,:)               = p * b(b_row_indices,:); % todo: transformed_b can be used
    last_rows = J_last_rows(:,:,i);
    for j=1:nphase*ncol
      for k = 1:2
        if last_rows(k,nphase + j) ~= 0
          scaling_factor       = last_rows(k,nphase + j) / upper(j,nphase + j);
          last_rows(k,:)       = last_rows(k,:)       - upper(j,:)                          * scaling_factor;
          condensed_b(end-2+k) = condensed_b(end-2+k) - transformed_b(b_row_indices(1)-1+j) * scaling_factor;
          bottom_rect(k,:)     = bottom_rect(k,:)     - transformed_last_cols(j,:)          * scaling_factor;
          if testing
            transformed_b(end-2+k) = transformed_b(end-2+k) - transformed_b(b_row_indices(1)-1+j) * scaling_factor;
            J(end-2+k,:) = J(end-2+k,:) - scaling_factor * J(j_row_indices(1)-1+j,:);
          end
        end
      end
    end
    condensed_system(end-1:end, cs_col_indices_A) = ...
                             last_rows(:, 1:nphase); %#ok<*SPRIX>
   
    J_last_rows(:,1:nphase,i+1) = last_rows(:, (1:nphase) + nphase * ncol);
    
    cs_col_indices_A = cs_col_indices_A + nphase;
    cs_col_indices_B = cs_col_indices_B + nphase;
    cs_row_indices   = cs_row_indices   + nphase;
    b_row_indices    = b_row_indices    + ncol * nphase;
    
     
    j_row_indices = j_row_indices + ncol * nphase;
    j_col_indices = j_col_indices + ncol * nphase;
    
  end
 
  
 
  
  condensed_system(end-1:end, end-nphase-1:end-2) = ...
                                      J_last_rows(:, 1:nphase, ntst+1);
	condensed_system(end-nphase-1:end-2, 1:nphase          ) =   eye(nphase);
  condensed_system(end-nphase-1:end-2, end-nphase-1:end-2) = - eye(nphase);
  condensed_system(end-1:end,end-1:end) = bottom_rect;
  
  
  if testing
    rows = [];
    cols = [];

    for i=1:ntst
      cols = [cols ...
        (1:nphase) + (i-1) * nphase * ncol]; %#ok<AGROW> 
        % ok, since it is only used for testing
      rows = [rows ...
        (ncol-1) * nphase + (1:nphase) + (i-1) * nphase * ncol]; %#ok<AGROW> 
        % ok, since it is only used for testing
    end

    cols = [cols ntst * ncol * nphase + (1:nphase+2) ];
    rows = [rows ncol * ntst * nphase + (1:nphase+2) ];


    alternative_condensed_system = J(rows, cols);
    alternative_condensed_b      = transformed_b(rows);

    assert(all(abs(alternative_condensed_b - condensed_b) < 1e-14));

    assert(all(all(abs(alternative_condensed_system - condensed_system)...
          < 1e-14)));
  end
 
%  spy(condensed_system)
  if ~ isempty(contopts.lsqminnorm_threshold)
    cond_sol = lsqminnorm(sparse(condensed_system),condensed_b, ...
                                                 contopts.lsqminnorm_threshold);
  else
    cond_sol = sparse(condensed_system)\ condensed_b;
  end

  
  full_solution = zeros((ntst * ncol + 1) * nphase + 2, 1);
  full_solution_indices      = (1:nphase);
  indices_middle             = (1:(ncol-1) * nphase) + nphase;
  indices_b                  = (1:(ncol-1) * nphase);
  
  J_block_row_indices        = (1:((ncol-1) * nphase));
  J_block_col_indices_A      = 1:nphase;
  J_block_col_indices_mid    =   nphase + 1 : ncol*nphase;
  J_block_col_indices_B      =                ncol*nphase + 1 : (ncol+1)*nphase;
  
  
  cs_indices_A = (1:nphase);
  cs_indices_B = (1:nphase) + nphase;
  
  if testing
    j_row_indices   = 1 :  ((ncol-1) * nphase);
    j_col_indices_A      = 1:nphase;
    j_col_indices_mid    =   nphase + 1 : ncol*nphase;
    j_col_indices_B      =                ncol*nphase + 1 : (ncol+1)*nphase;
  end
  
  for i = 1:ntst
    full_solution(full_solution_indices) = ...
      cond_sol(cs_indices_A);
    
   
    rhs = J_blocks(J_block_row_indices, J_block_col_indices_mid, i);
    lhs = b(indices_b) - ...
      J_blocks(J_block_row_indices, J_block_col_indices_A, i) * ...
      cond_sol(cs_indices_A) - ...
      J_blocks(J_block_row_indices, J_block_col_indices_B, i) * ...
      cond_sol(cs_indices_B) - ...
      J_last_cols(J_block_row_indices,:,i) * cond_sol(end-1:end);
    
    
    if testing
      alt_rhs = J_original(j_row_indices, j_col_indices_mid);
      alt_lhs = b(indices_b) - ...
        J_original(j_row_indices, j_col_indices_A) * cond_sol(cs_indices_A) - ...
        J_original(j_row_indices, j_col_indices_B) * cond_sol(cs_indices_B) - ...
        J_original(j_row_indices, end-1:end)       * cond_sol(end-1:end);

      assert(all(all(abs(rhs-alt_rhs) < 1e-14)))
      assert(all(all(abs(lhs-alt_lhs) < 1e-14)))
      j_row_indices           = j_row_indices     + ncol * nphase;
      j_col_indices_A         = j_col_indices_A   + ncol * nphase;
      j_col_indices_mid       = j_col_indices_mid + ncol * nphase;
      j_col_indices_B         = j_col_indices_B   + ncol * nphase;
    end
   
    if ~ isempty(contopts.lsqminnorm_threshold)
      full_solution(indices_middle) = lsqminnorm(rhs, lhs, ...
                                                 contopts.lsqminnorm_threshold);
    else
      full_solution(indices_middle) = rhs\lhs;
    end
    
    %full_solution(indices_middle) = lsqminnorm(rhs,lhs,1e-3);
    full_solution(indices_middle) = rhs\lhs;
      
   % assert(all(abs(full_solution(indices_middle)-alt_sol) < 1e-14));
    
    full_solution_indices = full_solution_indices + ncol * nphase;
    indices_middle        = indices_middle        + ncol * nphase;
    indices_b             = indices_b             + ncol * nphase;
  
    cs_indices_A     = cs_indices_A      + nphase;
    cs_indices_B     = cs_indices_B      + nphase;
    
   


  end
  
  full_solution(end-nphase-1:end) = cond_sol(end-nphase-1:end);
  
  
  
end