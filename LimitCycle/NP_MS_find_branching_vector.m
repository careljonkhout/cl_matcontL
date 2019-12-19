function v = NP_MS_find_branching_vector(x, v_orig)
  global cds
  m = cds.n_mesh_intervals;
  
  [V, reduced_jacobian, delta_q_gamma, delta_q_r, ~, ~, ~, ~, ~] = ...
    NP_MS_compute_reduced_jacobian(x);
  
  basis_size = size(V,2);

  % take the right singular vectors corresponding to the smallest two singular
  % values
  [~, s, delta_p__delta_T_and_delta_gamma] = ...
                                          svds(reduced_jacobian, 2, 'smallest');
  
  if max(abs(diag(s))) > 0.1
    warning(['The specified point does not appear to be ' ... 
          'a branching point of cycles, since the second smallest singular ' ...
          'value of the Jacobian is %.2f. At a branching point of cycles, ' ...
          'the Jacobian has two singular values equal to zero.'], max(diag(s)));
    fprintf('\nPress a key to continue or ctrl-c to abort\n');
    pause
  end
    
  % append _s to variable names to indicate that there are two delta_p's, two
  % delta_T's, etc..
  delta_p_s     = delta_p__delta_T_and_delta_gamma(1:end-2, :);
  delta_T_s     = delta_p__delta_T_and_delta_gamma(end-1  , :);
  delta_gamma_s = delta_p__delta_T_and_delta_gamma(end    , :);

  
  % We project the parts of the tangent vector that correspond the P_i -
  % components of the tangent vector at each mesh point to the full dimenional
  % space. Note that the P_i space is spanned the columns of V(:,:,i).
  V_delta_p_s = zeros(cds.n_phases * m, 2);
  for i=1:m % m == cds.n_mesh_intervals
    indices1 = (i-1) * cds.n_phases + (1:cds.n_phases);
    indices2 = (i-1) * basis_size  + (1:basis_size );
    V_delta_p_s(indices1, :) = V(:,:,i) * delta_p_s(indices2, :);
  end
  
  % change delta_q_r, and delta_q_gamma into column vectors
  delta_q_r     = delta_q_r(:);
  delta_q_gamma = delta_q_gamma(:); 
 
  % add the q_r and q_gamma components of the Q components of the part of the
  % null-vectors that corresponds to the tangent vectors at the mesh points.
  delta_q_s     = delta_q_r + delta_q_gamma * delta_gamma_s;
  
  % add the P and Q components of the part of the null-vectors that correspond
  % to the tangent vector at the mesh points
  delta_phases  = V_delta_p_s + delta_q_s;
  
  % assemble the null vectors 
  nullspace     = [delta_phases; delta_T_s; delta_gamma_s];
  
  inner_prod    = abs(nullspace' * v_orig);
  
  % we return the vector of the approximated nullspace which has the smallest
  % inner product with the original tangent vector v_orig.
  [~, i] = min(inner_prod);
  v = nullspace(:,i);
  
  % if the angle between v and v_orig is greater than 90 degrees, reverse the
  % direction of v
  if v' * v_orig < 0
    v = - v;
  end
end
    
