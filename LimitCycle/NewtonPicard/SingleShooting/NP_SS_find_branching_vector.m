function v = NP_SS_find_branching_vector(x, v_orig)
   [V, reduced_jacobian, delta_q_gamma, delta_q_r, ~, ~, ~, ~, ~] = ...
          NP_SS_compute_reduced_jacobian(x);

 
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
  
  % compute a 2 dimensional approximation of the nullspace of the Jacobian by
  % taking the 2 right singular vectors associated to the 2 smallest singular
  % values.
  
  delta_p_s      = delta_p__delta_T_and_delta_gamma(1:end-2, :); 
  delta_T_s      = delta_p__delta_T_and_delta_gamma(end-1  , :);
  delta_gamma_s  = delta_p__delta_T_and_delta_gamma(end    , :);
  delta_q_s      = delta_q_r + delta_q_gamma * delta_gamma_s;
  delta_phases   = V * delta_p_s + delta_q_s;
  nullspace      = [delta_phases; delta_T_s; delta_gamma_s];
  inner_prod     = abs(nullspace' * v_orig);
  
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
  