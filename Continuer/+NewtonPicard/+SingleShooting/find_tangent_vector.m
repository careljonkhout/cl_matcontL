function v = find_tangent_vector(x)

   [V, reduced_jacobian, delta_q_gamma, delta_q_r, ~, ~, ~, ~, ~] = ...
          NewtonPicard.SingleShooting.compute_reduced_jacobian(x);

	if (rank(reduced_jacobian) < min(size(reduced_jacobian)))
    error(['The reduced Jacobian is rank deficient. ' ...
           'Cannot find tangent vector. '])
  end
        
  delta_p__delta_T_and_delta_gamma = null(reduced_jacobian);
  delta_p     = delta_p__delta_T_and_delta_gamma(1:end-2);
  delta_T     = delta_p__delta_T_and_delta_gamma(end-1);
  delta_gamma = delta_p__delta_T_and_delta_gamma(end);

  
  delta_q        = delta_q_r + delta_gamma * delta_q_gamma;
  delta_phases   = V * delta_p + delta_q;
  v              = [delta_phases; delta_T; delta_gamma];
    
