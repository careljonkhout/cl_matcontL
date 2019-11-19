function v = find_tangent_vector(x)
  global cds
  m = cds.n_mesh_intervals;
  
  [V, reduced_jacobian, delta_q_gamma, ~, ~, ~, ~, ~, ~] = ...
    NewtonPicard.MultipleShooting.compute_reduced_jacobian(x);

%   [~, ~, delta_p__delta_T_and_delta_gamma] = ...
%           svds(reduced_jacobian, 1, 'smallest');
  delta_p__delta_T_and_delta_gamma = null(full(reduced_jacobian));
  delta_p__delta_T_and_delta_gamma = delta_p__delta_T_and_delta_gamma(:,1);
        
  delta_p     = delta_p__delta_T_and_delta_gamma(1:end-2);
  delta_T     = delta_p__delta_T_and_delta_gamma(end-1);
  delta_gamma = delta_p__delta_T_and_delta_gamma(end);

  V_delta_p = zeros(cds.n_phases * m, 1);
  col_offset = 0;
  for i = 1 : m
    indices_delta_p              = col_offset + ( 1 : size(V{i}, 2) );
    indices_V_delta_p            = (i-1) * cds.n_phases + (1:cds.n_phases);
    V_delta_p(indices_V_delta_p) = V{i} * delta_p(indices_delta_p);
    col_offset                   = col_offset + size(V{i}, 2);
  end
  
  delta_q      = delta_gamma * delta_q_gamma;
  delta_phases = V_delta_p + delta_q(:); % x = x(:) turns x into a column vector
  v = [delta_phases; delta_T; delta_gamma];
end
    
