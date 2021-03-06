function pd_vector = NP_SS_find_pd_vector(x)
  [phases_0,period,parameters] = ...
    NP_SS_extract_phases_period_and_parameters(x);

  global cds contopts
  
  % if no options have been specified yet, set default options
  if isempty(contopts)
    contopts = contset;
  end
    
  if cds.using_cvode
    cds.phases_0 = phases_0;
  else
    integration_opt = odeset(...
            'AbsTol',      contopts.integration_abs_tol / 10,    ...
            'RelTol',      contopts.integration_rel_tol / 10     ...
    );
    if ~ isempty(cds.jacobian_ode)
      integration_opt = odeset(integration_opt, ...
            'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}));
    end
    cds.cycle_orbit = cds.integrator(...
            @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
            [0 period], ...
            phases_0, integration_opt);
  end
  
  monodromy_map = ...
          @(x) NP_SS_monodromy_map(x, period, parameters);
  

  [eigenvectors, eigenvalues_matrix] = ...
          eigs(monodromy_map, cds.n_phases, cds.preferred_basis_size);
                                                    
  eigenvalues                             = diag(eigenvalues_matrix);
  distance_to_minus_one                   = abs(eigenvalues + 1);
  [~, index_of_eigenvalue_closest_to_one] = min(distance_to_minus_one);
                                                    
  if distance_to_minus_one(index_of_eigenvalue_closest_to_one) > 0.1
    warning(['The specified point does not appear to be ' ... 
            'a period doubling point of cycles, since the distance to ' ...
            'minus one of the eigenvalue closest to minus one ' ...
            'of the monodromy matrix is %.2f. ' ... 
            'At a period doubling point the monodromy matrix ' ...
            'has an eigenvalue equal to minus one.'], ...
            distance_to_minus_one(index_of_eigenvalue_closest_to_one));
    fprintf('\nPress a key to continue or ctrl-c to abort');
    pause
  end
  
  pd_vector = eigenvectors(:, index_of_eigenvalue_closest_to_one);
  
 
end
  