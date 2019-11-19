function pd_vector = find_pd_vector(x)
   [phases_0,period,parameters] = ...
    NewtonPicard.MultipleShooting.extract_phases_period_and_parameters(x);

  global cds contopts
  
  % if no options have been specified yet, set default options
  if isempty(contopts)
    contopts = contset;
  end
    
  if cds.using_cvode
    cds.phases_0 = phases_0;
  else
    NewtonPicard.MultipleShooting.compute_stitched_orbit()
  end
  
  monodromy_map = ...
   @(x) NewtonPicard.MultipleShooting.full_monodromy_map(x, period, parameters);
  

  [eigenvectors, eigenvalues_matrix] = ...
          eigs(monodromy_map, cds.n_phases, cds.preferred_basis_size);
                                                    
  eigenvalues                                   = diag(eigenvalues_matrix);
  distance_to_minus_one                         = abs(eigenvalues + 1);
  [~, index_of_eigenvalue_closest_to_minus_one] = min(distance_to_minus_one);
                                                    
  if distance_to_minus_one(index_of_eigenvalue_closest_to_minus_one) > 0.1
    warning(['The specified point does not appear to be ' ... 
            'a period doubling point of cycles, since the distance to one ' ...
            'of the eigenvalue closest to one ' ...
            'of the monodormy matrix is %.2f. ' ... 
            'At a period doubling point the monodromy matrix ' ...
            'has an eigenvalue equal to minus one'], ...
            distance_to_minus_one(index_of_eigenvalue_closest_to_one));
    fprintf('\nPress a key to continue or ctrl-c to abort');
    pause
  end
  
  pd_vector = eigenvectors(:, index_of_eigenvalue_closest_to_minus_one);
  
end
  