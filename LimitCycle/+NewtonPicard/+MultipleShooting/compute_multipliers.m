function multipliers = compute_multipliers(x, nMults_to_compute)
  global cds contopts;
  print_diag(3,'computing multipliers\n');
  
  [cds.phases_0, period, parameters] = ...
          NewtonPicard.MultipleShooting.extract_phases_period_and_parameters(x);
  
  if ( ~ contopts.monodromy_by_finite_differences) && ( ~ cds.using_cvode )
    NewtonPicard.MultipleShooting.compute_cycle_parts(x);
  end
  
  M                      = @(x) ...
     NewtonPicard.MultipleShooting.full_monodromy_map(1, x, period, parameters);
  nMults_to_compute      = min(nMults_to_compute, cds.n_phases);
  [~, multiplier_matrix] = eigs(M, cds.n_phases, nMults_to_compute);                                    
  multipliers            = diag(multiplier_matrix);
  
  print_diag(1, multipliers2str(multipliers));
end