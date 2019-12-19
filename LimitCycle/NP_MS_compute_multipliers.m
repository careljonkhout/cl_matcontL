function multipliers = NP_MS_compute_multipliers(x, nMults_to_compute)
  global cds contopts;
  print_diag(3,'computing multipliers\n');
  
  [cds.phases_0, period, parameters] = ...
          NP_MS_extract_phases_period_and_parameters(x);
  
  if ( ~ contopts.monodromy_by_finite_differences) && ( ~ cds.using_cvode )
    NP_MS_compute_cycle_parts(x);
  end
  
  M                      = @(x) ...
     NP_MS_full_monodromy_map(1, x, period, parameters);
  nMults_to_compute      = min(nMults_to_compute, cds.n_phases);
  [~, multiplier_matrix] = eigs(M, cds.n_phases, nMults_to_compute);                                    
  multipliers            = diag(multiplier_matrix);
  
  print_diag(1, multipliers2str(multipliers));
end
