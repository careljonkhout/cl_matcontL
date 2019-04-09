function x = do_one_correction(x0,x,v0)
  global cds;

  [~, reduced_jacobian, phases_0, phi, period, active_par_val] = ...
         PartialMonodromy.compute_reduced_jacobian(x);

  left_hand_side = [
    reduced_jacobian;
    v0';
  ];

  right_hand_side = [
    - (phi - phases_0);
    - cds.previous_dydt_0' * phases_0;
    - v0' * (x-x0)          
  ];
  
  delta_phases__delta_T_and_delta_gamma = left_hand_side \ right_hand_side;
  delta_phases    = delta_phases__delta_T_and_delta_gamma(1:end-2);
  delta_T         = delta_phases__delta_T_and_delta_gamma(end-1);
  delta_gamma     = delta_phases__delta_T_and_delta_gamma(end);
  
  phases          = phases_0 + delta_phases;
  period          = period + delta_T;
  active_par_val  = active_par_val + delta_gamma;
  x = [phases; period; active_par_val];



  
