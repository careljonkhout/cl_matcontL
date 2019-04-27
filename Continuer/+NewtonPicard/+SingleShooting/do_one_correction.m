function x = do_one_correction(x0,x,v0)
  global cds;

  [V, reduced_jacobian, delta_q_gamma, delta_q_r, M_delta_q_r, ...
          phases_0, phi, period, active_par_val] = ...
          NewtonPicard.SingleShooting.compute_reduced_jacobian(x);

  left_hand_side = [
    reduced_jacobian;
    v0(1:end-2)' * V     v0(end-1)     v0(end) - v0(1:end-2)' * delta_q_gamma;
  ];

  % plus of minus M_delta_q_r ?????
  % plus or minus delta_qr ???
  right_hand_side = -[
     V'                   * (phi - phases_0                + M_delta_q_r);
     cds.previous_dydt_0' * (phases_0 - cds.previous_phases  + delta_q_r);
     v0' * (x-x0)          + v0(1:end-2)' * delta_q_r
  ];
  
  delta_p__delta_T_and_delta_gamma = left_hand_side \ right_hand_side;
  delta_p         = delta_p__delta_T_and_delta_gamma(1:end-2);
  delta_T         = delta_p__delta_T_and_delta_gamma(end-1);
  delta_gamma     = delta_p__delta_T_and_delta_gamma(end);

  
  delta_q         = delta_q_r + delta_gamma * delta_q_gamma;
  phases          = phases_0 + V * delta_p + delta_q;
  period          = period + delta_T;
  active_par_val  = active_par_val + delta_gamma;
  x = [phases; period; active_par_val];



  
