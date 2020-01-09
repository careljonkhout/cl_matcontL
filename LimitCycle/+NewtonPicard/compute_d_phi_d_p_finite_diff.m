function d_phi_d_p = compute_d_phi_d_p_finite_diff(x0, delta_t, parameters)
  global cds
  ap = cds.ActiveParams;
  h = 1e-6;
  parameters{ap} = parameters{ap} - h;
  phi_1 = NewtonPicard.shoot(x0, delta_t, parameters);
  parameters{ap} = parameters{ap} + 2*h;
  phi_2 = NewtonPicard.shoot(x0, delta_t, parameters);
  d_phi_d_p = (phi_2 - phi_1)/h/2;
end