%-------------------------------------------------------------------------------
% Computes the derivative of the solution of the problem x'(t) = f(x), x(0) = x0
% w.r.t. the active parameter, evaluated at x.
function d_phi_d_p = compute_d_phi_d_p(x0, delta_t, parameters)
  global cds
  ap = cds.ActiveParams;
  h = 1e-6;
  parameters{ap} = parameters{ap} - h;
  phi_1 = NewtonPicard.shoot(x0, delta_t, parameters);
  parameters{ap} = parameters{ap} + 2*h;
  phi_2 = NewtonPicard.shoot(x0, delta_t, parameters);
  d_phi_d_p = (phi_2 - phi_1)/h/2;
  
  % d_phi_d_p_var = NewtonPicard.d_phi_d_p_variational(x0, delta_t, parameters);
  % print_diag(1,'%.6f\n', norm(d_phi_d_p - d_phi_d_p_var));
end