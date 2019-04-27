function d_phi_d_p = d_phi_d_p_variational(x, period, parameters)
  f          = @(t,w) dydt(t,w,parameters);
  [~, orbit] = ode15s(f, [0 period],zeros(length(x),1));
  d_phi_d_p  = orbit(end,:)';
end

function dydt = dydt(t, w, parameters)
  global cds
  phi_t      = deval(cds.cycle_orbit,t);
  jacobian   = cds.jacobian_ode  (t, phi_t, parameters{:});
  p_cols     = zeros(length(cds.P0),1);
  p_cols(cds.ActiveParams) = 1;
  jacobian_p = cds.jacobian_p_ode(t, phi_t, parameters{:}) * p_cols;
  dydt       = jacobian * w + jacobian_p;
end

