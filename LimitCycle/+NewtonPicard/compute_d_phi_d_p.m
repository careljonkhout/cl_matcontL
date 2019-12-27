function d_phi_d_p = compute_d_phi_d_p(x, delta_t, parameters)
  global cds contopts
  persistent warning_given
  
  if contopts.parameter_sensitivity_by_finite_diff
    d_phi_d_p = compute_d_phi_d_p_finite_diff(x, delta_t, parameters);
    return
  end
  
  if isempty(cds.jacobian_p_ode)
    if ~ warning_given
      warning(['Jacobian of ODE system w.r.t. parameters not found. ' ...
               'using finite differences for parameter sensitivity']);
      warning_given = true;
    end
    d_phi_d_p = compute_d_phi_d_p_finite_diff(x, delta_t, parameters);
    return
  end
    
  if cds.using_cvode
    [~, ~, d_phi_d_p] = cds.integrator( ...
      'initial_point',         x, ...
      'ode_parameters',        cell2mat(parameters), ...
      't_values',              [0 delta_t], ...
      'parameter_sensitivity', cds.ActiveParams - 1);

  else
    f          = @(t,w) dydt(t,w,parameters);
    [~, orbit] = ode15s(f, [0 delta_t],zeros(length(x),1));
    d_phi_d_p  = orbit(end,:)';
  end
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

function d_phi_d_p = compute_d_phi_d_p_finite_diff(x0, delta_t, parameters)
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

