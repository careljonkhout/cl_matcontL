%-------------------------------------------------------------------------------
% Computes the solution of the initial value problem x(0) = x, x'(t) = f(x),
% where f is defined by cds.dydt_ode, and returns x(period). The initial value
% problem is solved using solver from the matlab ode suite or a fully compatible
% alternative. The specific solver is specified by cds.integrator.
function x_end = shoot(x, delta_t, parameters)
  global cds contopts
  integration_opt = odeset(...
    'AbsTol',      contopts.integration_abs_tol,    ...
    'RelTol',      contopts.integration_rel_tol     ...
  );
  if ~ isempty(cds.jacobian_ode)
    integration_opt = odeset(integration_opt, ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}));
  end
  f = @(t, y) cds.dydt_ode(t, y, parameters{:});
  [~, orbit] = cds.integrator(f, [0 delta_t], x, integration_opt);
  x_end = orbit(end,:)';