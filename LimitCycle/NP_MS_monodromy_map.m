%
% i:                1-based index of the mesh point from where to start
%
% x:                vector to which the monodromy map is applied
%
% time_interval:    length of the time interval for time integration
%                   (time_interval should be equal to the period of the cycle
%                   if the full monodromy map is to be evaluated. In multiple
%                   shooting partial monodromy, the value of time_interval
%                   passed to this method can be less than the period of the
%                   cycle, since "partial" monodromy maps must also be
%                   evaluated during a multiple shooting continuation.)
%
% parameters:       cell array of parameters for the system of ODEs
%
% abs_tol:          absolute tolerance for the time integration. This argument
%                   is optional. However, either both abs_tol and rel_tol must
%                   be specified, or neither.
%
% rel_tol:          absolute tolerance for the time integration. This argument
%                   is optional. However, either both abs_tol and rel_tol must
%                   be specified, or neither.
function Mx  = NP_MS_monodromy_map(i, x, time_interval, parameters, abs_tol, rel_tol)
  global cds contopts
  if nargin == 4
    abs_tol = contopts.integration_abs_tol;
    rel_tol = contopts.integration_rel_tol;
  elseif nargin ~= 6
    error( ['The number of input arguments to ' ...
      'NP_MS_monodromy_map is not correct.\n' ...
      'The number of input arguments should be either 4 or 6.\n' ...
      'The actual number of input arguments is %d.\n'], nargin);
  end
  int_opt = odeset(...
    'AbsTol',       abs_tol,    ...
    'RelTol',       rel_tol,    ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode, ...
                      t, deval(cds.orbits(i),t), parameters{:}) ...
  );
  if cds.using_cvode
    [~,~,Mx] = feval(cds.integrator, ...
      't_values',                [0 time_interval], ...
      'initial_point',           cds.phases_0(:,i), ...
      'sensitivity_vector',      x, ...
      'ode_parameters',          cell2mat(parameters), ...
      'abs_tol',                 abs_tol, ...
      'rel_tol',                 rel_tol);
    return
  end

  if ~ contopts.monodromy_by_finite_differences
    dydt_mon = @(t, y) ...
      cds.jacobian_ode(t, deval(cds.orbits(i), t), parameters{:}) * y;

    [~,orbit] = cds.integrator(dydt_mon, [0 time_interval], x, int_opt);

    Mx = orbit(end,:)';
  else 
    % Below an alternative method of computing the action of the monodromy
    % matrix is implemented. Here, the action of the monodromy matrix is
    % computed by finite differences. Seems to be faster than the method above.
    % The error of the trivial multiplier is much larger when using finite
    % differences.
    %
    % This method of computing the action of the monodromy matrix also makes
    % continuation possible is the Jacobian of the system of ODEs is not
    % available
    
    integration_opt = odeset(...
      'AbsTol',       1e-13, ...
      'RelTol',       1e-13  ...
    ); 
    h = 5e-5;
    x_cycle = deval(cds.orbits(i),0);
    f  = @(t, x) cds.dydt_ode(0,x,parameters{:});
    ff = @(t, x1_and_x2) [f(0, x1_and_x2(1:cds.n_phases    ))
                          f(0, x1_and_x2(  cds.n_phases+1:end))];
    [~, orbit] = ode15s(ff, ...
      [0 time_interval], ...
      [x_cycle - h * x; x_cycle+h * x], ...
      integration_opt);
    
    phi_x1__and__phi_x2 = orbit(end,:)';
    phi_x1 = phi_x1__and__phi_x2(1:cds.n_phases);
    phi_x2 = phi_x1__and__phi_x2(cds.n_phases+1:end);
    
    Mx = (phi_x2 - phi_x1)/h/2; 
  end
end