function Mx = NP_SS_monodromy_map(x, period, parameters, abs_tol, rel_tol)
  global cds contopts
  if isfield(cds, 'mv_count')
    cds.mv_count = cds.mv_count + 1;
  else
    cds.mv_count = 1;
  end
  if nargin == 3
    abs_tol = contopts.integration_abs_tol;
    rel_tol = contopts.integration_rel_tol;
  elseif nargin ~= 5
    error( ['The number of input arguments to ' ...
      'NP_SS_monodromy_map is not correct.\n' ...
      'The number of input arguments should be either 3 or 5.\n' ...
      'The actual number of input arguments is %d.\n'], nargin);
  end
  if cds.using_cvode
    tic
    [~,~,Mx] = feval(cds.integrator, ...
      't_values',                [0 period], ...
      'initial_point',           cds.phases_0, ...
      'sensitivity_vector',      x, ...
      'ode_parameters',          cell2mat(parameters), ...
      'abs_tol',                 abs_tol, ...
      'rel_tol',                 rel_tol);
    toc
    pause
    return
  end

  if ~ contopts.monodromy_by_finite_differences
    integration_opt = odeset( ...
        'AbsTol',       abs_tol, ...
        'RelTol',       rel_tol,  ...
        'Jacobian',     @(t,y) feval(cds.jacobian_ode, ...
                                t, deval(cds.cycle_orbit,t), parameters{:}));
    dydt_mon = @(t, y) ...
      cds.jacobian_ode(t, deval(cds.cycle_orbit,t), parameters{:}) * y;
    [~, orbit] = ode15s(dydt_mon, [0 period], x, integration_opt);
    Mx = orbit(end, :)';
  else
    % By finite differences: very inaccurate, not reccomended.
    %
    % Below an alternative method of computing the action of the monodromy
    % matrix is implemented. Here, the action of the monodromy matrix is
    % computed by finite differences. 
    
    integration_opt = odeset(...
      'AbsTol',       1e-13, ...
      'RelTol',       1e-13  ...
    ); 
    
    x = x(:);
  
   
    h  = 1e-5;
    f  = @(t, x) cds.dydt_ode(0,x,parameters{:});
    ff = @(t, x1_and_x2) [f(0, x1_and_x2(1:cds.n_phases    ))
                          f(0, x1_and_x2(  cds.n_phases+1:end))];
    x_cycle = deval(cds.cycle_orbit,0);
   
    [~, orbit] = ode15s(ff, ...
      [0 period], ...
      [x_cycle - h * x; x_cycle + h * x], ...
      integration_opt);
    
    phi_x1__and__phi_x2 = orbit(end,:)';
    phi_x1 = phi_x1__and__phi_x2(1:cds.n_phases);
    phi_x2 = phi_x1__and__phi_x2(cds.n_phases+1:end);
    
    Mx = ((phi_x2 - phi_x1)/h)/2; 
  end
end