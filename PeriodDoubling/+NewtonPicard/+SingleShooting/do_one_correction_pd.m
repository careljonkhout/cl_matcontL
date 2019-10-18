function do_one_correction_pd(x0,x,v)
  [phases_0, v, period, parameters] = pd_ss_extract_data(x);
  [phases_end, monodromy] = compute_monodromy(phases_0, period, parameters);
  jacobian          = [monodromy-eye(cds.nphases); cds.previous_dydt_0'];
  % add d_phi__d_T and d_s__d_T
  d_phi_d_T         = cds.dydt_ode(0, phases_end, parameters{:});
  d_s_d_T           = cds.previous_dydt_0' * d_phi_d_T;
  jacobian          = [jacobian [d_phi_d_T; d_s_d_T]];
  
  rhs = [ phases_end(:) - phases_0(:); % r
          phases_0 - cds.previous_phases)' * cds.previous_dydt_0 ] %s
  
  delta_x_r__and__delta_T_r = - jacobian \ rhs;
  delta_x_r = delta_x_r__and__delta_T_r(1:end-1);
  detla_T_r = delta_x_r__and__delta_T_r(end);
  
  
  
  

    
  compute_d_phi_d_p = @NewtonPicard.compute_d_phi_d_p;
  d_phi__d_p        = compute_d_phi_d_p(phases, period, parameters, p_idx);
  d_s__d_p          = cds.previous_dydt_0' * d_phi__d_p;

  delta_x_gamma__and__delta_T_gamma = jacobian \ [d_phi__d_p; d_s__d_p];
  delta_x_gamma = delta_x_gamma + delta_x_gamma__and__delta_T_gamma(1:end-1);
  delta_T_gamma = delta_T_gamma + delta_x_gamma__and__delta_T_gamma(end);

  M_v_1 = compute_monodromy(phases_0, period, parameters) * v;
  M_v_2 = compute_monodromy(phases_0, period, parameters) * v;
  
  M_v_1 = compute_monodromy(phases_0 - h * delta_x_gamma, period, parameters) * v;
  M_v_2 = compute_monodromy(phases_0 + h * delta_x_gamma, period, parameters) * v;
  M_x_v_x = (M_v_2 - M_v_1) / (2 * h);
  
  M_v_1 = compute_monodromy(phases_0, period - h, parameters) * v;
  M_v_2 = compute_monodromy(phases_0, period + h, parameters) * v;
  M_T_v_T = (M_v_2 - M_v_1) / (2 * h) * delta_T_gamma;
  
  

 
end
function M_x_v_x = compute_M_x_v(x, phases_0, v, period, parameters);
  h = 1e-5;
  M_x_v_1 = compute_monodromy(phases_0 - h * x, period, parameters) * v;
  M_x_v_2 = compute_monodromy(phases_0 + h * x, period, parameters) * v;
  M_x_v_x = (M_x_v_2 - M_x_v_1) / (2 * h);
end
  
%-------------------------------------------------------------------------------
% Computes the monodromy matrix of the current approximation of the cycle 
% starting from the point x near the cycle.
function [y_end, monodromy] = compute_monodromy(x, period, parameters)
  global cds
  if cds.using_cvode
    [y_end, monodromy] = monodromy_cvode(x,period, parameters);
  elseif ~ cds.options.PartitionMonodromy
    [y_end, monodromy] = monodromy_full(x, period, parameters);
  else
    [y_end, monodromy] = monodromy_column_by_column(x, period, parameters);
  end
end 
%-------------------------------------------------------------------------------
% Computes the monodromy matrix of the current approximation of the cycle
% starting from the point x near the cycle. The matrix is computed by
% integrating one N + N^2 dimensional system defined by:
%
% x'(t) = f(x(t), parameters)
% M'(t) = f_x(x(t), parameters) * M(t)
% x(0)  = x_0
% M(0)  = I
%
% where M is the monodromy matrix and N is the number of 1 dimensional equations
% in the system of ODE's. My (Carel Jonkhout) estimation, based on a couple of
% experiments, is that this is only efficient for N upto about 10. (This is
% probably system dependent) Otherwise, is seems to be more efficient to
% integrate x and each column of M separately.
%
% This method is for testing purposes only. For actual continuation of cycles it
% is recomended to use NewtonPicard or orthogonal collocation.
function [y_end, monodromy] = monodromy_full(x_0, period, parameters)
  global cds contopts
  nphases = cds.nphases;
  f =@(t, y) dydt_monodromy_full(t, y, parameters);
  integration_opt = odeset(...
    'AbsTol',      contopts.integration_abs_tol,    ...
    'RelTol',      contopts.integration_abs_tol     ... % todo add JPattern
  );

  x_with_monodromy = [x_0; reshape(eye(nphases),[nphases^2 1])];
  [~, trajectory] = cds.integrator(...
    f, [0 period], x_with_monodromy, integration_opt);
  y_end = trajectory(end,1:nphases)';
  monodromy = trajectory(end,nphases+1:end);
  monodromy = reshape(monodromy, [nphases nphases]);
end
%-------------------------------------------------------------------------------
function dydt_mon = dydt_monodromy_full(t,y, parameters)
  global cds
  y_ode = y(1:cds.nphases);
  
  y_mon = reshape(y(cds.nphases+1:end),cds.nphases,cds.nphases);
  dydt_mon = [
      cds.dydt_ode(t, y_ode, parameters{:}); 
      reshape( ...
        cds.jacobian_ode(t, y_ode, parameters{:}) * y_mon, ...
        [cds.nphases^2 1]) 
  ];
end
%-------------------------------------------------------------------------------   
function [y_end, monodromy] = monodromy_column_by_column(x, period, parameters)
  global cds contopts;
  integration_opt = odeset(...
    'AbsTol',      contopts.integration_abs_tol,    ...
    'RelTol',      contopts.integration_rel_tol,    ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  f = @(t, y) cds.dydt_ode(t, y, parameters{:});
  cycle = cds.integrator(f, [0 period], x, integration_opt);
  y_end = deval(cycle,period);
  monodromy = eye(cds.nphases);
  integration_opt = odeset(integration_opt, 'Jacobian', ...
    @(t,y) feval(cds.jacobian_ode, t, deval(cycle,t), parameters{:}));
  f = @(t, y) cds.jacobian_ode(t, deval(cycle,t), parameters{:}) * y;
  integrator = cds.integrator;
  for i=1:cds.nphases
    fprintf('%d ',i);
    [~, monodromy_map_trajectory] = feval(integrator, ...
      f, [0 period], monodromy(:,i),integration_opt);
    monodromy(:,i) = monodromy_map_trajectory(end,:);
  end 
end
function [y_end, M] = monodromy_cvode(x, period, parameters)
  global cds contopts
  M = zeros(cds.nphases);
  for i=1:cds.nphases
    e_i = zeros(cds.nphases,1);
    e_i(i) = 1;
    [~,y,M(:,i)] = feval(cds.integrator, ...
      't_values',                [0 period], ...
      'initial_point',           x, ...
      'sensitivity_vector',      e_i, ...
      'ode_parameters',          cell2mat(parameters), ...
      'abs_tol',                 contopts.integration_abs_tol, ...
      'rel_tol',                 contopts.integration_rel_tol);
  end
  y_end = y(end,:)';
end

%-------------------------------------------------------------------------------
% Computes the derivative of the solution of the problem x'(t) = f(x), x(0) = x0
% w.r.t. the active parameter, evaluated at x.
% todo: use variational method instead for improved accuracy
% todo: implement computation of d_phi_d_p by cvodes
function d_phi_d_p = compute_d_phi_d_p(x0, delta_t, parameters)
  global cds
  p_idx_1 = cds.ActiveParams(1);
  p_idx_2 = cds.ActiveParams(2);
  h = 1e-5;
  parameters{p_idx_1} = parameters{p_idx_1} - h;
  parameters{p_idx_2} = parameters{p_idx_2} - h;
  phi_1 = NewtonPicard.shoot(x0, delta_t, parameters);
  parameters{p_idx_1} = parameters{p_idx_1} + 2*h;
  parameters{p_idx_2} = parameters{p_idx_2} + 2*h;
  phi_2 = NewtonPicard.shoot(x0, delta_t, parameters);
  d_phi_d_p = (phi_2 - phi_1)/h/2;
  
  % d_phi_d_p_var = NewtonPicard.d_phi_d_p_variational(x0, delta_t, parameters);
  % print_diag(1,'%.6f\n', norm(d_phi_d_p - d_phi_d_p_var));
end
