function x_corrected = do_one_correction_pd(x0,x,v_cont)
  global cds contopts
  [phi_0, v, T, parameters] = pd_ss_extract_data(x);
  % phi_0 is the current approximation of the reference point on the cycle.
  % v is the current approximation of the eigenvector associated the the
  % eigenvavlue -1
  % T is the current approximation of the period of the cycle.
  % parameters are the parameters value of the two active parameters
  [phi_end, M] = compute_phi_end_and_monodromy(phi_0, T, parameters);
  I            = eye(cds.nphases);
  jacobian     = [M - I; cds.previous_dydt_0'];
  % add d_phi__d_T and d_s__d_T
  d_phi_d_T         = cds.dydt_ode(0, phi_end, parameters{:});
  d_s_d_T           = cds.previous_dydt_0' * d_phi_d_T;
  jacobian          = [jacobian [d_phi_d_T; d_s_d_T]];
  
  rhs = [ phi_end(:) - phi_0(:); % r
         (phi_0 - cds.previous_phases)' * cds.previous_dydt_0 ]; %s
  
  D_x_r__and__D_T_r = - jacobian \ rhs;
  D_x_r = D_x_r__and__D_T_r(1:end-1);
  D_T_r = D_x_r__and__D_T_r(end);
  
  h = contopts.Increment;
  
  M_x_v_x_r = M_x_v(D_x_r, phi_0, v, T, parameters);
  
  
  
  M_v_1 = compute_monodromy(phi_0, T - h, parameters) * v;
  M_v_2 = compute_monodromy(phi_0, T + h, parameters) * v;
  M_T_v = (M_v_2 - M_v_1) / (2 * h);
  
  M_T_v_T_r = M_T_v * D_T_r;
  
 
  D_x_p      = zeros(cds.nphases, 2);
  D_T_p      = zeros(1          , 2);
  M_x_v_x_p  = zeros(cds.nphases, 2);
  M_p_v      = zeros(cds.nphases, 2);
  M_T_v_T_p  = zeros(cds.nphases, 2);
 
  
  
  for i=1:2
    ap_idx     = cds.ActiveParams(i);
    d_phi__d_p = compute_d_phi_d_p(phi_0, T, parameters, ap_idx);
    d_s__d_p   = cds.previous_dydt_0' * d_phi__d_p;
    
    D_x_p__and__D_T_p = jacobian \ [d_phi__d_p; d_s__d_p];
    
    D_x_p(:,i)     = D_x_p__and__D_T_p(1:end-1);
    D_T_p(:,i)     = D_x_p__and__D_T_p(end);
    M_x_v_x_p(:,i) = M_x_v(D_x_p(:,i), phi_0, v, T, parameters);
    
    p1 = parameters;
    p1{ap_idx} = p1{ap_idx} - h;
    M_v_1 = compute_monodromy(phi_0, T, p1) * v;
    p2 = parameters;
    p2{ap_idx} = p2{ap_idx} + h;
    M_v_2 = compute_monodromy(phi_0, T, p2) * v;
    
    M_p_v(:,i) = (M_v_2 - M_v_1) / 2 / h;
    
    M_T_v_T_p(:,i) = M_T_v * D_T_p(i);
  end

  


  
  lhs_1_2 = M_p_v(:,1) + M_x_v_x_p(:,1) + M_T_v_T_p(:,1);
  lhs_1_3 = M_p_v(:,2) + M_x_v_x_p(:,2) + M_T_v_T_p(:,2);
  
  lhs_3_1 = v_cont(cds.nphases + 1 : 2 * cds.nphases)'; 
  
  c_n     = v_cont(1:cds.nphases);
  lhs_3_2 = v_cont(end-1) + c_n' * D_x_p(:,1) + v_cont(end-2) * D_T_p(:,1);
  lhs_3_3 = v_cont(end  ) + c_n' * D_x_p(:,2) + v_cont(end-2) * D_T_p(:,2);
  
  lhs = [ M + I       lhs_1_2     lhs_1_3    ;
          cds.l'      0           0          ;
          lhs_3_1     lhs_3_2     lhs_3_3   ];
 
  

  n = v_cont' * (x - x0);
        
  rhs          = (M + I) * v + M_x_v_x_r + M_T_v_T_r;
  rhs(end + 1) = cds.l' * v - 1;
  rhs(end + 1) = n + c_n' * D_x_r + v_cont(end-2) * D_T_r;
  
  D_v__and__D_p = lhs \ (- rhs);
  
  D_v = D_v__and__D_p(1:cds.nphases);
  D_p = D_v__and__D_p(end-1:end);
  
  D_x = D_x_r + D_x_p(:,1) * D_p(1) + D_x_p(:,2) * D_p(2);
  D_T = D_T_r + D_T_p(:,1) * D_p(1) + D_T_p(:,2) * D_p(2);
  
  
  x_corrected = [
    phi_0    + D_x;
    v        + D_v;
    T        + D_T;
    x(end-1) + D_p(1);
    x(end  ) + D_p(2);];
 
end
function M_x_v_x = M_x_v(x, phi_0, v, period, parameters)
  global contopts
  h = contopts.Increment;
  M_x_v_1 = compute_monodromy(phi_0 - h * x, period, parameters) * v;
  M_x_v_2 = compute_monodromy(phi_0 + h * x, period, parameters) * v;
  M_x_v_x = (M_x_v_2 - M_x_v_1) / (2 * h);
end

%-------------------------------------------------------------------------------
% Computes the monodromy matrix of the current approximation of the cycle 
% starting from the point x near the cycle.
function monodromy = compute_monodromy(phi_0, period, parameters)
  [~, monodromy] = compute_phi_end_and_monodromy(phi_0, period, parameters);
end   

function [phi_T, M] = compute_phi_end_and_monodromy(phi_0, period, parameters)
  global cds
  if cds.using_cvode
    [phi_T, M] = monodromy_cvode(phi_0,period, parameters);
  elseif ~ cds.options.PartitionMonodromy
    [phi_T, M] = monodromy_full(phi_0, period, parameters);
  else
    [phi_T, M] = monodromy_column_by_column(phi_0, period, parameters);
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
function d_phi_d_p = compute_d_phi_d_p(x0, delta_t, parameters, ap_idx)
  h = 1e-5;
  parameters{ap_idx} = parameters{ap_idx} - h;
  phi_1 = NewtonPicard.shoot(x0, delta_t, parameters);
  parameters{ap_idx} = parameters{ap_idx} + 2*h;
  phi_2 = NewtonPicard.shoot(x0, delta_t, parameters);
  d_phi_d_p = (phi_2 - phi_1)/h/2;
  
  % d_phi_d_p_var = NewtonPicard.d_phi_d_p_variational(x0, delta_t, parameters);
  % print_diag(1,'%.6f\n', norm(d_phi_d_p - d_phi_d_p_var));
end
