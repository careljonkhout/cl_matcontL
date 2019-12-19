% Computes the Jacobian matrix with subspace reduced applied to the phase space
% also returns additional values that were used in the computation of the
% reduced Jacobian matrix, which are needed in do_one_correction.m and/or 
% find_tangent_vector.m. The subspace used is the subspace associated to the
% largest eigenvalues (in modulus) of the mondromy matrix.

function [V, reduced_jacobian, delta_q_gamma, delta_q_r, M_delta_q_r, ...
          phases_0, phi, period, active_par_val] = NP_SS_compute_reduced_jacobian(x)

  global cds contopts
  
  % if no options have been set yet, use default options
  if isempty(contopts)
    contopts = contset;
  end
  
  [phases_0, period, parameters] = ...
        NP_SS_extract_phases_period_and_parameters(x);
  active_par_val = parameters{cds.ActiveParams};
 
  using_cvode = cds.using_cvode;
  

  if using_cvode
 
    [~,y] = cds.integrator( ...
      't_values',                [0 period], ...
      'initial_point',           phases_0, ...
      'ode_parameters',          cell2mat(parameters), ...
      'abs_tol',                 contopts.integration_abs_tol, ...
      'rel_tol',                 contopts.integration_rel_tol);
    
    phi = y(end,:)';
    cds.phases_0 = phases_0;
  else
    % Since the result of the next time-integration will be reused many times,
    % we set the tolerances a bit tighter than the rest.
    integration_opt = odeset(...
      'AbsTol',      contopts.integration_abs_tol / 10,    ...
      'RelTol',      contopts.integration_rel_tol / 10     ...
    );
    if ~ isempty(cds.jacobian_ode)
      integration_opt = odeset(integration_opt, ...
      'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}));
    end
    cds.cycle_orbit = cds.integrator(...
      @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
      [0 period], ...
      phases_0, integration_opt);

    phi = deval(cds.cycle_orbit,period);
  end
  
 
  cds.V = NP_SS_compute_subspace(period, parameters);
  V = cds.V;
  basis_size = size(V,2);


  MV = zeros(cds.n_phases,basis_size);
  
  integrator = cds.integrator;
  rel_tol = contopts.integration_rel_tol;
  abs_tol = contopts.integration_abs_tol;
  if contopts.contL_ParallelComputing
    try 
      parfor i=1:basis_size
        % The function monodromy_map cannot be used here, since it depends on
        % the global variable cds, and global variables are not copied so the
        % the workspace of the workers that parfor uses.
        if using_cvode
          [~,~,y] = feval(integrator, ...
            't_values',                t_values, ...
            'sensitivity_vector',      V(:,i), ...
            'ode_parameters',          cell2mat(parameters), ...
            'abs_tol',                 abs_tol, ...
            'rel_tol',                 rel_tol);
          MV(:,i) = y(end,:)';
        else
     %     [~,orbit] = feval(integrator, ...
     %       dydt_monodromy_map, [0 period], V(:,i), int_opt);

     %     MV(:,i) = orbit(end,:)';
        end
      end
      parfor_failed = false;
    catch error
      if (strcmp(error.identifier,'MATLAB:remoteparfor:AllParforWorkersAborted'))
        % Something went wrong with the parfor workers.
        % We try again with ordinary for.
        fprintf('Parfor aborted, retrying with ordinary for.\n');
        parfor_failed = true;
      else
        % in case of some other error, we want to know about it
        rethrow(error)
      end
    end
  end
  
  if (~ contopts.contL_ParallelComputing) || parfor_failed
    for i=1:basis_size
      MV(:,i) = NP_SS_monodromy_map(V(:,i), period, ...
                    parameters);
    end
  end
  
  

  S = V'*MV;
  
  
  % the r in q_r means residual
  [delta_q_r,     M_delta_q_r] = ...
    NP_SS_solve_Q_system(V, ( phi - phases_0), ...
    period, parameters);

  d_phi_d_gamma_val = NP_compute_d_phi_d_p( ...
                            phases_0,period,parameters);
                          
  [delta_q_gamma, M_delta_q_gamma] = ...
    NP_SS_solve_Q_system(V, d_phi_d_gamma_val, ...
    period, parameters);
  

  
  d_phi_d_T = cds.dydt_ode(0,phi,parameters{:});
  

  d_s_d_T      = 0;
  d_s_d_gamma  = 0;

  
  reduced_jacobian = [
    S - eye(basis_size)         V'* d_phi_d_T      V' * (d_phi_d_gamma_val + M_delta_q_gamma);
    cds.previous_dydt_0' * V    d_s_d_T            d_s_d_gamma + cds.previous_dydt_0' * delta_q_gamma                               ;
  ];
  
