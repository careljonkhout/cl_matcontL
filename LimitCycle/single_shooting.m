% Curve file of cycle continuation by single shooting
function out = single_shooting
  out{1}  = @curve_function;
  out{2}  = @default_processor;
  out{3}  = @options;
  out{4}  = @jacobian;
  out{5}  = [];%@hessians;
  out{6}  = @testfunctions;
  out{7}  = [];%@userf;
  out{8}  = @process_singularity;
  out{9}  = @cycle_singularity_matrix;
  out{10} = @locate;
  out{11} = @init;
  out{12} = [];%@done;
  out{13} = @adapt;
  out{14} = @curve_CIS_first_point;
  out{15} = @curve_CIS_step;
end
%-------------------------------------------------------------------------------
% Computes the curve function. if curve_function(x) == 0, then x corresponds to
% a point on the limitcycle. This point will be near the point from the previous
% continuation step, due to the phase condition. The phase condition is that the
% next point must lie on the plane P through the current point x, where P is
% perpendicular to the tangent vector to the limit cycle at x. The phase
% condition is necessary to uniquely define the next point, and is used to limit
% the phase shift from one step to the next.
function f = curve_function(varargin)
  global cds
  x = varargin{1};
  active_par_val               = x(end);
  period                       = x(end-1);
  phases_0                     = x(1:end-2);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = active_par_val;
  parameters                   = num2cell(parameters);
  shoot                        = @NewtonPicard.shoot;
  phases_end                   = shoot(phases_0, period, parameters);
  f = [phases_end - phases_0; 
      (phases_0 - cds.previous_phases)' * cds.previous_dydt_0  ]; 
end
%-------------------------------------------------------------------------------
% Computes the Jacobian matrix of the curvefunction at evaluated at varargin{1}.
function jacobian = jacobian(varargin)
  global cds
  cont_state                   = varargin{1};
  active_par_val               = cont_state(end);
  period                       = cont_state(end-1);
  phases                       = cont_state(1:end-2);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = active_par_val;
  parameters                   = num2cell(parameters);
  
  
  % 
  
  global contopts lds;
  integration_opt = odeset(...
    'AbsTol',      contopts.integration_abs_tol,    ...
    'RelTol',      contopts.integration_rel_tol,    ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  f = @(t, y) cds.dydt_ode(t, y, parameters{:});
  cycle = cds.integrator(f, [0 period], phases, integration_opt);
%  x_for_monodromy = zeros((lds.ntst*lds.ncol + 1) * lds.nphase + 2,1);
%  mesh = linspace(0,period,lds.ntst*lds.ncol + 1);
%  x_indices = 1:lds.nphase;
%  for t = mesh
%    x_for_monodromy(x_indices) = deval(cycle, t);
%    x_indices = x_indices + lds.nphase;
%  end
%  x_for_monodromy(end-1) = period;
%  x_for_monodromy(end  ) = active_par_val;
  % an experimental way of computing the monodromy matrix:
 % monodromy = compute_monodromy_oc(x_for_monodromy);
  %y_end = deval(cycle,period);
  
  [y_end, monodromy] = compute_monodromy(phases, period, parameters);
  
  
  
  
  
  jacobian     = [monodromy-eye(cds.nphases); cds.previous_dydt_0'];
  % add d_phi__d_T and d_s__d_T
  d_phi_d_T         = cds.dydt_ode(0,y_end,parameters{:});
  d_s_d_T           = cds.previous_dydt_0' * d_phi_d_T;
  jacobian          = [jacobian [d_phi_d_T; d_s_d_T]];
  compute_d_phi_d_p = @NewtonPicard.compute_d_phi_d_p;
  d_phi__d_p        = compute_d_phi_d_p(phases, period, parameters);
  d_s__d_p          = cds.previous_dydt_0' * d_phi__d_p;
  jacobian          = [jacobian [d_phi__d_p; d_s__d_p]];
end

%-------------------------------------------------------------------------------
% Computes the monodromy matrix of the current approximation of the cycle 
% starting from the point x near the cycle.
function [y_end, monodromy] = compute_monodromy(x, period, parameters)
  global cds
  if ~ cds.options.PartitionMonodromy
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
%-------------------------------------------------------------------------------
function init(~,~); end
%-------------------------------------------------------------------------------
function out = default_processor(varargin)
  global cds contopts
  point = varargin{1};
  
  basis_size_changed = false;
  
  if contopts.NewtonPicard
    update_multipliers_if_needed(point.x)
    if abs(cds.multipliers(end)) > contopts.basis_grow_threshold
      basis_size_changed = true;
      print_diag(2, 'expanding basis\n');
      nMults_to_compute = cds.preferred_basis_size + 10;
      nMults_to_compute = min(nMults_to_compute, cds.nphases);
      cds.multipliersX = point.x;
      cds.multipliers = NewtonPicard.SingleShooting.compute_multipliers(...
        point.x, nMults_to_compute);
      
      i = length(cds.multipliers);
     
      while abs(cds.multipliers(i)) < contopts.basis_grow_threshold / 3
        i = i - 1;
      end
      if i < length(cds.multipliers)
        i = i + 1;
      end
      cds.preferred_basis_size = i;
      cds.p                    = i;
      cds.multipliers = cds.multipliers(1:i);
    elseif abs(cds.multipliers(end)) < contopts.basis_shrink_threshold
      basis_size_changed = true;
      i = length(cds.multipliers);
      while abs(cds.multipliers(i)) < contopts.basis_shrink_threshold
        i = i - 1;
      end
      cds.preferred_basis_size = i;
      cds.p = i;
      cds.multipliers = cds.multipliers(1:i);
    end
  end
  
  if contopts.Multipliers || contopts.NewtonPicard
    update_multipliers_if_needed(point.x)
    point.multipliers = cds.multipliers;
  end
  
  if contopts.Singularities && basis_size_changed
    point.tvals = testfunctions(cds.ActTest,point.x,point.v,[]);
    point.tvals = testfunctions(cds.ActTest,point.x,point.v,[]);
    test_function_labels = {'BPC', 'PD', 'LPC', 'NS'};
    print_diag(1,'testfunctions: [')
    for i=1:length(point.tvals)
      print_diag(1,' %s:',test_function_labels{i});
      print_diag(1,'%+.5e',point.tvals(i))
    end
    print_diag(1,']\n')
  end
  
  
  x = point.x;
  % needed for phase condition:
  
  
  [phases, ~, parameter_values] = getComponents(x);
  
  cds.previous_phases           = phases;
  % needed for phase condition:
  cds.previous_dydt_0           = cds.dydt_ode(0, ...
                                    cds.previous_phases,parameter_values{:});
  
  out                           = point;
  point.paramters_values        = parameter_values;
  savePoint(point, varargin{2:end});
end
%-------------------------------------------------------------------------------
function update_multipliers_if_needed(x)
  global cds
  if ~ isfield(cds,'multipliersX') || all(cds.multipliersX ~= x)
    cds.multipliersX = x;
    cds.multipliers = NewtonPicard.SingleShooting.compute_multipliers(x, ...
                       cds.preferred_basis_size);
  end
end
%-------------------------------------------------------------------------------
% Test functions are used for detecting AND locating singularities by bisection.
% When detecting ids_testf_requested will be cds.ActTest, and when locating
% ids_testf_requested will contain only those ids of the testfunctions relevant
% to the bifurcation that is being located.
function [out, failed] = testfunctions(ids_testf_requested, x0, v, ~) 
  % unused arguments are v and CISdata
  global cds;
  failed = false;
  
  const = Constants;
  
  if any(ismember([const.BPC_id const.PD_id const.NS_id], ids_testf_requested))
    update_multipliers_if_needed(x0)
  end
  
  out = cycle_testfunctions(ids_testf_requested, cds.multipliers, v);
end
%-------------------------------------------------------------------------------
% After a singularity is detected and located, the contL calls this function.
% Usually a normal form coefficient (nfc) is computed, via a function call from
% this function. However, for continuation of cycles by single shooting or
% multiple shooting, the computation of normal form coefficients is not
% implemented yet.
function [failed,s] = process_singularity(id,point,s)
  global cds
  x = point.x;
  switch id
  case 1
    format_string = 'Branch Point cycle(period = %e, parameter = %e)\n'; 
    print_diag(0, format_string, x(end-1), x(end));
    s.msg  = sprintf('Branch Point cycle'); 
  case 2
    format_string = 'Period Doubling (period = %e, parameter = %e)\n';
    print_diag(0, format_string, x(end-1), x(end));
    s.msg  = 'Period Doubling';
  case 3
    s.msg = 'Limit point cycle';
    format_string = 'Limit point cycle (period = %e, parameter = %e)\n';
    print_diag(0, format_string, x(end-1), x(end));
  case 4
    d = cds.multipliers;
    smallest = Inf;
    % we find the pair of multipliers whose product is closest to one, and check
    % to see these the multipliers have a nonzero imaginary part.
    for jk=1:length(d)-1
      [val,idx] = min(abs(d(jk+1:length(d))*d(jk)-1));
      if val < smallest
        idx2 = jk+idx;
        smallest = val;
      end
    end
    threshold = cds.deviation_of_trivial_multiplier;
    singularity_is_neutral_saddle = abs(imag(d(idx2))) < threshold;
    if singularity_is_neutral_saddle
      s.msg = 'Neutral saddle cycle';
      format_string = 'Neutral Saddle Cycle (period = %e, parameter = %e)\n';
      % A neutral saddle is not really a bifurcation
      print_diag(0, format_string, x(end-1), x(end));
    else
      s.msg = 'Neimark Sacker';
      format_string = 'Neimark-Sacker (period = %e, parameter = %e)\n';
      print_diag(0, format_string, x(end-1) ,x(end));
      %print_diag(0, 'Normal form coefficient = %d\n', s.data.nscoefficient);
    end
  end
  failed = 0;
end  
%-------------------------------------------------------------------------------
function options
end
%-------------------------------------------------------------------------------
function CISdata = curve_CIS_first_point(~) % unused argument is x
  CISdata = 1;
end
%-------------------------------------------------------------------------------
function CISdata = curve_CIS_step(~, ~) 
  % unused arguments are x an CISdata_in
  CISdata = 1;
end
%-------------------------------------------------------------------------------
function [has_changed,x2,v2,CISdata] = adapt(x,v,~,~)
  % unused arguments are CISdata, and tfUpdate
  x2 = x;
  v2 = v;
  CISdata = [];
  has_changed = false;
end
%-------------------------------------------------------------------------------
function p_out = locate(id, p1, p2) %#ok<STOUT,INUSD>
  switch id   
    otherwise
      error('No locator defined for singularity %d', id);
  end
end
%-------------------------------------------------------------------------------
function [y,period,parameters] = getComponents(x)
  global cds
  y                            = x(1:cds.nphases);
  period                       = x(end-1);
  parameter_value              = x(end);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = parameter_value;
  parameters                   = num2cell(parameters);
end
