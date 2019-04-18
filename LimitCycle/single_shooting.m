% Curve file of cycle continuation by single shooting
function out = single_shooting
  out{1}  = @curve_function;
  out{2}  = @defaultprocessor;
  out{3}  = @options;
  out{4}  = @jacobian;
  out{5}  = [];%@hessians;
  out{6}  = @testfunctions;
  out{7}  = [];%@userf;
  out{8}  = @process_singularity;
  out{9}  = @singularity_matrix;
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
% Computes the jacobian of the curvefunction at evaluated at varargin{1}.
% This function is not accesible from outside this file.
function jacobian = jacobian(varargin)
  global cds
  cont_state                   = varargin{1};
  active_par_val               = cont_state(end);
  period                       = cont_state(end-1);
  phases                       = cont_state(1:end-2);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = active_par_val;
  parameters                   = num2cell(parameters);
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
% x'    = f(x, parameters)
% x(0)  = x
% M'(x) = f_x(x, parameters) * M
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
function [y_end, monodromy] = monodromy_full(x, period, parameters)
  global cds contopts
  nphases = cds.nphases;
  f =@(t, y) dydt_monodromy_full(t, y, parameters);
  integration_opt = odeset(...
    'AbsTol',      contopts.integration_abs_tol,    ...
    'RelTol',      contopts.integration_abs_tol     ... % todo add JPattern
  );

  x_with_monodromy = [x; reshape(eye(nphases),[nphases^2 1])];
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
  cycle = cds.integrator(f, linspace(0,period,cds.nDiscretizationPoints), ...
    x, integration_opt);
  y_end = deval(cycle,period);
  monodromy = eye(cds.nphases);
  integration_opt = odeset(integration_opt, 'Jacobian', ...
    @(t,y) feval(cds.jacobian_ode, t, deval(cycle,t), parameters{:}));
  f = @(t, y) cds.jacobian_ode(t, deval(cycle,t), parameters{:}) * y;
  integrator = cds.integrator;
  parfor i=1:cds.nphases
    fprintf('%d ',i);
    [~, monodromy_map_trajectory] = feval(integrator, ...
      f, [0 period], monodromy(:,i),integration_opt);
    monodromy(:,i) = monodromy_map_trajectory(end,:);
  end 
end
%-------------------------------------------------------------------------------
function init(~,~); end
%-------------------------------------------------------------------------------
function out = defaultprocessor(varargin)
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
  
  if basis_size_changed
    point.tvals = testfunctions(cds.ActTest,point.x,point.v,[]);
    print_diag(1,'Test Functions: [')
    print_diag(1,' %+.5e',point.tvals)
    print_diag(1,']\n')
  end
  
  
  x = point.x;
  % needed for phase condition:
  cds.previous_phases          = x(1:end-2);
  
  active_par_val               = x(end);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = active_par_val;
  parameters                   = num2cell(parameters);
  
  % needed for phase condition:
  cds.previous_dydt_0          = cds.dydt_ode(0, ...
                                    cds.previous_phases,parameters{:});
  
  out                          = point;
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
% Test functions are used for detecting AND location singularities by bisection.
% When detecting ids_testf_requested will be cds.ActTest, and when locating
% ids_testf_requested will contain only those ids of the testfunctions relevant
% to the bifurcation that is being located.
function [out, failed] = testfunctions(ids_testf_requested, x0, v, ~) 
  % unused arguments are v and CISdata
  global cds contopts
  
  failed = false;
  
  const = Constants;
  
  if any(ismember([const.BPC_id const.PD_id const.NS_id], ids_testf_requested))
    update_multipliers_if_needed(x0)
  end
  if any(ismember(const.BPC_id, ids_testf_requested))
    mults = cds.multipliers;
    distance_to_one = abs(mults - 1);
    [~,index]       = min(distance_to_one);
    mults(index) = 0;
    out(const.BPC_id) = real(prod(mults - ones(size(mults))));
  end
  if ismember(const.PD_id, ids_testf_requested)
    % real is needed to ensure that small complex parts induced by roundoff
    % errors do not make the result complex.
    % A complex valued test function would cause false positives when detecting 
    % bifurcations.
    out(const.PD_id) = real(prod(cds.multipliers + ones(size(cds.multipliers))));
  end
  if ismember(const.LPC_id, ids_testf_requested)
    out(const.LPC_id) = v(end);
  end
  if any(ismember(const.NS_id, ids_testf_requested))
    mults                  = cds.multipliers;
    threshold              = contopts.real_v_complex_threshold;
    complex_mults          = mults(abs(imag(mults)) > threshold);
    unstable_complex_mults = complex_mults(abs(complex_mults) > 1); 
    out(const.NS_id)       = length(unstable_complex_mults);
  end
  
end
%-------------------------------------------------------------------------------
% defines which changes in testfunctions correspond to which singularity type
% 0 == require sign-change
% 1 == require sign-non-change
% 2 == require change
% anything else: no requirement
% columns correspond to testfunctions
% rows correspond to singularities BPC, PD, LPC, and NS respectively
function [S,L] = singularity_matrix

  
  sign_change    = Constants.sign_change;
  sign_constant  = Constants.sign_constant;
  value_change   = Constants.value_change;
  o              = Constants.ignore; 

  S = [ 
  % BPC_testfunc PD_testfunc  LPC_testfunc   NS_testfunc
    sign_change  o            sign_constant  o            % BPC    
    o            sign_change  o              o            % PD   
    o            o            sign_change    o            % LPC   
    o            o            o              value_change % NS
  ];
  L = [ 'BPC';'PD '; 'LPC'; 'NS ' ];
  
  assert(size(S,1) == length(L))
end
%-------------------------------------------------------------------------------
% After a singularity is detected and located, the contL calls this function.
% Usually a normal form coefficient (nfc) is computed, via a function call from
% this function. However, for continuation of cycles by single shooting or
% multiple shooting, the computation of normal form coefficients is not
% implemented yet.
function [failed,s] = process_singularity(id,point,s)
  global cds contopts
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
    threshold = contopts.real_v_complex_threshold;
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
  global contopts
  [sing_mat, ~]                 = singularity_matrix();
  use_locators                  = false(size(sing_mat,1),1);
  use_locators(Constants.NS_id) = true;
  contopts = contset(contopts, 'Locators', use_locators);
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
function p_out = locate(id, p1, p2)
  switch id   
    case Constants.NS_id
      p_out = locate_NS(p1, p2, @testfunctions);
    otherwise
      error('No locator defined for singularity %d', id);
  end
end
%-------------------------------------------------------------------------------
