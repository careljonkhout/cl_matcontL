function out = single_shooting
%
% Curve file of cycle continuation with single shooting
%
    out{1}  = @curve_func;
    out{2}  = @defaultprocessor;
    out{3}  = @options;
    out{4}  = @jacobian;
    out{5}  = [];%@hessians;
    out{6}  = @testf;
    out{7}  = [];%@userf;
    out{8}  = @process_singularity;
    out{9}  = @singmat;
    out{10} = [];%@locate;
    out{11} = @init;
    out{12} = [];%@done;
    out{13} = @adapt;
    out{14} = @curve_CIS_first_point;
    out{15} = @curve_CIS_step;
end

function func = curve_func(varargin)
  global cds
  x = varargin{1};
  active_par_val               = x(end);
  period                       = x(end-1);
  phases_0                     = x(1:end-2);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = active_par_val;
  parameters                   = num2cell(parameters);
  phases_end                   = shoot(phases_0, period, parameters);
  func = [phases_end - phases_0; 
          (phases_0 - cds.previous_phases)' * cds.previous_dydt_0  ]; 
end

function x_end = shoot(x, period, parameters)
  global cds contopts
  f =@(t, y) cds.dydt_ode(t, y, parameters{:});
  integration_opt = odeset(...
    'AbsTol',      contopts.integration_abs_tol,    ...
    'RelTol',      contopts.integration_rel_tol,    ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  [~, trajectory] = cds.integrator(f, [0 period], x, integration_opt);
  x_end = trajectory(end,:)';
end

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
  d_phi_d_T    = cds.dydt_ode(0,y_end,parameters{:});
  d_s_d_T      = cds.previous_dydt_0' * d_phi_d_T;
  jacobian     = [jacobian [d_phi_d_T; d_s_d_T]];
  d_phi__d_p   = compute_d_phi_d_p(phases, period, parameters);
  d_s__d_p     = cds.previous_dydt_0' * d_phi__d_p;
  jacobian     = [jacobian [d_phi__d_p; d_s__d_p]];
end

function dphidp = compute_d_phi_d_p(x, period, parameters)
  global cds
  ap = cds.ActiveParams;
  h = 1e-6;
  parameters{ap} = parameters{ap} - h;
  phi_1 = shoot(x, period, parameters);
  parameters{ap} = parameters{ap} + 2*h;
  phi_2 = shoot(x, period, parameters);
  dphidp = (phi_2 - phi_1)/h/2;
end
  
function [y_end, monodromy] = compute_monodromy(x, period, parameters)
  global cds
  if ~ cds.options.PartitionMonodromy
    [y_end, monodromy] = monodromy_full(x, period, parameters);
  else
    [y_end, monodromy] = monodromy_column_by_column(x, period, parameters);
  end
end 
    
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

function init(~,~); end

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
    point.tvals = testf(1:8,point.x,point.v,[]);
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
  
function update_multipliers_if_needed(x)
  global cds
  if ~ isfield(cds,'multipliersX') || all(cds.multipliersX ~= x)
    cds.multipliersX = x;
    cds.multipliers = NewtonPicard.SingleShooting.compute_multipliers(x, ...
                       cds.preferred_basis_size);
  end
end


% test functions are used for detecting AND location singularities by bisection
% when detecting ids will be cds.ActTest, and when locating ids will contain
% only those ids relevant to the bifurcation that is being located.
function [out, failed] = testf(ids, x0, v, ~) 
  % unused arguments are v and CISdata
  global cds
  
  failed = false;
  
  PD_id  = 6; % id for period doubling test function
  LPC_id = 7; % id for limit point of cycles test function
  NS_id  = 8; % id for Neimarck-Sacker test function
  
  if any(ismember([PD_id NS_id],ids))
    update_multipliers_if_needed(x0)
  end
  if any(ismember(1:5,ids))
    % detection of Branching points of cycles is currently not implemented
    out(1:5) = ones(5,1);
  end
  if ismember(PD_id, ids)
    % real is needed to ensure that small complex parts induced by roundoff
    % errors do not make the result complex.
    % A complex valued test function would cause false positives when detecting 
    % bifurcations.
    out(6) = real(prod(cds.multipliers + ones(size(cds.multipliers))));
  end
  if ismember(LPC_id, ids)
    out(7) = v(end);
  end
  if ismember(NS_id, ids)
    mults = cds.multipliers;
    psi_ns = 1;
    for i = 1 : (length(mults) - 1)
      for j = (i + 1) : length(mults)
        psi_ns = psi_ns * (real(mults(i) * mults(j)) - 1);
      end
    end
    out(8) = psi_ns;
  end
end

function [S,L] = singmat
  % defines which changes in testfunctions correspond to which singularity type
  % 0 == require sign-change
  % 1 == require sign-non-change
  % 2 == require change
  % anything else: no requirement
  % columns correspond to testfunctions
  % rows correspond to singularities BPC, PD, LPC, and NS respectively
  S = [ 0 0 0 0 8 8 8 8
        8 8 8 8 8 0 8 8
        8 8 8 8 8 8 0 8
        8 8 8 8 0 8 1 0];


  L = [ 'BPC';'PD '; 'LPC'; 'NS ' ];
end
  
function [failed,s] = process_singularity(id,point,s)
  global cds
  x = point.x;
  switch id
  case 1
    % note: as of march 2019 detection for BPC not yet implemented 
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
    smallest_sum = Inf;
    for jk=1:cds.nphases-1
      [val,idx] = min(abs(d(jk+1:cds.nphases)*d(jk)-1));
      if val < smallest_sum
        idx2 = jk+idx;
        smallest_sum = val;
      end
    end
    singularity_is_neutral_saddle = imag(d(idx2)) == 0;
    if singularity_is_neutral_saddle
      s.msg = 'Neutral saddle cycle';
      format_string = 'Neutral Saddle Cycle (period = %e, parameter = %e)\n';
      % A neutral saddle is not really a bifurcation, therefore we use priority 
      % 1 instead of 0, so that it is only logged if 
      % contopts.contL_DiagnosticsLevel is set higher than the default value
      % which is zero.
      print_diag(1, format_string, x(end-1), x(end));
    else
      s.msg = 'Neimark Sacker';
      format_string = 'Neimark-Sacker (period = %e, parameter = %e)\n';
      print_diag(0, format_string, x(end-1) ,x(end));
      %print_diag(0, 'Normal form coefficient = %d\n', s.data.nscoefficient);
    end
  end
  failed = 0;
end  

function options; end

function CISdata = curve_CIS_first_point(x) %#ok<INUSD>
  CISdata = 1;
end
  
function CISdata = curve_CIS_step(x, CISdata_in) %#ok<INUSD>
  CISdata = 1;
end

function [has_changed,x2,v2,CISdata] = adapt(x,v,~,~)
  x2 = x;
  v2 = v;
  CISdata = [];
  has_changed = false;
  % unused arguments are v, CISdata, and tfUpdate

end
