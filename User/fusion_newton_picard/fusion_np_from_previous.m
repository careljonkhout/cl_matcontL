% continuation of cycles in fusion system
close all
run_init_if_needed

N = 25;
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
handles = feval(odefile);
a = -1; b = -0.3; q_inf = -0.72;
parameters = {a;b;q_inf};
title_string_args = {N,a,b,q_inf};

clear global cds
clear global lds

global cds
cds.preferred_basis_size = 6;
cds.poincare_tolerance = 1e-1;
cds.minimum_period = 1;
cds.nphases = 3*(N-1);
cds.nap = 1;
cds.ActiveParams = 3;
cds.ndim = cds.nphases + cds.nap + 1;
cds.P0 = cell2mat(parameters);
cds.options = contset();
cds.options.PartitionMonodromy = true;
cds.nDiscretizationPoints = 400;
cds.symjac = true;
cds.usernorm = [];
cds.probfile = odefile;
cds.dydt_ode = handles{2};
cds.jacobian_ode = handles{3};
cds.ncoo = cds.nphases;
cds.integrator = @ode15s;


a = -1;
b = -0.3;
q_inf = -0.72;
parameters = {a;b;q_inf};
%load('/home/carel/Documents/cl_matcontL/User/fusion_newton_picard/Data/fusion_np_from_previous_19-Mar-2019_18_22_19/point 16.mat')
%load('/home/carel/Documents/cl_matcontL/User/fusion_newton_picard/Data/fusion_np_from_previous_19-Mar-2019_20_35_42/point 61.mat')
load('/home/carel/Documents/cl_matcontL/User/fusion_newton_picard/Data/fusion_np_from_previous_20-Mar-2019_20_03_57/point 10.mat')
parameters{3} = point.x(end);
cds.previous_phases = point.x(1:end-2);
cds.previous_dydt_0 = cds.dydt_ode(0,cds.previous_phases,parameters{:});

%% Continue limit cycle from orbit



opt = contset();
opt = contset(opt, 'InitStepsize',   7.5e-3);
opt = contset(opt, 'MinStepsize',    1e-10);
opt = contset(opt, 'MaxStepsize',    7.5e-3);
opt = contset(opt, 'MaxNewtonIters', 8);
opt = contset(opt, 'MaxCorrIters',   10);
opt = contset(opt, 'MaxTestIters',   10);
opt = contset(opt, 'VarTolerance',   1e-6);
opt = contset(opt, 'FunTolerance',   1e-6);
% we don't want to adapt
% since it is not implemented
opt = contset(opt, 'Adapt',          1000*1000*1000);
opt = contset(opt, 'MaxNumPoints',   500);
opt = contset(opt, 'contL_SmoothingAngle', 10);
opt = contset(opt, 'CheckClosed',    50000);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       false);
opt = contset(opt, 'Singularities',  false);
opt = contset(opt, 'CIS_UsingCIS',   false);
opt = contset(opt, 'NewtonPicard',   true);
opt = contset(opt, 'console_output_level',   5);
opt = contset(opt, 'contL_DiagnosticsLevel', 5);
opt = contset(opt, 'every_point_in_separate_mat_file', true);
opt = contset(opt, 'contL_ParallelComputing', false);
preferredNumworkers = 2;
if opt.contL_ParallelComputing
  pool = gcp('nocreate');
  if isempty(pool)
    parpool(preferredNumworkers);
  elseif pool.NumWorkers ~= preferredNumworkers
    % shutdown existing pool
    delete(pool);
    % create new pool
    parpool(preferredNumworkers);
  end
end
opt = contset(opt, 'PicardTolerance', 1e-6);


tol = 1e-13;

opt = contset(opt, 'orbit_abs_tol'            , tol);
opt = contset(opt, 'orbit_rel_tol'            , tol);
opt = contset(opt, 'MV_abs_tol'               , tol);
opt = contset(opt, 'MV_rel_tol'               , tol);
opt = contset(opt, 'shoot_abs_tol'            , tol);
opt = contset(opt, 'shoot_rel_tol'            , tol);
opt = contset(opt, 'monodromy_map_abs_tol'    , tol);
opt = contset(opt, 'monodromy_map_rel_tol'    , tol);
opt = contset(opt, 'continue_subspace_rel_tol', tol);
opt = contset(opt, 'continue_subspace_abs_tol', tol);
opt = contset(opt, 'jacobian_rel_tol'         , tol);
opt = contset(opt, 'jacobian_abs_tol'         , tol);




[s, datafile] = contL(@single_shooting, point.x,point.v,opt); 
return;

figure
hold on;
coordinate1 = 1;
coordinate2 = 2;
 title(sprintf( ...
    'fusion N:%d a:%.2f n:%.2f q_{inf}:%.2f', ...
     N,a,b,q_inf));
 xlabel('x_1')
 ylabel('y_1')
 [x, v, h, mult] = loadPoint(datafile); % DV: load computed cycles
%load('Data\testbruss_Orb_LC.mat')    % DV: load singular points

close all
figure
hold on
for i=1:opt.MaxNumPoints
  xx = x(1:end,i);
  period                       = xx(end-1);
  phases_0                     = xx(1:end-2);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = xx(end);
  parameters                   = num2cell(parameters);
  [t,y] = compute_cycle(phases_0,period,parameters);
  coord1_vals = y(:,1);
  coord2_vals = y(:,3);
  plot(coord1_vals,coord2_vals,'b');
  title(sprintf(title_format_string, title_format_args{:}));
  xlabel('n_1')
  ylabel('U_1')
end


function [value, isterminal, direction] = returnToPlane(t, x)
  global cds;
  % x and should be a column vector
  value = cds.previous_dydt_0'*(x-cds.previous_phases);
  isterminal = t > cds.minimum_period ...
    && max(abs(x-cds.previous_phases)) < cds.poincare_tolerance;
  direction = -1;
end

function [t,x] = compute_cycle(x, period, parameters)
  global cds
  f =@(t, y) cds.dydt_ode(t, y, parameters{:});
  integration_opt = odeset(...
    'AbsTol',      1e-10,    ...
    'RelTol',      1e-10,    ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  [t,x] = ode15s(f, [0 period], x, integration_opt);
end

