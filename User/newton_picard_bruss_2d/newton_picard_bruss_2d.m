% continuation of cycles in brusselator
format long
run_init_if_needed
% continuation of cycles cycles in brusselator
odefile = @bruss_2d2; %@brusselator_N_2;
N=5;
L = 1; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
parameters = {N; L; A; B; Dx; Dy};%parameters = {L; A; B; Dx; Dy};
clear global cds
clear global lds
global cds
cds.ActiveParams = 2; %4;
handles = feval(odefile);
title_format_string = ...
  'Brusselator 2D N:%d  L:%.2f  A:%.2f  B:%.2f  Dx:%.3f  Dy:%.3f';
title_format_args = {N; L; A; B; Dx; Dy;};
cds.preferred_basis_size = 15;
global contopts;
contopts = contset();
print_diag(0,[title_format_string '\n'], title_format_args{:});


% cds.poincare_tolerance = 1e-1;
% cds.minimum_period = 1;
% cds.dydt_ode = handles{2};
% cds.jacobian_ode = handles{3};
% 
% cds.probfile = odefile;
% cds.nap = 1;
% 
 cds.nphases = 2*N*N;
% cds.ndim = cds.nphases + cds.nap + 1;
% cds.P0 = cell2mat(parameters);
% cds.options = contset();
% cds.options.PartitionMonodromy = cds.nphases > 30;
% cds.nDiscretizationPoints = 1000;
% cds.symjac = false;
% cds.usernorm = [];
% cds.probfile = odefile;
% cds.ncoo = cds.nphases;

    
int_opt = odeset( ...
  'AbsTol',      1e-10,    ...
  'RelTol',      1e-13    ...
);
%'Jacobian',     @(t,y) feval(handles{3},t,y,parameters{:}) ...


x0 = ones(cds.nphases,1);
dydt = handles{2};
f =@(t, y) dydt(t, y, parameters{:});
[t1, x1] = ode15s(f, [0 150], x0, int_opt);

initial_continuation_data = init_single_shooting( ...
  'point_on_limitcycle',     x1(end,:)', ...
  'odefile',                 odefile,  ...
  'ode_parameters',          parameters, ...
  'active_parameter_index',  cds.ActiveParams, ...
  'lower_bound_period',      1, ...
  'upper_bound_period',      30, ...
  'poincare_tolerance',      1e-1, ...
  'subspace_size',           8);




%% Continue limit cycle from orbit


opt = contset();
opt = contset(opt, 'MaxNumPoints',   100);
opt = contset(opt, 'InitStepsize',   1e-1);
opt = contset(opt, 'MinStepsize',    1e-6);
opt = contset(opt, 'MaxStepsize',    1e-1);
opt = contset(opt, 'MaxNewtonIters', 3);
opt = contset(opt, 'MaxCorrIters',   4);
opt = contset(opt, 'MaxTestIters',   10);
opt = contset(opt, 'VarTolerance',   1e-6);
opt = contset(opt, 'FunTolerance',   1e-6);
opt = contset(opt, 'NewtonPicardBasisTolerance',   1e-6);
opt = contset(opt, 'contL_SmoothingAngle',   3);
opt = contset(opt, 'CheckClosed',    1000);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       false);
opt = contset(opt, 'Singularities',  true);
opt = contset(opt, 'CIS_UsingCIS',   false);
opt = contset(opt, 'NewtonPicard',   true);
opt = contset(opt, 'console_output_level',   4);
opt = contset(opt, 'contL_DiagnosticsLevel', 4);


initial_continuation_tangent_vector = [];
[s, datafile] = contL(@single_shooting, ...
  initial_continuation_data, ...
  initial_continuation_tangent_vector, opt); 


figure
hold on;
coordinate1 = 1;
coordinate2 = 2;
 title(sprintf( ...
    'brusselator'))
 ylabel('y_1')
x = loadPoint(datafile); % DV: load computed cycles
%load('Data\testbruss_Orb_LC.mat')    % DV: load singular points


close all
figure
hold on
for i=1:size(x,2)
  xx = x(1:end,i);
  period                       = xx(end-1);
  phases_0                     = xx(1:end-2);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = xx(end);
  parameters                   = num2cell(parameters);
  [t,y] = compute_cycle(phases_0,period,parameters);
  coord1_vals = y(:,1);
  coord2_vals = y(:,end);
  plot(coord1_vals,coord2_vals,'b');
  title(sprintf(title_format_string, title_format_args{:}));
  xlabel('n_1')
  ylabel('U_1')
end

return
figure
plot(x(end,:),x(end-1,:))
title(sprintf(title_format_string, title_format_args{:}));
xlabel('TODO: set xlabel')
ylabel('period')

figure
hold on;
nMults = size(mult,1);
for i=nMults-10:nMults
  plot(x(end,:),mult(i, :))
end

function [value, isterminal, direction] = returnToPlane(t, x)
  global cds;
  % x and should be a column vector
  value = cds.previous_dydt_0'*(x-cds.previous_phases);
  isterminal = t > cds.minimum_period ...
    && max(abs(x-cds.previous_phases)) < cds.poincare_tolerance;
  direction = -1;
end
function [time_values, trajectory] ...
  = compute_cycle(initial_condition, period, parameters)
  global cds
  f =@(t, y) cds.dydt_ode(t, y, parameters{:});
  integration_opt = odeset(...
    'AbsTol',      1e-10,    ...
    'RelTol',      1e-10,    ...
    'BDF',         'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1,      ...
    'Events',      [],...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  time_values = linspace(0,period, 1000);
  [time_values, trajectory] ...
    = ode15s(f, time_values, initial_condition, integration_opt);
end