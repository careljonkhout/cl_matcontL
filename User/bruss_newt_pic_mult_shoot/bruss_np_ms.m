% continuation of cycles in brusselator using multiple shooting with Newton
% Picard. For a description of the algorithm see (bibtex citation follows):
% @phdthesis{lust-phd,
%	  author={Lust, Kurt},
%	  title={Numerical bifurcation analysis 
%        of periodic solutions of partial differential equations},
%	  school={K.U.Leuven},
%	  year={1997},
% }
% and continuer/Newton_Picard_Multiple_Shooting.m
format long
run_init_if_needed
% continuation of cycles cycles in brusselator
odefile = @bruss_1d; %@brusselator_N_2;
N=10;
L = 1.1; A = 1; B = 2.2; Dx = 0.008; Dy = 0.004;
parameters = {N; L; A; B; Dx; Dy};%parameters = {L; A; B; Dx; Dy};
clear global cds
clear global lds
global cds
cds.ActiveParams = 2; %4;
cds.preferred_basis_size = 3;
handles = feval(odefile);
title_format_string = ...
  'Brusselator N:%d  L:%.0f  A:%.0f  B:%.1f  Dx:%.3f  Dy:%.3f';
title_format_args = {N; L; A; B; Dx; Dy;};

cds.poincare_tolerance = 1e-2;
cds.minimum_period = 1;
cds.dydt_ode = handles{2};
cds.jacobian_ode = handles{3};

cds.probfile = odefile;
cds.nap = 1;
cds.q_systems_tolerance = 1e-10;
cds.nphases = 2*N;
cds.ndim = cds.nphases + cds.nap + 1;
cds.P0 = cell2mat(parameters);
cds.options = contset();
cds.options.PartitionMonodromy = cds.nphases > 30;
cds.nDiscretizationPoints = 1000;
cds.symjac = false;
cds.usernorm = [];
cds.probfile = odefile;
cds.ncoo = cds.nphases;
cds.nShootingPoints = 4;

    
int_opt = odeset( ...
  'AbsTol',      1e-10,    ...
  'RelTol',      1e-13,    ...
  'Jacobian',     @(t,y) feval(handles{3},t,y,parameters{:}) ...
);


x0 = ones(cds.nphases,1);
dydt = handles{2};
f =@(t, y) dydt(t, y, parameters{:});
% we compute the cycle from which we start the continuation form by integration:
[t1, x1] = ode15s(f, [0 150], x0, int_opt);

draw_plots = false;
if draw_plots
  figure(1) %#ok<UNRCH>
  plot(t1,x1)
  title(sprintf(title_format_string, title_format_args{:}));
  xlabel('t');
  ylabel('x_1,x_2,y_1,y_2');
end



approximate_period = 30;

cds.previous_phases = x1(end,:)';
cds.previous_dydt_0 = f(0,x0);
int_opt = odeset(int_opt, 'Events', @returnToPlane);

sol = ode15s(f, linspace(0,approximate_period,1000), x1(end,:), int_opt); 
period = sol.x(end);
initial_continuation_data = zeros(cds.nphases * cds.nShootingPoints + 2,1);
for i=0:cds.nShootingPoints-1
  indices = (1:cds.nphases) + i * cds.nphases;
  initial_continuation_data(indices) = ...
    deval(sol, i / cds.nShootingPoints * period);
end

int_opt = odeset(int_opt, 'Events', []);

if draw_plots
  figure(2)
  plot(t2,x2)
  title(sprintf(title_format_string, title_format_args{:}));
  xlabel('t')
  ylabel('x_1,x_2,y_1,y_2');
end



if draw_plots
  figure(4)
  plot(x2(:,1),x2(:,N+1))
  title(sprintf(title_format_string, title_format_args{:}));
  xlabel('x_1')
  ylabel('y_N')

  drawnow
end

%% Continue limit cycle from orbit


opt = contset();
opt = contset(opt, 'MaxNumPoints',   1000);
opt = contset(opt, 'InitStepsize',   1e-1);
opt = contset(opt, 'MinStepsize',    1e-6);
opt = contset(opt, 'MaxStepsize',    1e-1);
opt = contset(opt, 'MaxNewtonIters', 3);
opt = contset(opt, 'MaxCorrIters',   10);
opt = contset(opt, 'MaxTestIters',   10);
opt = contset(opt, 'VarTolerance',   1e-6);
opt = contset(opt, 'FunTolerance',   1e-9);
opt = contset(opt, 'NewtonPicardBasisTolerance',   1e-1);
opt = contset(opt, 'contL_SmoothingAngle',   3);
% we don't want to adapt
% since it is not implemented
opt = contset(opt, 'Adapt',          1000*1000*1000);
opt = contset(opt, 'CheckClosed',    1000);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       true);
opt = contset(opt, 'Singularities',  false);
opt = contset(opt, 'CIS_UsingCIS',   false);
opt = contset(opt, 'NewtonPicard',   false);
opt = contset(opt, 'console_output_level',   5);
opt = contset(opt, 'contL_DiagnosticsLevel', 5);
opt = contset(opt, 'MoorePenrose', false);


initial_continuation_data(end-1) = period;
initial_continuation_data(end  ) = cds.P0(cds.ActiveParams);
initial_continuation_tangent_vector = [];
[s, datafile] = contL(@multiple_shooting, ...
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
  direction = 1;
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