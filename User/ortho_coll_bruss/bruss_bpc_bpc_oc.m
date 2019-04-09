% continuation of cycles in brusselator
format long
run_init_if_needed
% continuation of cycles cycles in brusselator
odefile = @brusselator_1d_no_jac; %@brusselator_N_2;
N=16;
%L = 1; A = 1; B = 2.2; Dx = 0.008; Dy = 0.004;
L = 0.8; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
parameters = {N; L; A; B; Dx; Dy};%parameters = {L; A; B; Dx; Dy};

global cds
handles = feval(odefile);
cds.dydt_ode = handles{2};
active_parameter_index = 2;
handles = feval(odefile);

nphases = 2*N;




%% Continue limit cycle from orbit


opt = contset();
opt = contset(opt, 'MaxNumPoints',   200);
opt = contset(opt, 'InitStepsize',   0.2);
opt = contset(opt, 'MinStepsize',    1e-6);
opt = contset(opt, 'MaxStepsize',    0.2);
opt = contset(opt, 'MaxNewtonIters', 3);
opt = contset(opt, 'MaxCorrIters',   4);
opt = contset(opt, 'MaxTestIters',   30);
opt = contset(opt, 'VarTolerance',   1e-6);
opt = contset(opt, 'FunTolerance',   1e-6);
opt = contset(opt, 'contL_SmoothingAngle',   3);
% we don't want to adapt
% since it is not implemented
opt = contset(opt, 'Adapt',          1000*1000*1000);
opt = contset(opt, 'CheckClosed',    1000);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       false);
opt = contset(opt, 'Singularities',  true);
opt = contset(opt, 'CIS_UsingCIS',   false);
opt = contset(opt, 'console_output_level',   3);
opt = contset(opt, 'every_point_in_separate_mat_file', true, ...
                   'enable_nf_ns', false, ...
                   'enable_nf_pd', false, ...
                   'enable_nf_lpc', false, ...
                   'always_save_s', true, ...
                   'every_point_in_separate_mat_file', true);

ntst = 20; % number of mesh intervals
ncol = 4;  % number of collocation points
init_orb_lc_tolerance = 1e-2;

load('s');
x_bpc = cds.sout(3).data.x;
v_bpc = cds.sout(3).data.v;
s_bpc = cds.sout(3);

global lds
lds.ActiveParams = cds.ActiveParams;

[x0_for_contL,~] = init_BPC_LC(odefile, x_bpc, v_bpc, s_bpc, ntst, ncol, nphases, 1e-2);

figure
hold on
title_format_string = ...
  'Brusselator N:%d  A:%.0f  B:%.1f  Dx:%.3f  Dy:%.3f';
title_format_args = {N,A,B,Dx,Dy};
xlabel('L')
ylabel('period')
title(sprintf(title_format_string, title_format_args{:}));

[sout, datafile] = contL(@limitcycleL,x0_for_contL,[],opt, @plot_T_versus_param); %#ok<*ASGLU>

function plot_T_versus_param(currpoint, trialpoint)
  curr_T      = currpoint.x(end-1);
  trial_T     = trialpoint.x(end-1);
  curr_param  = currpoint.x(end);
  trial_param = trialpoint.x(end);
  plot([curr_param trial_param],[curr_T trial_T],'b');
  drawnow
end
% 
% 
% hold on;
% coordinate1 = 1;
% coordinate2 = 2;
%  title(sprintf( ...
%     'brusselator'))
%  ylabel('y_1')
% x = loadPoint(datafile); % DV: load computed cycles
% %load('Data\testbruss_Orb_LC.mat')    % DV: load singular points
% 
% 
% close all
% figure
% hold on
% for i=1:size(x,2)
%   xx = x(1:end,i);
%   period                       = xx(end-1);
%   phases_0                     = xx(1:2*N);
%   parameters                   = cds.P0;
%   parameters(cds.ActiveParams) = xx(end);
%   parameters                   = num2cell(parameters);
%   [t,y] = compute_cycle(phases_0,period,parameters);
%   coord1_vals = y(:,1);
%   coord2_vals = y(:,end);
%   plot(coord1_vals,coord2_vals,'b');
%   title(sprintf(title_format_string, title_format_args{:}));
%   xlabel('')
%   ylabel('')
% end
% 
% 
% 
% figure
% plot(x(end,:),x(end-1,:))
% title(sprintf(title_format_string, title_format_args{:}));
% xlabel('L')
% ylabel('period')
% 
% return
% figure
% hold on;
% nMults = size(mult,1);
% for i=nMults-10:nMults
%   plot(x(end,:),mult(i, :))
% end
% 
% 
% function [time_values, trajectory] ...
%   = compute_cycle(initial_condition, period, parameters)
%   global cds
%   f =@(t, y) cds.dydt_ode(t, y, parameters{:});
%   integration_opt = odeset(...
%     'AbsTol',      1e-10,    ...
%     'RelTol',      1e-10,    ...
%     'BDF',         'off',   ...
%     'MaxOrder',     5,      ...
%     'NormControl',  'off',  ...
%     'Refine',       1,      ...
%     'Events',      [] ...
%   );
%   time_values = linspace(0,period, 1000);
%   [time_values, trajectory] ...
%     = ode15s(f, time_values, initial_condition, integration_opt);
% end