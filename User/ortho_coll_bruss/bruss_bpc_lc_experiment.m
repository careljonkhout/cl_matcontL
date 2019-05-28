load('s')
handles = limitcycleL;
jac_handle = handles{4};
jac = feval(jac_handle,cds.sout(3).data.x);
[V,~] = eigs(jac(:,1:end-1),1,'smallestabs');
v0 = [V; 0];


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
opt = contset(opt, ...
                   'enable_nf_ns', false, ...
                   'enable_nf_pd', false, ...
                   'enable_nf_lpc', false, ...
                   'always_save_s', true, ...
                  );

ntst = 20; % number of mesh intervals
ncol = 4;  % number of collocation points
init_orb_lc_tolerance = 1e-2;

figure
hold on
title_format_string = ...
  'Brusselator N:%d  A:%.0f  B:%.1f  Dx:%.3f  Dy:%.3f';
title_format_args = {N,A,B,Dx,Dy};
xlabel('L')
ylabel('period')
title(sprintf(title_format_string, title_format_args{:}));

[sout, datafile] = contL(@limitcycleL,cds.sout(3).data.x+v0,v0,opt, @plot_T_versus_param); %#ok<*ASGLU>

function plot_T_versus_param(currpoint, trialpoint)
  curr_T      = currpoint.x(end-1);
  trial_T     = trialpoint.x(end-1);
  curr_param  = currpoint.x(end);
  trial_param = trialpoint.x(end);
  plot([curr_param trial_param],[curr_T trial_T],'b');
  drawnow
end