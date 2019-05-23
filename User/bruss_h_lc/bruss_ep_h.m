% continuation of cycles in brusselator
close all
format long
format compact
clear global
% continuation of cycles cycles in brusselator
odefile = @brusselator_1d_no_jac; %@brusselator_N_2;
N=400;
L = 0.5; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
ode_parameters = [N L A B Dx Dy];

x0=[ones(N,1)*A; ones(N,1)*B/A];
active_parameter = 2;

[x0,v0] = init_EP_EP_L(odefile,x0,ode_parameters,active_parameter);

opts_ep_ep = contset();

opts_ep_ep = contset(opts_ep_ep, ...
  'MaxNumPoints', 4, ...
  'Singularities', true);


s = contL(@equilibriumL,x0,v0, opts_ep_ep);

struct2table(s)
x = s(2).data.x(1:end-1);
ode_parameters(active_parameter) = s(3).data.x(end);

ntst = 20;
ncol = 4;
h = 1e-2;


[x0, v0] = ...
  init_H_LC_L(odefile, x, ode_parameters, active_parameter, h, ntst, ncol);

opts_h_lc = contset();
opts_h_lc = contset(opts_h_lc, ...
  ...
  'MaxStepsize',            0.5, ...
  'contL_SmoothingAngle',   100, ...
  'newtcorrL_use_max_norm', true, ...
  'Singularities',          true, ...
  'enable_nf_ns',           false, ...
  'enable_nf_lpc',          false, ...
  'enable_nf_pd',           false, ...
  'contL_DiagnosticsLevel', 3, ...
  'console_output_level',   3);

hold on
title_format_string = ...
  'Brusselator N:%d  A:%.0f  B:%.1f  Dx:%.3f  Dy:%.3f';
title_format_args = {N,A,B,Dx,Dy};
xlabel('L')
ylabel('period')
title(sprintf(title_format_string, title_format_args{:}));

contL(@limitcycleL, x0, v0, opts_h_lc, 'callback', @plot_T_versus_param);


