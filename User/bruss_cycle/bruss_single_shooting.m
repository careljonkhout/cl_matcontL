% continuation of cycles in brusselator

format long
format compact
clear global
% continuation of cycles cycles in brusselator
odefile = @brusselator_1d;
N=248;
L = 0.5; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
ode_parameters = [N L A B Dx Dy];

x0=[ones(N,1)*A; ones(N,1)*B/A];
active_parameter = 2;

[x0,v0] = init_EP_EP_L(odefile,x0,ode_parameters,active_parameter);

opts_ep_ep = contset();

opts_ep_ep = contset(opts_ep_ep, ...
  'MaxNumPoints', 23, ...
  'Singularities', true);


s = contL(@equilibriumL,x0,v0, opts_ep_ep);

struct2table(s)
hopf  = s(2);
x = hopf.data.x;
ode_parameters(active_parameter) = hopf.data.x(end);

ntst = 20;
ncol = 4;
h = 2e-1;

[x0, v0] = ...
  init_single_shooting_from_hopf( ...
    odefile, x, ode_parameters, active_parameter, h, 7);

opts_h_lc = contset();
opts_h_lc = contset(opts_h_lc, ...
  ...
  'InitStepsize',           0.05, ...
  'MaxStepsize',            0.05, ...
  'MaxNumPoints',           30, ...
  'contL_DiagnosticsLevel', 4, ...
  'console_output_level',   4, ...
  'NewtonPicard',           true, ...
  'contL_SmoothingAngle',   100, ...
  'newtcorrL_use_max_norm', true, ...
  'Singularities',          true, ...
  'enable_nf_ns',           false, ...
  'enable_nf_lpc',          false, ...
  'enable_nf_pd',           false);

hold on
title_format_string = ...
  'Brusselator N:%d  A:%.2f  B:%.2f  Dx:%.3f  Dy:%.3f';
title_format_args = {N,A,B,Dx,Dy};
xlabel('L')
ylabel('period')
title(sprintf(title_format_string, title_format_args{:}));


singularities = ...
  contL(@single_shooting, x0, v0, opts_h_lc, @plot_T_versus_param);



for singularity = singularities
  plot_sing(singularity)
end

function plot_sing(s, varargin)
  parameter_value = s.data.x(end);
  T = s.data.x(end-1);
  plot(parameter_value,T,'r*')
  text(parameter_value,T,s.label,varargin{:})
end



