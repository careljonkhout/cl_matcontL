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
x = hopf.data.x(1:end-1);
ode_parameters(active_parameter) = hopf.data.x(end);

ntst = 20;
ncol = 4;
h = 1e-3;


[x0, v0] = ...
  init_H_LC_L(odefile, x, ode_parameters, active_parameter, h, ntst, ncol);

opts_h_lc = contset();
opts_h_lc = contset(opts_h_lc, ...
  ...
  'MaxNumPoints',           70, ...
  'InitStepsize',           0.2, ...
  'MaxStepsize',            0.2, ...
  'contL_SmoothingAngle',   100, ...
  'newtcorrL_use_max_norm', true, ...
  'Singularities',          true, ...
  'enable_nf_ns',           false, ...
  'enable_nf_lpc',          false, ...
  'enable_nf_pd',           false, ...
  'contL_DiagnosticsLevel', 3, ...
  'console_output_level',   3, ...
  'every_point_in_separate_mat_file', true);
figure
hold on
title_format_string = ...
  'Brusselator N:%d  A:%.0f  B:%.1f  Dx:%.3f  Dy:%.3f';
title_format_args = {N,A,B,Dx,Dy};
xlabel('L')
ylabel('period')
title(sprintf(title_format_string, title_format_args{:}));

singularities = contL(@limitcycleL, x0, v0, opts_h_lc, ...
            'callback',           @plot_T_versus_param, ...
            'stopping_condition', @(point) point.x(end) > 2 );
for singularity = singularities
  plot_sing(singularity)
end

function plot_sing(s, varargin)
  L = s.data.parametervalues(2);
  T = s.data.T;
  plot(L,T,'r*')
  text(L,T,s.label,varargin{:})
end



