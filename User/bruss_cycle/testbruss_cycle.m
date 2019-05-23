% continuation of cycles in brusselator
% using orthogonal collocation

format long
format compact
clear global
odefile = @brusselator_1d;
% parameters for the system of ODEs:
N = 400;
L = 0.5; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
ode_parameters = [N L A B Dx Dy];

% x0 will be the trivial equilibrium which will be continued to find Hopf points
x0=[ones(N,1)*A; ones(N,1)*B/A];
% the continuation will be with respect to the second parameter:
active_parameter = 2;
% we run the initializer for an equilibrium continuation:
[x0,v0] = init_EP_EP_L(odefile,x0,ode_parameters,active_parameter);


% we set the options for the equilibrium continuation:
opts_ep_ep = contset();
opts_ep_ep = contset(opts_ep_ep, ...
  'MaxNumPoints', 3, ...
  'Singularities', true);

% we run the equlibrium continuation:
singularities = contL(@equilibriumL,x0,v0, opts_ep_ep);

% we print a table of the Hopf points
struct2table(singularities)

% we set the value of the active parameter (L) (the parameter in which we
% continued the equlibrium) to the value of L at the first Hopf point:
global cds
disp(cds.sout)
hopf  = singularities(2);
x = hopf.data.x;
ode_parameters(active_parameter) = hopf.data.x(end);

% h will be the amplitude of the initial cycle
h = 1;

% ntst is the number of mesh intervals
% That is, the cycle will be represented using ntst polynomials of degree ncol
ntst = 20;
% ncol is the number of collocation point in each mesh interval
ncol = 4;

% We run the initializer for continuation of cycles by collocation from a Hopf
% point:
[x0, v0] = ...
  init_H_LC_L(odefile, x, ode_parameters, active_parameter, h, ntst, ncol);

% We specify the option for the cycle continuation.
opts_h_lc = contset();
opts_h_lc = contset(opts_h_lc, ...
  ...
  'MaxNumPoints',           1000, ...
  'InitStepsize',           0.1, ...
  'MaxStepsize',            0.1, ...
  'contL_SmoothingAngle',   100, ...
  'newtcorrL_use_max_norm', true, ...
  'Singularities',          true, ...
  'enable_nf_ns',           false, ...
  'enable_nf_lpc',          false, ...
  'enable_nf_pd',           false, ...
  'contL_DiagnosticsLevel', 3, ...
  'console_output_level',   3, ..., 
  'every_point_in_separate_mat_file', true);
figure
hold on
title_format_string = ...
  'Brusselator N:%d  A:%.0f  B:%.1f  Dx:%.3f  Dy:%.3f';
title_format_args = {N,A,B,Dx,Dy};
xlabel('L')
ylabel('period')
title(sprintf(title_format_string, title_format_args{:}));

% We run the cycle continuation:
singularities = contL(@limitcycleL, x0, v0, opts_h_lc, 'callback', ...
                                          @plot_T_versus_param);

% We plot the singularities from the cycle continuation
for singularity = singularities
  plot_sing(singularity)
end

function plot_sing(s, varargin)
  L = s.data.x(end);
  T = s.data.x(end-1);
  plot(L,T,'r*')
  text(L,T,s.label,varargin{:})
end



