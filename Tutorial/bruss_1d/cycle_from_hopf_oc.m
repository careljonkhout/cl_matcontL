% continuation of cycles in the Brusselator using orthogonal collocation

problem_file = @Brusselator_1d.homogeneous_x0;
% N is the number of mesh points of the discretization of the system of PDE's
N = 31; 
% other parameters of the system of ODEs:
L = 0.5; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
ode_parameters = [N L A B Dx Dy];


% the continuation will be with respect to the second parameter:
active_parameter = 2;
% we run the initializer for an equilibrium continuation:
[x0,v0] = init_EP_EP_L(problem_file,[],ode_parameters,active_parameter);


% we set the options for the equilibrium continuation:

opts_ep_ep = contset(...
  'MaxNumPoints', 23, ...
  'Singularities', true);

% we run the equlibrium continuation:
singularities = contL(@equilibriumL,x0,v0, opts_ep_ep);

% we set the value of the active parameter (L) (the parameter in which we
% continued the equlibrium) to the value of L at the first Hopf point:
hopf  = singularities(2);
x = hopf.data.x;
ode_parameters(active_parameter) = hopf.data.x(end);

% h will be the amplitude of the initial cycle
h = 0.01;


% ntst is the number of mesh intervals. That is, the cycle will be represented
% using a piecewise polynomial function of ntst pieces.
ntst = 20;
% ncol is the number of collocation points in each mesh interval
ncol = 4;

% We run the initializer for continuation of cycles by collocation from a Hopf
% point:
[x0, v0] = init_collocation_from_hopf(...
            problem_file, x, ode_parameters, active_parameter, h, ntst, ncol);

% We specify the options for the cycle continuation.
opts_h_lc = contset( ...
  ...
  'MaxNumPoints',           70, ...
  'InitStepsize',           0.1, ...
  'MaxStepsize',            0.1, ...
  'contL_SmoothingAngle',   100, ...
  'newtcorrL_use_max_norm', true, ...
  'Singularities',          true, ...
  'enable_nf_ns',           false, ...
  'enable_nf_lpc',          false, ...
  'enable_nf_pd',           false, ...
  'contL_DiagnosticsLevel', 0, ...
  'console_output_level',   0, ...
  'singularity_callback',   @plot_singularity_of_cycles);

% we open a plot window:
figure
% each line segment of the approximation of the curve is plotted separately
% therefore, the plot has to "hold" the previously plotted line segments
hold on
% we set the axes labels
xlabel('L')
ylabel('period')
title_format_string = 'Brusselator N:%d  A:%.0f  B:%.1f  Dx:%.3f  Dy:%.3f';
title_format_args = {N,A,B,Dx,Dy};
% we set the title of the plot.
% the title has two lines
title({'Cycle from Hopf - orthogonal collocation', ...
       sprintf(title_format_string, title_format_args{:})});


% we run the cycle continuation:
contL(@limitcycleL, x0, v0, opts_h_lc, 'callback', @plot_T_versus_param);





