% Continuation of cycles in the Brusselator 
% using multiple shooting with Newton Picard

% note: the current cl_matcontL implementation (june 2019) of multiple shooting
% with (or without) Newton-Picard is quite slow compared to orthogonal
% collocation for systems of ODEs of up to 350 equations. This file is just
% intended to show how to set up a multiple shooting continuation.

% note: The cycle we will continue here is relatively stable. Therefore, we
% could use single shooting instead of multiple shooting. Single shooting will
% probably always be slightly faster than multiple shooting (regardless of the
% size of the system of ODEs). This file is just intended to show how to set up
% a multiple shooting continuation.

% Note that we define the demo as a function to facilitate testing of the demo
% during the development of cl_matcontL. (Defining a demo as a function prevents
% demo's from being accidentally dependent on some variable in the workspace.)
% It is perfectly fine to run continuations in cl_matcontL using Matlab scripts
% without defining functions.
function cycle_from_hopf_ms
  % parameters for the system of ODEs:
  N = 15;
  L = 0.5; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
  ode_parameters = [N L A B Dx Dy];

  % x0 will be the trivial equilibrium which will be continued to find Hopf points
  % the continuation will be with respect to the second parameter:
  active_parameter = 2;
  % we run the initializer for an equilibrium continuation:
  problem_file = @Brusselator_1d.homogeneous_x0;
  [x0,v0]      = init_EP_EP_L(problem_file,[],ode_parameters,active_parameter);


  % we set the options for the equilibrium continuation:
  opts_ep_ep = contset();
  opts_ep_ep = contset(opts_ep_ep, ...
    'MaxNumPoints', 23, ...
    'Singularities', true);

  % we run the equlibrium continuation:
  singularities = contL(@equilibriumL,x0,v0, opts_ep_ep);

  struct2table(singularities)
  % we set the value of the active parameter (L) (the parameter in which we
  % continued the equlibrium) to the value of L at the first Hopf point:
  hopf  = singularities(2);
  x = hopf.data.x;
  ode_parameters(active_parameter) = hopf.data.x(end);

  % h will be the amplitude of the initial cycle
  h = 1e-1;
  % subspace size is the size of the subspace used in the Newton-Picard algorithm
  subspace_size = floor(N/2);
  nMeshIntervals = 2;

  % we run the initializer for continuation of cycles using single shooting:
  [x0, v0] = init_multiple_shooting_from_hopf(problem_file, x, ...
              ode_parameters, active_parameter, h, nMeshIntervals, subspace_size);

  % we specify the options for the continuation of cycles using single shooting:
  opts_h_lc = contset();
  opts_h_lc = contset(opts_h_lc, ...
    ...
    'MaxNumPoints',           30, ...
    'InitStepsize',           0.1, ...
    'MaxStepsize',            0.1, ...
    'contL_SmoothingAngle',   100, ...
    'newtcorrL_use_max_norm', true, ...
    'Singularities',          true, ...
    'enable_nf_ns',           false, ...
    'enable_nf_lpc',          false, ...
    'enable_nf_pd',           false, ...
    'contL_DiagnosticsLevel', 5, ...
    'console_output_level',   5, ...
    'integration_abs_tol',    1e-9, ...
    'integration_rel_tol',    1e-9, ...
    'FunTolerance',           1e-5, ...
    'VarTolerance',           1e-3, ...
    'NewtonPicard',           true, ...
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
  title_format_args = num2cell([N,A,B,Dx,Dy]);
  % we set the title of the plot.
  % the title has two lines
  title({'Cycle from Hopf - multiple shooting', ...
         sprintf(title_format_string, title_format_args{:})});

  % we run the cycle continuation:
  contL(@multiple_shooting, x0, v0, opts_h_lc, ...
                                              'callback', @plot_T_versus_param);
end


