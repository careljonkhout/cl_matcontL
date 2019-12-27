% Continuation of cycles in the Brusselator 
% using single shooting with Newton Picard

% note: the current cl_matcontL implementation (june 2019) of single shooting
% with (or without) Newton-Picard is quite slow compared to orthogonal
% collocation for systems of ODEs of up to 350 equations. This file is just
% intended to show how to set up a single shooting continuation.

% Note that we define the demo as a function to facilitate testing of the demo
% during the development of cl_matcontL. (Defining a demo as a function prevents
% demo's from being accidentally dependent on some variable in the workspace.)
% It is perfectly fine to run continuations in cl_matcontL using Matlab scripts
% without defining functions.
function cycle_from_hopf_ss
  try
    path_to_this_file = get_path();
    my_file = get_latest_singularity_file(path_to_this_file, 'hopf_for_cycle');
    load(my_file, 's');
  catch
    fprintf(['Could not find a file with a Hopf point. ' ...
             'You must run hopf_for_cycle.m first\n']);
    return
  end
  
  singularities        = s;
  hopf                 = singularities(2);
  x                    = hopf.data.x;
  ode_parameters       = hopf.data.P0;
  ode_parameters_cell  = num2cell(ode_parameters);
  [N, ~, A, B, Dx, Dy] = deal(ode_parameters_cell{:});
  % we set the value of the active parameter (L) (the parameter in which we
  % continued the equlibrium) to the value of L at the first Hopf point:
  active_parameter_index                 = 2;
  ode_parameters(active_parameter_index) = hopf.data.x(end);

  % h will be the amplitude of the initial cycle
  h = 1e-1;
  % subspace size is the size of the subspace used in the Newton-Picard algorithm
  subspace_size = 7;

  % we run the initializer for continuation of cycles using single shooting:
  [x0, v0] = init_single_shooting_from_hopf(@brusselator_1d, ...
                   x, ode_parameters, active_parameter_index, h, subspace_size);

  % we specify the options for the continuation of cycles using single shooting:
  opts = contset( ...
    'MaxNumPoints',           30, ...
    'InitStepsize',           0.1, ...
    'MaxStepsize',            0.1, ...
    'contL_SmoothingAngle',   100, ...
    'newtcorrL_use_max_norm', true, ...
    'Singularities',          true, ...
    'enable_nf_ns',           false, ...
    'enable_nf_lpc',          false, ...
    'enable_nf_pd',           false, ...
    'contL_DiagnosticsLevel', 3, ...
    'console_output_level',   3, ...
    'integration_abs_tol',    1e-9, ...
    'integration_rel_tol',    1e-9, ...
    'FunTolerance',           1e-5, ...
    'VarTolerance',           1e-3, ...
    'NewtonPicard',           true, ...
    'parameter_sensitivity_by_finite_diff', true, ...
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
  title({'Cycle from Hopf - single shooting', ...
         sprintf(title_format_string, title_format_args{:})});


  % we run the cycle continuation:
  contL(@single_shooting, x0, v0, opts, 'callback', @plot_T_versus_param);
end




