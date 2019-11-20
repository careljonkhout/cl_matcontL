% continuation of cycles in the Brusselator using orthogonal collocation

% Note that we define the demo as a function to facilitate testing of the demo
% during the development of cl_matcontL. (Defining a demo as a function prevents
% demo's from being accidentally dependent on some variable in the workspace.)
% It is perfectly fine to run continuations in cl_matcontL using Matlab scripts
% without defining functions.
function cycle_from_hopf_oc
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
  h = 0.01;


  % ntst is the number of mesh intervals. That is, the cycle will be represented
  % using a piecewise polynomial function of ntst pieces.
  ntst = 20;
  % ncol is the number of collocation points in each mesh interval
  ncol = 4;

  % We run the initializer for continuation of cycles by collocation from a Hopf
  % point:
  [x0, v0] = init_collocation_from_hopf(...
          @brusselator_1d, x, ode_parameters, active_parameter_index, ...
          h, ntst, ncol);

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
end




