N=248;
L = 0.5; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
ode_parameters = {N L A B Dx Dy};

load(fullfile(get_path(), 'Data', 'most_recent_point.mat'));


init_single_shooting_extend_curve( ...
  'initial_continuation_state',             point.x, ...
  'initial_continuation_tangent_vector',    point.v, ... 
  'ode_parameters',                         ode_parameters, ...
  'active_parameter_index',                 2, ...
  'odefile',                                @brusselator_1d,  ...
  'subspace_size',                          length(point.multipliers) ...
);

opts = contset();
opts = contset(opts, ...
  ...
  'InitStepsize',           0.05, ...
  'MaxStepsize',            0.05, ...
  'MaxNumPoints',           10, ...
  'contL_DiagnosticsLevel', 4, ...
  'console_output_level',   4, ...
  'NewtonPicard',           true, ...
  'contL_SmoothingAngle',   100, ...
  'newtcorrL_use_max_norm', true, ...
  'Singularities',          true, ...
  'enable_nf_ns',           false, ...
  'enable_nf_lpc',          false, ...
  'enable_nf_pd',           false, ...
  'every_point_in_separate_mat_file', true, ...
  'Filename',               session_filename);

%title_format_string = ...
%  'Brusselator N:%d  A:%.2f  B:%.2f  Dx:%.3f  Dy:%.3f';
%title_format_args = {N,A,B,Dx,Dy};


singularities = ...
  contL(@single_shooting, point.x, point.v, opts, ...
            'stopping_condition', @(point) point.x(end) > 20 );


