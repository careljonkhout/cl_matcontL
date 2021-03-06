% continuation of cycles in fusion system using orthogonal collcocation
function fusion_cycles
  N = 25;                     
  odefile = str2func(sprintf('fusion_N_%d', N));
  a = -1;
  b = -0.3;  
  q_inf = -0.72;

  ode_parameters = {a ; b; q_inf};


  initial_continuation_data = init_collocation_find_stable_cycle( ...
    'initial_point',             ones(3*(N-1),1), ...
    'time_to_converge_to_cycle', 150, ... 
    'odefile',                   odefile,  ...
    'ode_parameters',            ode_parameters, ...
    'active_parameter_index',    3, ...
    'lower_bound_period',        1, ...
    'upper_bound_period',        10, ...
    'n_mesh_intervals',            40, ...
    'poincare_tolerance',        1e-1, ...
    'show_plots',                true ...
  );

  options = contset( ...
    'MaxNumPoints',            70, ...
    'InitStepsize',            0.25, ...
    'MinStepsize',             1e-6, ...
    'MaxStepsize',             0.25, ...
    'MaxNewtonIters',          8, ...
    'MaxCorrIters',            10, ...
    'Backward',                false, ...
    'VarTolerance',            1e-7, ...
    'FunTolerance',            1e-7, ...
    'Adapt',                   3, ...
    'MaxTestIters',            40, ...
    'Multipliers',             true, ...
    'console_output_level',    3, ...
    'Singularities',           true, ...
    'MoorePenrose',            false, ...
    'contL_SmoothingAngle',    1, ...
    'enable_nf_pd',            false, ...
    'enable_nf_lpc',           true, ...
    'enable_nf_ns',            false, ...
    'singularity_callback',    @plot_singularity_of_cycles);


  global cds
  cds.singularity_callback = @plot_singularity_of_cycles;

  figure
  hold on
  title('fusion - continuation of cycles from stable cycle');
  xlabel('q_{inf}')
  ylabel('period')

  contL(@limitcycleL, initial_continuation_data, [], options, ...
                     'callback', @plot_T_versus_param);
end