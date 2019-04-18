
% continuation of cycles in brusselator
odefile = @brusselator_1d_no_jac;
N = 31;
L = 0.9; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
parameters = {N; L; A; B; Dx; Dy};
handles = feval(odefile);

 
initial_continuation_state = init_collocation_stable_cycle( ...
  'time_to_converge_to_cycle',                 150, ...
  'initial_point',                             ones(2*N,1), ...
  'odefile',                                   odefile, ....
  'ode_parameters',                            parameters, ...
  'active_parameter_index',                    2, ...
  'lower_bound_period',                        1, ...
  'upper_bound_period',                        4, ...
  'nMeshIntervals',                            20);

opt = contset();
opt = contset(opt,'contL_LogFile',          true); 
opt = contset(opt,'contL_DiagnosticsLevel',    3);  
opt = contset(opt,'Backward',              false);
opt = contset(opt,'InitStepsize',            0.1);
opt = contset(opt,'MinStepsize',            1e-5);
opt = contset(opt,'MaxStepsize',              0.1);
opt = contset(opt,'MaxCorrIters',             12);
opt = contset(opt,'MaxNewtonIters',            5);
opt = contset(opt,'FunTolerance',           1e-6);
opt = contset(opt,'VarTolerance',           1e-5);
opt = contset(opt,'contL_SmoothingAngle',   pi/30); 
opt = contset(opt,'Singularities',             1);
opt = contset(opt,'Userfunctions',             0);
opt = contset(opt,'MaxNumPoints',              50);

opt = contset(opt,'MaxTestIters',             10);
%opt = contset(opt,'contL_Testf_VarTolerance', 1e-4); 
%opt = contset(opt,'contL_Testf_FunTolerance', 1e-5); 
opt = contset(opt,'CIS_SparseSolvers',         1);
opt = contset(opt,'CIS_NStableRef',            4);
opt = contset(opt,'CIS_MaxUnstable',           5); %new
opt = contset(opt,'CIS_Ric_Cayley_Shift',      1); 
opt = contset(opt,'contL_EQ_BranchingMethod',  2); 
opt = contset(opt, 'CIS_UsingCIS',   false);
    % disable smoothing by angle:
opt = contset(opt, 'contL_SmoothingAngle', pi/2);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename',     'testbruss_Orb_LC');
opt = contset(opt, 'Multipliers',   1);
opt = contset(opt, 'contL_DiagnosticsLevel', 5, ...
                   'console_output_level',   4, ...
                   'newtcorrL_use_max_norm', 1, ...
                   'enable_nf_ns',           false, ...
                   'enable_nf_lpc',          false, ...
                   'enable_nf_pd',           false);

opt.TestTolerance = 1e-5;



[sout,datafile] = contL(@limitcycleL,initial_continuation_state,[],opt, ...
  @plot_T_versus_param);
s = sout;

