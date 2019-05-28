clear global

path_to_this_script = get_path;

load(fullfile(path_to_this_script,'Data','bruss_oc_orb_lc'));

N= 31;
 
init_collocation_stable_cycle( ...
  'time_to_converge_to_cycle',                 150, ...
  'initial_point',                             ones(2*N,1), ...
  'odefile',                                   odefile, ....
  'ode_parameters',                            parameters, ...
  'active_parameter_index',                    2, ...
  'lower_bound_period',                        1, ...
  'upper_bound_period',                        4, ...
  'nMeshIntervals',                            20);

odefile = @brusselator_1d_no_jac;

bpc_data = s(2).data;
x = bpc_data.x;
v = bpc_data.v;
s = s(2);
ntst = 20;
ncol = 4;
nphase = 2*N;
h      = 1e-2;

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
opt = contset(opt,'MaxNumPoints',              1000);

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
opt = contset(opt, 'Multipliers',   1);
opt = contset(opt, 'contL_DiagnosticsLevel', 5, ...
                   'console_output_level',   4, ...
                   'newtcorrL_use_max_norm', 1, ...
                   'enable_nf_ns',           false, ...
                   'enable_nf_lpc',          false, ...
                   'enable_nf_pd',           false);

initial_continuation_state = ...
  init_BPC_LC(odefile, x, v, s, ntst, ncol, nphase, h);


[sout,datafile] = contL(@limitcycleL, initial_continuation_state, [], opt, ...
  @plot_T_versus_param);
s = sout;
