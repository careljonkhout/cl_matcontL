function opt = contset(opt, name, value) %#ok<INUSD>

if nargin == 0
    opt = defaultOptions();
    return
end

if ~ischar(name)
    error('2nd argument must be a string');
end

if ~isfield(defaultOptions(), name)
    if isUnsupportedOption(name)
        warning(['Option "', name, '" currently not supported by cl_matcontL'])
    else
        error(['Unrecognized option ' name])
    end
else
  eval(['opt.', name, '= value;']);
end

%----------------------------------------------------
function options = defaultOptions()
                            %% BVP Options
 
 %% fprintf('          AbsTol: [ positive scalar or vector {1e-6} ]\n');
 %%fprintf('          RelTol: [ positive scalar {1e-3} ]\n');  
 %%fprintf('    SingularTerm: [ matrix ]\n'); 
 %%fprintf('       FJacobian: [ function_handle ]\n');
 %%fprintf('      BCJacobian: [ function_handle ]\n');
 %%fprintf('           Stats: [ on | {off} ]\n');
 %%fprintf('            Nmax: [ nonnegative integer {floor(10000/n)} ]\n'); 
 %%fprintf('      Vectorized: [ on | {off} ]\n'); 

  %%options.BVP_AbsTol         = 1e-6;  % Tolerance for point Correction ???
  %%options.BVP_RelTol         = 1e-3;  % Tolerance for point Correction ??? 
  

  
                               %% CIS Options
                               
options.CIS_UsingCIS        =    1;         % Does the problem use CIS                             
options.CIS_SparseSolvers   =    0;         % Whether sparse solvers should be used
options.CIS_DetectOverlap   =    1;         % Whether to detect overlap
options.CIS_AdaptOnOverlap  =    1;         % Whether to adapt space on overlap
options.CIS_MaxAngle        = 1e-1;         % Maximal angle between subspaces
options.CIS_resizeable_flag =    1;         % Whether to resize subspace   MP 2018
options.CIS_NUnstable       =   -1;         % Number of initial unstable eigs (negative means contCIS should compute it). MP 2018 
%options.CIS_InitSubSize    =    4;         % MP 2018 Initial subspace size (NUnstable + NStableRef)                       
options.CIS_NSub            =    4;         % Initial subspace size (if not NUnstable + NStableRef, then error msg)   MP 2018 
options.CIS_NExtra          =    4;         % Number of extra eigs to compute
options.CIS_NStableRef      =    4;         % Required stable reference eigs
%options.CIS_NUnstableGuess  =    6;         % do not need 
options.CIS_MaxUnstable     =   10;         % too big?

                            %% CIS: RIC solver options
options.CIS_Ric_Cayley_Shift   = 1;        % Cayley transform shift parameter (sparse solvers only)
options.CIS_Ric_Euler          = 0;        % use Euler predictor?
options.CIS_Ric_FunTolerance   = 1e-5;     % Newton residual convergence tolerance
options.CIS_Ric_VarTolerance   = 1e-5;     % Newton step convergence tolerance
options.CIS_Ric_MaxNewtonIters = 15;       % maximum allowed Newton steps
options.CIS_Ric_PartialQ       = 0;        % Keep only Q1 even in dense case?
options.CIS_Ric_SubspaceSelect = 'ric';    % Method used to select subspace (or 'eig')    MP 2018
options.CIS_Ric_Transform      = 'cayley'; % transformation (none, invert, cayley) (sparse solvers only)
options.CIS_Ric_schur_blsz     = 3;        % Block size used for Schur decomposition (sparse solvers only)
options.CIS_Ric_schur_nbls     = 10;       % Number of blocks used in computation (sparse solvers only)
options.CIS_Ric_schur_maxit    = 100;      % Maximum number of iterations for Schur decomposition (sparse solvers only)
options.CIS_Ric_schur_tol      = 1e-6;     % Tolerance use in Schur decomposition (sparse solvers only)

                            %% Continuer Options
% DV: Continuer options mostly match standard matcont options (compare to list of names in contid.m)
% A few standard matcont options are currently not supported, and
% cl_matcontL has additional options
% A warning will be printed when an unsupported options is set
options.InitStepsize   = 0.01;        % Initial Stepsize    DV: old name 'Cont_InitStepsize'
options.MinStepsize    = 1e-5;        % Minimum stepsize    DV: old name 'Cont_MinStepsize'
options.MaxStepsize    =  0.1;        % Maximum stepsize    DV: old name 'Cont_MaxStepsize'

options.h_inc_fac      = 1.3;         % Factor to increase stepsize  MP: old name 'Stepsize_Inc_fac'
options.h_dec_fac      = 0.5;         % Factor to decrease stepsize  MP: old name 'Stepsize_Dec_fac'   
options.MaxCorrIters   =  10;         % Maximum number of correction steps   DV: old name 'Cont_MaxCorrIters'
options.MaxNewtonIters =   3;         % Maximum number of times that the jacobian is computed per continuation step DV: old name 'Cont_MaxNewtonIters'
% MaxTestIters (SUPPORTED: see locator options below)
options.MoorePenrose   =   1;         % Solver to be used   DV: old name 'Cont_Solver' 
options.SymDerivative  = [1 1 1 1 1]; % ADDED Boolen Array indicating whether symbolic derivaties are used (when supplied in problem file)
options.SymDerivativeP = [1 1];       % ADDED Boolen Array indicating whether symbolic derivaties are used (when supplied in problem file)
options.Increment      = 1e-5;        % Increment for finite difference approximations    DV: old name 'Cont_IncrFinDiff'
options.FunTolerance   = 1e-6;        % residual tolerance for curve   DV: old name     'Cont_FunTolerance'
options.VarTolerance   = 1e-6;        % Tolerance for curve            DV: old name     'Cont_VarTolerance'
options.NewtonPicardBasisTolerance = 1e-6; % Tolerance for the basis in Newton Picard corrections, see continuer/Newton_Picard_Correction.m
options.NewtonPicardMaxSubspaceIterations = 10; % number of iterations before subspace iteration in Newton Picard is aborted, see continuer/Newton_Picard_Correction.m
% TestTolerance (NOT SUPPORTED: cl_matcontL now makes a distinction between
% fun and var tolerance, see locator options below) MP why??????
options.Singularities   =    0;       % Whether or not to detect Singularities     DV: old name 'Cont_Singularities'
options.MaxNumPoints    =  300;       % Maximum Number of Steps to be taken        DV: old name 'Cont_MaxNumPoints'
options.Backward        =    0;       % Direction of continuation {0,1} = {Forward, Backwards}   DV: old name 'Cont_Direction'
options.CheckClosed     =   50;       % Number of steps to start checking for closed loops   DV: old name 'Cont_CheckClosed '
% Testfunctions
% Workspace
options.Locators        =   [];       % Boolean Array: Toggles which locators to use   DV: old name 'Loc_UseLocators'
options.Adapt           =    3;       % Number of steps before automatic adapt         DV: old name 'Cont_AdaptSteps'
% IgnoreSingularity (SUPPORTED: see locator options below)
% ActiveParams
options.Multipliers     =    0;       % Whether or not to calculate multipliers DV 2018
options.Eigenvalues     =    0;       % Whether or not to calculate and save eigenvalues     DV: old name 'Cont_Eigenvalues'
options.Userfunctions   =    0;       % Whether or not to detect userfunctions      DV: old name 'Cont_Userfunctions'
options.UserFuncInfo    =    [];      % Struct:                                        DV: old name 'Cont_UserfuncInfo'
options.PRC             =    0;       % DV 2018
options.dPRC            =    0;       % DV 2018
options.Input           =    0;       % DV 2018
options.NewtonPicard    =    false;

                                                
% ActiveUParams
% ActiveSParams
% ActiveSParam
% TSearchOrder
% note that cl_matcontL allows to set more options
options.contL_DiagnosticsLevel  =    0;   % Diagnostic Level, -inf to suppress all output
options.contL_LogFile           =    1;   % Whether a logfile will be created
options.contL_SmoothingAngle    = pi/50;  % Minimum allowed change in angle between points
options.contL_ParallelComputing = 0;      % Whether or not to use the parallel computing toolbox


                            %% Locator options
options.contL_EQ_BranchingMethod = 0;    % 0: Normal Form Method, 1: Perpendicular Guess, 2: Bisection
options.contL_EQ_SingLinSystTol  = 0.01; % Bordered System Options
options.IgnoreSingularity        = [];   % Boolean Array: Toggles which singularities to ignore  DV: old name 'Loc_IgnoreSings'

options.MaxTestIters     =   20;         % Tolerances for bisection to detect singularities (without locators)   DV: old name 'Loc_Testf_MaxIters'
options.contL_Testf_FunTolerance = 1e-5; % DV: similar to TestTolerance
options.contL_Testf_VarTolerance = 1e-4; % DV: similar to TestTolerance

options.contL_Userf_MaxIters     = 10;   % Tolerances for user function location
options.contL_Userf_FunTolerance = 1e-5;
options.contL_Userf_VarTolerance = 1e-4;

options.contL_Loc_MaxCorrIters       = 5; % Tolerances used in locators
options.contL_Loc_FunTolerance   = 1e-5;
options.contL_Loc_VarTolerance   = 1e-4;

                            %% limitcycle options
                            
options.Multipliers               = 0;             % MP?
options.Workspace                 = 1;             % MP?
options.enable_nf_lpc   =    true;    
% Set enable_nf_lpc to false to disable computation of the normal form 
% for limit points of cycle. Normal form computations can take a long time to
% run.
% Note: nf_lpc is only implemented for limitcycle.m
options.enable_nf_pd    =    true;
% Set enable_nf_pd to false to disable computation of the normal form 
% for period doubling point of cycle. Normal form computations on large systems
% take a long time to run.
% Note: nf_pd is only implemented for limitcycle.m
options.enable_bpc     = true;
options.bpc_tolerance = options.contL_Userf_FunTolerance;
options.console_output_level = 0; % set to 5 to see all debug info.
options.newtcorrL_use_max_norm = false;
options.always_save_s = true;
options.SingularTestFunction = false;
options.every_point_in_separate_mat_file = false;
options.MaxPicardIterations = 15;
options.PicardTolerance     = 1e-10;

                            %% Determine testpath
options.Filename = [];           
ST = dbstack('-completenames');
if length(ST) > 2
    name = ST(3).file;
    % remove .m (or other extension)
    ind = find(name == '.', 1, 'last');
    if ~isempty(ind) && ind > length(name) - 3
        options.TestPath = name(1:ind-1);
    else
        options.TestPath = name;
    end
else
    options.TestPath = []; 
end

function b = isUnsupportedOption(name)
name = deblank(name);   % remove spaces
if strcmp(name, 'TestTolerance') || ...
        strcmp(name, 'TestFunctions') || ...
        ...
        strcmp(name, 'ActiveParams') || ...
         ...
        strcmp(name, 'PRC') || ...
        strcmp(name, 'dPRC') || ...
        strcmp(name, 'Input') || ...
        strcmp(name, 'ActiveUParams') || ...
        strcmp(name, 'ActiveSParams') || ... 
        strcmp(name, 'ActiveSParam') || ...
        strcmp(name, 'TSearchOrder')
    b = 1;
else
    b = 0;
end
% strcmp(name, 'Multipliers') || ...
% strcmp(name, 'Workspace') || ...