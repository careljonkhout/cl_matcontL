function opt = contset(varargin)

  if nargin == 0
    opt = defaultOptions();
    return
  end
  if isstruct(varargin{1})
    if is_options_struct(varargin{1})
      opt = varargin{1};
      i = 2;
    else
      error(['The function conset received an invalid options structure as ' ...
             'it''s first argument'])
    end
  else
    opt = defaultOptions();
    i = 1;
  end

  while i < nargin
    name = varargin{i};
    if ~ischar(name)
      error('Argument %d must be a string.',i);
    end

    if ~isfield(defaultOptions(), name)
      if isUnsupportedOption(name)
        warning(['Option "', name, '" currently not supported by cl_matcontL'])
      else
        error(['Unrecognized option ' name])
      end
    else
      % todo: give error on missing value
      value = varargin{i+1};
      opt.(name) = value;
    end
    i = i + 2;
  end
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
  options.max_rel_funcnorm_increase = 20;    % maximum relative increase in function norm that does not lead to corrections being aborted
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
  options.PRC             =    0;       % PRC means phase response curve. is related to limitcycleL
  options.dPRC            =    0;       % DV 2018
  options.Input           =    0;       % DV 2018
  options.NewtonPicard    =    false;
  options.monodromy_by_finite_differences = false;
  % in case an analytic Jacobian of the system of ODEs is not available, a
  % continuation by single shooting with Newton-Picard can still be done by
  % enabling this option.

  % ActiveUParams
  % ActiveSParams
  % ActiveSParam
  % TSearchOrder
  % note that cl_matcontL allows to set more options
  options.contL_DiagnosticsLevel  =    0;   % Diagnostic Level, -inf to suppress all output
  options.contL_LogFile           =    1;   % Whether a logfile will be created
  options.contL_SmoothingAngle    = pi/50;  % Minimum allowed change in angle between points
  options.contL_ParallelComputing = 0;      % Whether or not to use the parallel computing toolbox
  % note: while trying the parallel computing toolbox matlab 2019a to
  % parellelize part of the comptutations for mutliple shooting with
  % Newton-Picard, I (Carel Jonkhout) found that memory use increased about
  % tenfold, and matlab workers repeatedly abborted. Hence, I do not expect the
  % current version of the parallel computing toolbox to be of much use, unless
  % the computations are extremely parallelizable and really simple, or if one
  % has an enormous amount of RAM available.


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
  options.lsqminnorm_threshold      = [];
  % if lsqminnorm_threshold is set, the linear_solver_collocation will use
  % lsqminnorm instead of "\" to solve linear systems, ignoring singular values
  % less than threshold in the matrix whose inverse action is being computed.
  % This regularizes the solution, and may cure convergence problems.
  options.Multipliers               = false; % enables the computation of multipliers of a cycle
  options.Workspace                 = 1;             % MP?
  options.enable_nf_lpc   =    true;    
  % Set enable_nf_lpc to false to disable computation of the normal form 
  % for limit points of cycle. Normal form computations can take a long time to
  % run. However, nf_lpc does not take as long as nf_pd and nf_ns.
  % Note: nf_lpc is only implemented for limitcycleL.m i.e. cycle continuation by
  % collocation
  options.enable_nf_pd    =    false;
  % Set enable_nf_pd to false to disable computation of the normal form 
  % for period doubling point of cycle. Normal form computations on large systems
  % take a long time to run.
  % Note: nf_pd is only implemented for limitcycleL.mi.e. cycle continuation by
  % collocation
  options.enable_nf_ns    =    false;
  % Set enable_nf_ns to false to disable computation of the normal form 
  % for Neimark Sacker. Normal form computations on large systems
  % take a long time to run.
  % Note: nf_ns is only implemented for limitcycleL.m i.e. cycle continuation by
  % collocation
  options.console_output_level             = 0; % set to 6 to see all debug info.
  
  options.newtcorrL_use_max_norm           = false; 
  % for limitcycleL, this switches the norm used to normalize the tangent vector
  % for the Euclidian norm to the max-norm (a.k.a. infinity norm)
  
  options.MaxPicardIterations              = 15;
  % For Newton-Picard methods. The maximum number of iterations in the Q-system
  % solver, before the solver aborts.
  
  options.PicardTolerance                  = 1e-8;
  % Tolerances for time integration in single shooting, multiple shooting and
  % continuation of single shooting and multiple shooting with Newton Picard
  % Integration tolerances smaller than 1e-13 may give 'unable to meet tolerances'
  % warnings.
  
  options.num_cores                        = feature('numcores');
  % set the number of cores to be used in parallel computing to the number of
  % physical cores on the local machine. Only relevant if
  % contL_ParallelComputing is true.
  
  
  options.integration_abs_tol               = 1e-9;
  options.integration_rel_tol               = 1e-9;
  % the absolute and relative tolerance used in time-integration in the
  % Newton-Picard methods (and some initializers).
  % note that 
  
  options.multipliers_abs_tol               = 1e-9;
  options.multipliers_rel_tol               = 1e-9;
  % the absolute and relative tolerance used in time-integration in the
  % Newton-Picard methods when computing multipliers.
  
  
  options.singularity_callback              = [];
  % One can specify a function that is to by called when a singularity is saved.
  % for instance to plot a singularity right after is it located. For examples
  % see plot_singularity of cycles.


  %% Newton-Picard options
  
  options.basis_grow_threshold              = 7e-1;
  
  % If the monodromy matrix has an eigenvector (that is not already in the
  % basis) with eigenvalue with absolute value larger than basis_grow_threhold,
  % the basis will be expanded. basis_grow_thresold must be stricly less than
  % one.
  
  options.basis_shrink_threshold            = options.basis_grow_threshold/1.4;
  % If there is a eigenvector in the basis with norm less than
  % basis_shrink threshold, the basis will be shrunk. basis_shrink_thresold must
  % be strictly less than basis_grow_threshold.
  
  options.minimum_basis_size                = 3;
  % the minimum_basis_size. Currently only affects single shooting
  
  options.multiplier_print_threshold = 0.8;
  % if a multiplier is larger than multiplier_print_threshold, it is printed, if
  % console_output_level is set high enough, and logged, if
  % contL_DiagnosticsLevel is set high enough
  
  
  %% Continuer options
  
  options.pause                      = false;
  % enable pausing during continuations. Useful for debugging fast running
  % continuations.
  
  options.nsteps_before_pause        = 10;
  % number of continuation steps before a pause, if pausing is enables, see
  % options.pause.
  
  %% Continuer options useful for extending continuations
  
  options.initial_point_index        = 1;
  % the initial value of cds.i. This allows an extended continuation to produce
  % point files that extend the sequence of point files.
  
  options.set_direction              = true;
  % the direction of the continuation is set in contL.m to make the continuation
  % go in the direction corresponding to an increase in the active parameters
  % (if Backward is false) or decrease (if Backward is true) ( if the change in
  % the active parameter is above a small threshold). Setting set_direction to
  % false will disable this "setting of the direction", and go in the direction
  % of the tangent vector that was passed to contL. When extending
  % continuations, this is needed, since otherwise the extended continuation
  % might reverse direction ( for instance, if the first part of the
  % continuation has an uneven number of limit points)
  
  options.is_extension               = false;
  % causes previously saved singularities to by loaded, if the options.Filename
  % is set to the Filename of a previously saved continuation.

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
end

function is_opt_strct = is_options_struct(s)
  default_options = contset();
  default_names = fieldnames(default_options);
  names = fieldnames(s);
  is_opt_strct = true;
  for i=1:length(names)
    is_opt_strct = is_opt_strct && isequal(names{i}, default_names{i});
  end
end