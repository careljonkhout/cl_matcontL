function [sout, datafile] = contL(curvefile, x0, v_cont, opts, varargin)
  %
  % CONTINUE(curvefile, x0, v0, options)
  %
  % Continues the curve from x0 with optional directional vector v0
  % options is a option-vector created with CONTSET
  % The first three parameters are mandatory.
  global cds contopts

  %% I. Initialization
  if size(x0, 2) > 1
    error('Initial point must be a column vector')
  elseif isempty(x0)
    error('Initial point must be non-empty')
  end
  if nargin < 4 || isempty(opts)
    contopts = contset();
  else
    contopts = opts;
  end

  additional_arguments.callback           = [];
  additional_arguments.stopping_condition = [];
  if nargin > 4 
    i = 1;
    while i <= length(varargin)
      if ~ ischar(varargin{i})
        error('Please specify additional arguments as key-value pairs')
      end
      if ~ isfield(additional_arguments, varargin{i})
        error([varargin{i} ' is not a valid option.'])
      end
      additional_arguments.(varargin{i}) = varargin{i+1};
      i = i+2;
    end
  end
  if isfield(cds, 'curve') && ~ isequal(cds.curve, curvefile)
    warning(['The field cds.curve and the specified curvefile do no match. ', ...
           'cds.curve is %s and the specified curvefile is %s. ', ...
           'using %s'], ...
           func2str(cds.curve), func2str(curvefile), func2str(curvefile));
    cds.curve = curvefile;

  end
  % curvefile must be loaded before opening data files
  loadCurveFile(curvefile);
  feval(cds.curve_options);

  [datafile, ~]  = openFiles();

  AdaptSteps         = contopts.Adapt;
  CheckClosed        = contopts.CheckClosed;
  MaxNumPoints       = contopts.MaxNumPoints;
  Singularities      = contopts.Singularities;
  SmoothingAngle     = contopts.contL_SmoothingAngle;
  Userfunctions      = contopts.Userfunctions;
  IgnoreSings        = contopts.IgnoreSingularity;
  UseLocators        = contopts.Locators;
  UserInfo           = contopts.UserFuncInfo;
  UsingNewtonPicard  = contopts.NewtonPicard;
  if contopts.contL_ParallelComputing && isempty(gcp('NoCreate'))
    % initialize new parallel pool when none is available
    parpool(contopts.num_cores);   
  end
  
  if contopts.PoincareIterations && Singularities
    Singularities = false;
    contopts.Singularities = false;
    warning(['It is not possible to detect singularities while using ' ...
             'Poincare Iterations. ' ...
             'Singularity detection has been switched off'])
  end
  
  if UsingNewtonPicard && ( ~ isequal(curvefile, @single_shooting   ) ...
                       &&   ~ isequal(curvefile, @multiple_shooting ) ...
                       &&   ~ isequal(curvefile, @perioddoubling_ss))
    warning(['Newton-Picard not implemented for the curve file %s. ' ...
             'Newton-Picard will not be used.'], func2str(curvefile))
    contopts.NewtonPicard = false;
    UsingNewtonPicard = false;
  end

  cds.h          = contopts.InitStepsize;
  cds.h_max      = contopts.MaxStepsize;
  cds.h_min      = contopts.MinStepsize;
  cds.h_inc_fac  = contopts.h_inc_fac;  
  cds.h_dec_fac  = contopts.h_dec_fac; 
  cds.tfUpdate       = 0;
  cds.i              = contopts.initial_point_index;
  cds.lastpointfound = 0;

  %% determine active singularities and testfunctions
  cds.nActTest = 0;
  if Singularities
    [cds.S , cds.SingLables] = feval(cds.curve_singmat);
    nSing                    = size(cds.S,1);
    cds.S(IgnoreSings,:)     = Constants.ignore;
    cds.ActSing              = setdiff(1:nSing, IgnoreSings);
    cds.nActSing             = length(cds.ActSing);
    cds.ActTest              = find( sum((cds.S ~= Constants.ignore),1) > 0 );
  end

  %% userfunctions
  if Userfunctions
    if length(contopts.UserFuncInfo) == length(cds.userf)
       cds.nUserf = length(cds.userf);     % DV: Does this work properly??
    else
       error('Wrong user info specified');
    end
  end

  cds.StartTime = clock;
  cds.first_point = true;

  feval(cds.curve_init, x0, v_cont);
  cds.newtcorrL_needs_CISdata = 0;

  if isequal(curvefile, @limitcycleL) || ...
    isequal(curvefile, @limitpointcycle)
     if isempty(v_cont)
       v_cont = find_initial_tangent_vector(x0, v_cont, 1);
       print_diag(1,'found tangent vector\n')
     end
     first_point = newtcorrL(x0, v_cont, 1);
  elseif UsingNewtonPicard
    if isempty(v_cont)
      v_cont = NP_find_tangent_vector(curvefile, x0);
    end
    first_point = NP_do_corrections(x0,v_cont);
    if isempty(first_point)
      print_diag(0,'Correction of first point does not converge\n');
      return;
    end
    first_point.v =  NP_find_tangent_vector(curvefile, x0, v_cont);
  elseif contopts.PoincareIterations
    first_point = poincare_iterations(x0);
    if ~ isempty(first_point)
      first_point.v = [];
    end
  else
    try 
      cds.curve_func    (x0); 
      cds.curve_jacobian(x0);
    catch
      cds.newtcorrL_needs_CISdata = 1;
    end
    CISdata0 = [];
    if isempty(v_cont)
      if cds.newtcorrL_needs_CISdata
        CISdata0 = feval(cds.curve_CIS_first_point, x0);
        if isempty(CISdata0)
          print_diag(0,'contL: failed to intialize CIS algorithm.\n');
          sout = [];
          return;
        end
      end
      v_cont = find_initial_tangent_vector(x0, v_cont, CISdata0);
      if isempty(v_cont)
        print_diag(0,'contL: failed to find initial tangent vector.\n');
        sout = [];
        return;
      end
    end
    first_point = newtcorrL(x0, v_cont, CISdata0);
  end

  if isempty(first_point)
    print_diag(0,'contL: no convergence at x0.\n');
    sout = [];
    return;
  end
  first_point.h = cds.h;

  %% CIS data algorithm
  first_point.CISdata = feval(cds.curve_CIS_first_point, x0);
  if isempty(first_point.CISdata)
    print_diag(0,'contL: failed to intialize CIS algorithm.\n'); 
    sout = []; 
    return; 
  end


  %% Direction Vector Determination
  if contopts.set_direction && ~ contopts.PoincareIterations
    Backward       = contopts.Backward;
    if abs(first_point.v(end)) < 1e-6
      Vdir = sign(sum(first_point.v(1:end-1)));
    else
      Vdir = sign(first_point.v(end));
    end
    if (Backward && Vdir > 0) || (~Backward && Vdir < 0)
        first_point.v = -first_point.v;
    end
  end

  %% Test and user functions

  first_point.tvals = [];
  if Singularities
    % We calculate all testfunctions at once
    [first_point.tvals, failed] = EvalTestFunc(0, first_point);
    if failed
      message = 'contL: Evaluation of test functions failed at start point.\n';
      print_diag(0, message);
      sout     = [];
      datafile = [];
      return
    end
  end

  first_point.uvals = [];
  if Userfunctions
    [first_point.uvals,failed] = feval(cds.curve_userf, 0, UserInfo, x0);
    if failed
      message = 'contL: Evaluation of user functions failed at start point.\n';
      print_diag(0, message);
      sout     = [];
      datafile = [];
      return
    end
  end

  first_point.angle = 0;
  first_point = DefaultProcessor(first_point); 
  cds.first_point = false;
  %% II. Main Loop
  currpoint = first_point;

  while cds.i < MaxNumPoints && ~cds.lastpointfound
    corrections = 1;
    print_diag(1,'\n --- Step %d ---\n',cds.i);
    step_start_time = tic;
    while true

      %% A. Predict
      if UsingNewtonPicard
        xpre = currpoint.x + cds.h * currpoint.v;
      elseif contopts.PoincareIterations
        xpre = currpoint.x;
        if contopts.Backward
          xpre(end) = xpre(end) - cds.h;
        else
          xpre(end) = xpre(end) + cds.h;
        end
      else
        xpre = currpoint.x + cds.h * currpoint.v(1:cds.ndim);
      end
      reduce_stepsize = 0;

      %% B. Correct
      try
        if UsingNewtonPicard
           trailpoint = NP_do_corrections(xpre, currpoint.v);
          if ~ isempty(trailpoint)
            trailpoint.v = NP_find_tangent_vector(...
                                    curvefile, trailpoint.x, trailpoint.v);
          end
        elseif contopts.PoincareIterations
          trailpoint = poincare_iterations(xpre);
          if ~ isempty(trailpoint)
            trailpoint.v = [];
          end
        else 
          trailpoint = newtcorrL(xpre, currpoint.v, currpoint.CISdata);
        end
      catch my_error
        switch my_error.identifier
          case {'cvode:integrator_error', 'npms:wrong_basis_size'}
            print_diag(0, my_error.message);
            trailpoint = [];
          otherwise
            rethrow(my_error)
        end
      end
      
      

      if isempty(trailpoint)
        if UsingNewtonPicard 
          print_diag(1, 'contL: NP_do_corrections failed\n')
        elseif contopts.PoincareIterations
          print_diag(1, 'contL: poincare_iterations failed\n')
        else
          print_diag(1, 'contL: newtcorrL failed\n')
        end
        reduce_stepsize = 1;
      end

      %% C. curve smoothing
      if ~reduce_stepsize
        trailpoint.h = cds.h;
        if ~ isempty(currpoint.v)
          trailpoint.angle = innerangle(currpoint.v, trailpoint.v);
        end
        if isfield(trailpoint, 'angle') && trailpoint.angle > SmoothingAngle ... 
                && cds.i > 1
          print_diag(1, ...
            'contL: Innerangle too large, innerangle is %f degrees\n', ...
            trailpoint.angle / pi * 180);
          reduce_stepsize = 1;
        end

        % In single shooting, the corrector might "correct" the period down to
        % (almost) zero. This is of course not correct
        if isequal(curvefile, @single_shooting)
          if trailpoint.x(end-1)/currpoint.x(end-1) < 0.5
            print_diag(0, 'contL: period decreasing too much\n');
            reduce_stepsize = 1;
          end
        end
      end

      %% D. CIS Processing
      if ~reduce_stepsize
        trailpoint.CISdata = feval(...
          cds.curve_CIS_step, trailpoint.x, currpoint.CISdata);
        if isempty(trailpoint.CISdata)
          print_diag(1, 'contL: Candidate step failed\n');
          reduce_stepsize = 1;
        end
      end

      %% E. Test Function Evaluation
      trailpoint.tvals = [];
      if ~reduce_stepsize && Singularities
        cds.previous_tvals = currpoint.tvals;
        [trailpoint.tvals, failed] = EvalTestFunc(0, trailpoint);
        if failed
          print_diag(1, ...
            'contL: Unable to evaluate Test Functions at Point %d', cds.i);
          reduce_stepsize = 1;
        end
      end

      %% F. Detect singularities
      NeedToLocate = 0;
      if ~reduce_stepsize && Singularities

        % WM: the testvals arrays are not copied anymore, instead
        % WM: use sign function and compare instead of multiply (for speed).

        print_diag(4,'comparing values of test functions:\n')
        print_diag(4,'currpoint :');
        print_diag(4,'%.4f ', currpoint.tvals);
        print_diag(4,'\n');
        print_diag(4,'trailpoint:');
        print_diag(4,'%.4f ', trailpoint.tvals);
        print_diag(4,'\n');

        signchanges = sign(trailpoint.tvals) ~= sign(currpoint.tvals);
        changes     =      trailpoint.tvals  ~=      currpoint.tvals;
        if any(signchanges) || any(changes)
          % Every crossing that is required occurs
          % Every crossing that is not required does not occur
          % Required crossings matrix:
          S_true   = double(cds.S(:,cds.ActTest)' == Constants.sign_change);
          % Required noncrossings matrix:
          S_false  = double(cds.S(:,cds.ActTest)' == Constants.sign_constant);
          % DV: Required to change value (not sign):
          S_change = double(cds.S(:,cds.ActTest)' == Constants.value_change);  
          all_sings_detected = ...
              (  signchanges * S_true       == sum(S_true)       ) & ...
              ( ~signchanges * S_false      == sum(S_false)      ) & ...
              (  changes     * S_change     == sum(S_change)     );
          singsdetected = [];
          singsdetected(cds.ActSing) = ...
            all_sings_detected(cds.ActSing); %#ok<AGROW>
          % The AGROW warning is a false positive, hence it is ignored with
          % %#ok<AGROW>. The array singsdetected does not grow.

          if sum(singsdetected) > 1 || more_than_one_Neimark_Sacker( ...
                                          currpoint.tvals, trailpoint.tvals)
            print_diag(3, 'More than one singularity detected: ')
            reduce_stepsize = 1;
          elseif sum(singsdetected) == 1 % exactly one singularity was detected
            si = find(singsdetected==1);  
            NeedToLocate = 1;
          end
        end
      end

      %% G. User Function Evaluation
      trailpoint.uvals = [];
      if ~reduce_stepsize && Userfunctions
        [trailpoint.uvals, failed] = feval( ...
                                      cds.curve_userf, 0, UserInfo, trailpoint.x);
        failed = failed || isempty(trailpoint.uvals);
        if failed
          format_str = 'contL: Unable to evaluate User Functions at Point %d: ';
          print_diag(1, format_str, cds.i);
          reduce_stepsize = 1;
        end
      end

      %% I. Reduce stepsize if needed
      if reduce_stepsize
        print_diag(3, 'Reducing Stepsize\n')
        if cds.h > cds.h_min
          cds.h = max(cds.h_min, cds.h * cds.h_dec_fac);
          corrections = corrections + 1;
          continue;
        else
          print_diag(0, 'Current Stepsize is too small: %d.\n', cds.h);
          cds.lastpointfound = 1;
          break;
        end
      else
        % accept step
        break;
      end

    end % end of while true, this loop runs until a continuation step is 
        % accepted or the continuation stops because the stepsize cannot be
        % reduced anymore (because the minimum stepsize has been reached).



    if cds.lastpointfound % Carel Jonkhout
      % occurs if stepsize too small
      continue
    end

    %% Location of singularities
    if NeedToLocate
      % At this point in the continuation loop, we only have to deal with one
      % singularity at a time. If multiple singularities had been detected in
      % one step, the stepsize has been decreased in the previous parts of the
      % continuation loop until only one singularity is detected.

      format_str = 'contL: %s detected at step %d\n';
      print_diag(1, format_str, cds.SingLables(si,:),cds.i);

      % do we have a locator?
      locatorAvailable = ismember( si, find(UseLocators == 1) );  

      if locatorAvailable

        p1 = currpoint;
        p2 = trailpoint;

        % Use locator
        print_diag(3,'using locator\n');
        singpoint = feval(cds.curve_locate, si, p1, p2);

        if isempty(singpoint) % Singularity not located
          print_diag(3, ['contL: Unable to locate %s with locator.\n' ...
            'Try increasing contL_Loc_MaxCorrIters, ' ...
            'contL_Loc_MaxCorrIters, contL_Loc_FunTolerance' ...
            'contL_Loc_VarTolerance\n'], cds.SingLables(si,:));              
          % attempt to locate with bisection:
          singpoint = LocateSingularity(si, p1, p2);          
        end
      else
        singpoint = LocateSingularity(si, currpoint, trailpoint);
      end

      if ~ isempty(singpoint)
        % process singularity
        singpoint.CISdata = feval(...
          cds.curve_CIS_step, singpoint.x, currpoint.CISdata);
        if isempty(singpoint.CISdata)
          print_diag(0, 'CIS algorithm failed at singular point')
        end
        singpoint.R = normU(cds.curve_func(singpoint.x, singpoint.CISdata));

        cds.i = cds.i + 1;
        s.label = cds.SingLables(si,:);

        [singpoint.tvals, ~] = EvalTestFunc(0, singpoint); 
        [failed ,s] = feval(cds.curve_process, si, singpoint, s);

        if ~ failed
          singpoint.h = norm(currpoint.x - singpoint.x);
          singpoint.uvals = [];
          if Userfunctions
            [singpoint.uvals, ~] = cds.curve_userf(0, UserInfo, singpoint.x);
          end
          DefaultProcessor(singpoint, s);
        end
      else
        print_diag(0, 'Unable to locate %s\n', cds.SingLables(si,:));
      end
    end

    %% D. User Functions
    userchanges = sign(trailpoint.uvals) ~= sign(currpoint.uvals);
    if Userfunctions && any(userchanges)
      useridx = find(userchanges);
      for ti=1:size(useridx,2)
        id = useridx(ti);
        print_diag(1,'\nStep %d: User Function %s detected ... \n', ...
          cds.i+1, UserInfo{id}.label);
        [userpoint] = LocateUserFunction(id, UserInfo, currpoint, trailpoint);

        if isempty(userpoint)
          print_diag(0, ...
            'Unable to locate User Function %s \n', UserInfo{id}.label);
        else
          % Point found
          cds.i = cds.i+1; s = [];
          s.label = UserInfo{id}.label;
          userpoint.h = norm(currpoint.x - userpoint.x);
          [userpoint.uvals, ~] = feval(cds.curve_userf, 0,UserInfo,userpoint.x);
          s.data.userfunctions = userpoint.uvals;

          if Singularities
            userpoint.CISdata = feval(...
                    cds.curve_CIS_step, userpoint.x, currpoint.CISdata);
            [userpoint.tvals, ~] = EvalTestFunc(0, userpoint);
            s.data.testfunctions = userpoint.tvals;
          end

          s.msg  = sprintf('%s',UserInfo{id}.name);

          DefaultProcessor(userpoint, s);
        end
      end
    end

    %%      E. Process

    % DV: check if testfunctions should be updated
    tfUpdate = isfield(trailpoint, 'CISdata');
    tfUpdate = tfUpdate && ~isempty(trailpoint.CISdata);
    tfUpdate = tfUpdate && isstruct(trailpoint.CISdata);
    tfUpdate = tfUpdate && ...
      (        (contopts.CIS_AdaptOnOverlap && trailpoint.CISdata.overlap) ...
          ||   (Singularities && sum(singsdetected) > 0) ...
      );

    if (mod(cds.i,AdaptSteps)==0) || tfUpdate

      [has_changed,x2,v2,trailpoint.CISdata] = feval(cds.curve_adapt, ...
          trailpoint.x, trailpoint.v, trailpoint.CISdata, tfUpdate);
      trailpoint.x = x2;
      trailpoint.v = v2;

      if has_changed && Singularities
        % recompute testvals            
        [trailpoint.tvals,~] = EvalTestFunc(0,trailpoint);
      end
    end
    
    cds.i = cds.i + 1;
    if ~ isempty(additional_arguments.stopping_condition) && ...
                 additional_arguments.stopping_condition(currpoint)
      cds.lastpointfound = 1;
      DefaultProcessor(trailpoint);
      print_diag(1,'Stopping condition reached\n');
      break;
    end
    
    % closed curve check
    if CheckClosed>0 && cds.i>= CheckClosed && norm(trailpoint.x-x0) < cds.h
      cds.lastpointfound = 1;
      currpoint = first_point;
      print_diag(1, '\nClosed curve detected at step %d\n', cds.i);
      DefaultProcessor(currpoint);
      break;
    end
    
    
    trailpoint  = DefaultProcessor(trailpoint);
    if ~ isempty(additional_arguments.callback)
      additional_arguments.callback(currpoint, trailpoint)
    end
    if trailpoint.v'*currpoint.v < 0,  trailpoint.v = -trailpoint.v; end
    currpoint = trailpoint;

    % stepsize control
    if cds.h < cds.h_max && corrections==1
      cds.h = min(cds.h*cds.h_inc_fac, cds.h_max);
    end

    if contopts.pause && mod(cds.i, contopts.nsteps_before_pause) == 0
      pause
    end

    print_diag(1,'time to compute step %.3f\n',toc(step_start_time));

  end % end of main continuation loop

  %% III. Finalization
  sout = cds.sout;
  cds.sout = [];

  fclose(cds.dataFID);
  if(cds.logFID~=1)
      fclose(cds.logFID);
  end
end % end of contL

function v0 = find_initial_tangent_vector(x0, v0, CISdata0)
  global cds contopts
  if isempty(v0)
    A = contjac(x0, CISdata0);
    if isempty(A); v0 = []; return; end
    if contopts.switch
      kernel = null(A);
      v0 = kernel(:,1);
    else
      v = zeros(length(x0),1);
      v(end) = 1;
      B = [A; v'];
      C = zeros(length(x0),1);
      C(end) = -1;
      if isequal(cds.curve, @limitcycleL)
        v0 = linear_solver_collocation(B,C);
      else
        v0 = bordCIS1(B,C,1);
      end
    end
    if contopts.newtcorrL_use_max_norm
      v0 = v0 / max(abs(v0));
    else
      v0 = v0 / norm(v0);
    end
  end
end

%------------------------------------------
%
%  Evaluate testfunctions
%
%------------------------------------------

function [out,failed] = EvalTestFunc(id, point)
  global cds

  if id == 0
    % WM: evaluate all testfunctions at once
    [out,failed] = feval(...
            cds.curve_testf, cds.ActTest, point.x, point.v, point.CISdata);
  elseif isempty(id)  % DV
    out = [];
    failed = 0;
  else
    [out,failed] = feval(cds.curve_testf, id, point.x, point.v, point.CISdata);
  end
  print_diag(5,'tf eval at period: %.16f param: %.16f\n', ...
                                             point.x(end-1),  point.x(end))
  print_diag(1,'Test Functions: [ ')
  if length(out) == size(cds.SingLables, 1)
    for i = 1:length(out)
      print_diag(1,'%s:%+.5e ',cds.SingLables(i,:),out(i))
    end
  else
    print_diag(1,'%+.5e ',out)
  end
  print_diag(1,']\n')
end 


%------------------------------------------------------------
%
% Locate singularity xs between x1 and x2
% First locating zeros of testfunctions, then glue them
%
%------------------------------------------------------------
function singpoint = LocateSingularity(si, p1, p2)
  global cds contopts

  % zero locations of testf i is stored in testzero(:,i)

  % if 1 zero/tf check if nonzero/tf _is_ nonzero at that point if more zero/tf
  % then glue, nonzero/tf is kind of a problem because can always be found
  idx = find( cds.S(si,:)==0 );
  nzs = find( cds.S(si,:)==1 );
  nconst = find( cds.S(si,:)==2 );  % DV

  len = length(idx);

  xs = zeros(cds.ndim,len);
  vs = zeros(cds.ndim,len);
  Rs = zeros(1       ,len);

  SingLables    = cds.SingLables;

  % Find zeros of all needed testfunction zeros
  for ii = 1:len
    [psii, ~, ~] = LocateTestFunction(idx(ii), p1, p2);
    if isempty(psii)  
      format_string = 'Testfunction %d for Singularity %s failed to converge\n';
      print_diag(3, format_string, idx(ii), SingLables(si,:));
      singpoint = [];
      return;
    else
      xs(:, ii) = psii.x;
      vs(:, ii) = psii.v;
      Rs(:, ii) = psii.R;
    end
  end

  % Check that the zeros found are close enough
  switch len
  case 0
    % Oops, we have detected a singularity without a vanishing testfunction
    error('Internal error: trying to locate a non-detected singularity');      
  case 1
    % return xs and vs      
  otherwise
    nm = zeros(1,len);    
    for i=1:len
      nm(i) = norm(xs(:,i));
    end    
    if max(nm)-min(nm) < contopts.contL_Testf_FunTolerance
      xs = mean(xs,2);
      vs = mean(vs,2);
      Rs = max(Rs);
    else
      format_string = ['LocateSingularity:' ...
        ' Testfunctions for %s vanish at different points \n'];
      print_diag(3, format_string, SingLables(si,:));
      singpoint = [];
      return;
    end
  end

  % check non-zero testfunctions
  tval1 = EvalTestFunc(nzs, p1);
  tval2 = EvalTestFunc(nzs, p2);
  if any(sign(tval1(nzs)) ~= sign(tval2(nzs)))
    format_string = 'LocateSingularity: Nonzero testfunction for %s vanishes \n';
    print_diag(3, format_string, SingLables(si,:));
    singpoint = [];
  end

  % check for required changes
  tval1 = EvalTestFunc(nconst, p1);
  tval2 = EvalTestFunc(nconst, p2);
  if any(sign(tval1(nconst)) == sign(tval2(nconst)))
    format_string = ['LocateSingularity:' ...
      ' Testfunction for %s does not change its value \n'];
    print_diag(3, format_string, SingLables(si,:));
    singpoint = [];
  end

  singpoint.x = xs;
  singpoint.v = vs;
  singpoint.R = Rs;
  singpoint.tvals = [];
  singpoint.uvals = [];
  singpoint.angle = 0;
end

%----------------------------
%
% DefaultProcessor
%
%----------------------------

function point = DefaultProcessor(varargin)
  global cds
  point = feval(cds.curve_defaultprocessor, varargin{:});
end

%----------------------------------------------
%
%LocateTestFunction(id,x1,v1,x2,v2)
%
%----------------------------------------------
function [pout, p1, p2] = LocateTestFunction(id, p1, p2)
  bisection_start_time = tic;
  global contopts cds

  pout = [];
  MaxTestIters  = contopts.MaxTestIters;

  FunTolerance  = contopts.contL_Testf_FunTolerance;
  VarTolerance  = contopts.contL_Testf_VarTolerance;


  % default locator: bisection
  %print_diag(3,'Locating by test function %d\n', id);

  [t1, failed1]   = EvalTestFunc(id, p1);
  [t2, failed2]   = EvalTestFunc(id, p2);

  if ((~isempty(failed1)) || (~isempty(failed2))) && (failed1 || failed2)
    print_diag(3, 'Evaluation of testfunctions failed in bisection');
    pout = []; p1 = []; p2 = 0;
    return;
  end

  tmax = 10 * max(abs(t1), abs(t2));
  
  p = 1;
  is_within_tolerances = false;
  for i = 1 : MaxTestIters
    if tmax < Inf
      % WM: make educated guess of where the zero point might be
      r = abs(t1(id)/(t1(id)-t2(id)))^p;
      if r <= 0.1 || r >= 0.9
        r=0.5;
      end
    else
      r=0.5;          % r = 0.5; %  -> 'normal' way
    end
    %     r = 0.5; % DV:TEST

    x3 = p1.x + r*(p2.x-p1.x);
    v3 = p1.v + r*(p2.v-p1.v);
    v3 = v3/norm(v3);
    CISdata3 = [];
    if cds.newtcorrL_needs_CISdata 
      CISdata3 = feval(cds.curve_CIS_step, x3, p1.CISdata);
      if isempty(CISdata3)
        print_diag(3, 'CIS algorithm failed during bisection')
        return
      end
    end
    if contopts.NewtonPicard
      p3 = NP_do_corrections(x3,v3);
      if isempty(p3)
        print_diag(3, 'Newton-Picard corrections failed during bisection\n')
        return        
      end
      if id == Constants.LPC_id 
        % if we are locating a limit point of cycles (LPC)
        p3.v = NP_find_tangent_vector(cds.curve, p3.x, p3.v);
      end
    else
      p3 = newtcorrL(x3,v3, CISdata3);
      if isempty(p3)
        print_diag(3, 'newtcorrL algorithm failed during bisection\n')
        return
      end
    end

    p3.CISdata = feval(cds.curve_CIS_step, p3.x, p1.CISdata);
    [tval, failed] = EvalTestFunc(id,p3);
    if failed
      print_diag(3, 'Evaluation of testfunctions failed in bisection\n');
      return;
    end

    %MP: Changed to make the check work for very small values of x 07/2019
    % dist1 = 2*norm(p3.x-p1.x)/(norm(p1.x)+norm(p3.x));
    % dist2 = 2*norm(p3.x-p2.x)/(norm(p2.x)+norm(p3.x));
    dist1 = 2*norm(p3.x-p1.x)/(norm(p1.x)+norm(p3.x) + 1);
    dist2 = 2*norm(p3.x-p2.x)/(norm(p2.x)+norm(p3.x) + 1);

    if abs(tval(id)) > tmax
      print_diag(3, ['Could not locate the zero of testfunction %d, ' ... 
              'because the testfunction is not monotonous between the two' ...
              'points where the zero of the testfunction was the detected.\n']);
      break;
    end

    if (abs(tval(id)) <= FunTolerance && min(dist1,dist2) <= VarTolerance)
      is_within_tolerances = true;
      pout = p3;
      break;
    elseif sign(tval(id)) == sign(t2(id))
      p2 = p3;
      t2 = tval;
      p  = 1.02;
    else
      p1 = p3;
      t1 = tval;
      p  = 0.98;
    end
  end % end of for i = 1 : MaxTestIters

  if ~ is_within_tolerances
    if abs(tval(id)) > FunTolerance
      print_diag(3, ['Could not reach the specified FunTolerance. ' ...
                     'FunTolerance: %.3e test function value: %.3e ' ...
                     'MaxTestIters: %d\n'], ...
                     FunTolerance, tval(id), MaxTestIters);
    end
    if min(dist1,dist2) > VarTolerance
       print_diag(3, ['Could not reach the specified VarTolerance. ' ...
                     'VarTolerance: %.3e ', ...
                     'upper bound on distance to singularity: %.3e ' ...
                     'MaxTestIters: %d\n'], ...
                      VarTolerance, min(dist1,dist2), MaxTestIters);
    end
    pout = [];
  end 

  print_diag(1,'Time spent in bisection: %f\n', toc(bisection_start_time));
end % of LocateTestFunction
%----------------------------------------------
%
%LocateUserFunction(id,x1,v1,x2,v2)
%
%----------------------------------------------
function pout = LocateUserFunction(id, userinf, p1, p2)
  % default locator: bisection
  global cds contopts
  
  pout = [];

  i = 1;
  t1 = feval(cds.curve_userf, id, userinf, p1.x);
  t2 = feval(cds.curve_userf, id, userinf, p2.x);
  tmax = 10*max(abs(t1),abs(t2));
  p = 1;

  MaxTestIters  = contopts.contL_Userf_MaxIters;
  FunTolerance  = contopts.contL_Userf_FunTolerance;
  VarTolerance  = contopts.contL_Userf_VarTolerance;

  while i <= MaxTestIters
    if tmax < Inf
      % WM: make educated guess of where the zero point might be
      r = abs(t1/(t1-t2))^p;
    else
      r = 0.5;
    end
    %  r = 0.5; %  -> 'normal' way
    x3 = p1.x + r*(p2.x-p1.x);
    v3 = p1.v + r*(p2.v-p1.v);
    v3 = v3/norm(v3);CISdata3 = [];
    if cds.newtcorrL_needs_CISdata 
      CISdata3 = feval(cds.curve_CIS_step, x3, p1.CISdata);
      if isempty(CISdata3)
        print_diag(3, 'CIS algorithm failed during bisection')
        return
      end
    end
    if contopts.NewtonPicard
      error(['The use of Userfunctions in combination with Newton-Picard ' ...
             'is not supported']);
    else
      p3 = newtcorrL(x3,v3, CISdata3);
    end
    if isempty(p3)
      return
    end
    
    p3.CISdata = feval(cds.curve_CIS_step, p3.x, p1.CISdata);
    [tval, failed] = feval(cds.curve_userf, id, userinf, p3.x);
    if failed
      print_diag(3, 'Evaluation of testfunctions failed in bisection');
      return;
    end
    
    %JH: Changed to make the check for relative difference. 9/25/06
    %dist1 = norm(x-x1);
    %dist2 = norm(x-x2);
    dist1 = 2*norm(p3.x-p1.x)/(norm(p1.x)+norm(p3.x));
    dist2 = 2*norm(p3.x-p2.x)/(norm(p2.x)+norm(p3.x));
    
    if abs(tval(id)) > tmax
      print_diag(3,'Testfunction behaving badly.\n');
      break;
    end
    
    if (abs(tval(id)) <= FunTolerance && min(dist1,dist2) < VarTolerance)
      pout = p3;
      break;
    elseif sign(tval(id))==sign(t2(id))
      p2 = p3;
      t2 = tval;
      p  = 1.02;
    else
      p1 = p3;
      t1 = tval;
      p  = 0.98;
    end
  end % end of while i <= MaxTestIters
end % end of LocateUserFunction