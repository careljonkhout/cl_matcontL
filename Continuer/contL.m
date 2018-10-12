function [sout,datafile]=contL(curvefile, x0, v0, opts)
%
% CONTINUE(cds.curve, x0, v0, options)
%
% Continues the curve from x0 with optional directional vector v0
% options is a option-vector created with CONTSET
% The first two parameters are mandatory.
global cds contopts

%% I. Initialization
if size(x0, 2) > 1
    error('Initial point must be a column vector')
elseif isempty(x0)
    error('Initial point must be non-empty')
end
if isempty(opts); contopts = contset(); else; contopts = opts; end

loadCurveFile(curvefile);
feval(cds.curve_options);

[datafile, ~]  = openFiles();

AdaptSteps     = contopts.Adapt;
CheckClosed    = contopts.CheckClosed;
%Eigenvalues    = contopts.Cont_Eigenvalues;
MaxNumPoints   = contopts.MaxNumPoints;
Singularities  = contopts.Singularities;
SmoothingAngle = contopts.contL_SmoothingAngle;
Userfunctions  = contopts.Userfunctions;
IgnoreSings    = contopts.IgnoreSingularity;
UseLocators    = contopts.Locators;
UserInfo       = contopts.UserFuncInfo;
if contopts.contL_ParallelComputing && isempty(gcp('NoCreate'))
    parpool;   % initialize new parallel pool when no is available
end

cds.h          = contopts.InitStepsize;
cds.h_max      = contopts.MaxStepsize;
cds.h_min      = contopts.MinStepsize;
cds.h_inc_fac  = contopts.h_inc_fac;  
cds.h_dec_fac  = contopts.h_dec_fac;  
cds.tfUpdate       = 0;
cds.i              = 1;
cds.lastpointfound = 0;

%% determine active singularities and testfunctions
cds.nActTest = 0;
if Singularities
    [cds.S , cds.SingLables] = feval(cds.curve_singmat);
    [nSing , nTest         ] = size(cds.S);
    cds.S(IgnoreSings,:) = 8;
    cds.ActSing = setdiff(1:nSing, IgnoreSings);
    cds.nActSing = length(cds.ActSing);
    cds.ActTest  = find( sum((cds.S~=8),1) > 0 );
end

%% userfunctions
if Userfunctions
    if length(contopts.UserFuncInfo) == length(cds.userf)
        cds.nUserf = length(cds.userf);     % DV: Does this work properly??
    else
        error('Wrong user info specified');
    end
end

%% Algorithm starts here
cds.StartTime = clock;

feval(cds.curve_init, x0, v0); % DV 2018
cds.newtcorrL_needs_CISdata = 0;

if (~ isequal(curvefile,@limitcycleL))

    try feval(cds.curve_func    , x0); catch; cds.newtcorrL_needs_CISdata = 1; end
    try feval(cds.curve_jacobian, x0); catch; cds.newtcorrL_needs_CISdata = 1; end

    CISdata0 = [];
    if isempty(v0)
        if cds.newtcorrL_needs_CISdata
            CISdata0 = feval(cds.curve_CIS_first_point, x0);
            if isempty(CISdata0); print_diag(0,'contL: failed to intialize CIS algorithm.\n'); sout = []; return; end
        end
        v0       = find_initial_tangent_vector(x0, v0, CISdata0);
        if isempty(v0); print_diag(0,'contL: failed to find initial tangent vector.\n'); sout = []; return; end
    end


    %% Newton corrections for first point
    firstpoint = newtcorrL(x0, v0, CISdata0);
    
else
    if isempty(v0)        
        [x0, v0] = CorrectStartPoint(x0, v0);
    end
    firstpoint.x = x0;
    firstpoint.v = v0;
    firstpoint.R = 0;
    firstpoint.tvals = [];
    firstpoint.uvals = [];
    if isempty(x0)
        print_diag(0,'contL: no convergence at x0.\n');  
        return;            
    end
end

if isempty(firstpoint)
    print_diag(0,'contL: no convergence at x0.\n');
    sout = [];
    return;
end
firstpoint.h = cds.h;

%% CIS data algorithm
% if contopts.CIS_UsingCIS                    % MP 2018
    firstpoint.CISdata = feval(cds.curve_CIS_first_point, x0);
    if isempty(firstpoint.CISdata)
        print_diag(0,'contL: failed to intialize CIS algorithm.\n'); 
        sout = []; 
        return; 
    end
% end                                          % MP 2018


%% Direction Vector Determination
Backward       = contopts.Backward;
if abs(v0(end)) < 1e-6; Vdir = sign(sum(v0(1:end-1))); else; Vdir = sign(v0(end)); end
if (Backward && Vdir > 0) || (~Backward && Vdir < 0)
    firstpoint.v = -firstpoint.v;
end

 
firstpoint = DefaultProcessor(firstpoint); 

%% Test and user functions
if Singularities
    % WM: calculate all testfunctions at once
    [firstpoint.tvals,failed] = EvalTestFunc(0, firstpoint);
    if failed
        print_diag(0,'contL: Evaluation of test functions failed at start point.\n');
        sout = [];  datafile=[];
        return
    end
end
if Userfunctions
    [firstpoint.uvals,failed]   = feval(cds.curve_userf, 0, UserInfo, x0);
    if failed
        print_diag(0,'contL: Evaluation of user functions failed at start point.\n');
        sout = [];  datafile=[];
        return
    end
end



%% II. Main Loop
currpoint = firstpoint;
while cds.i < MaxNumPoints && ~cds.lastpointfound
    corrections = 1;
    print_diag(1,'\n --- Step %d ---\n',cds.i);
    while 1
        
        %% A. Predict
        xpre = currpoint.x + cds.h * currpoint.v(1:cds.ndim);
        reduce_stepsize = 0;
        
        %% B. Correct
        trialpoint = newtcorrL(xpre, currpoint.v, currpoint.CISdata);
        if isempty(trialpoint)
            %print_diag(3, 'contL: newtcorrL failed\n ')
            reduce_stepsize = 1;
        end
        
        % curve smoothing
        if ~reduce_stepsize
            trialpoint.h = cds.h;
            trialpoint.angle = innerangle(currpoint.v,trialpoint.v);
            if trialpoint.angle > SmoothingAngle
                print_diag(2, 'contL: Innerangle too large\n ');
                reduce_stepsize = 1;
            end
        end
        
        % CIS Processing
        if ~reduce_stepsize
            special_step = 0;
            if ~isempty(cds.curve_CIS_step)
                trialpoint.CISdata = feval(cds.curve_CIS_step, trialpoint.x, currpoint.CISdata);
                if isempty(trialpoint.CISdata)
                    print_diag(1, 'contL: Candidate step failed\n ');
                    reduce_stepsize = 1;
                end
            end
        end
        
        % Test Function Evaluation
        if ~reduce_stepsize && Singularities
            [trialpoint.tvals,failed] = EvalTestFunc(0, trialpoint);
            if failed
                print_diag(1, 'contL: Unable to evaluate Test Functions at Point %d: ',cds.i)
                reduce_stepsize = 1;
            end
        end
        
        % Detect singularities
        NeedToLocate = 0;
        if ~reduce_stepsize && Singularities
            singsdetected = [];% WM: the testvals arrays are not copied anymore, instead
            % WM: use sign function and compare instead of multiply (for speed).
            testchanges = sign(trialpoint.tvals) ~= sign(currpoint.tvals);
            testchanges2= trialpoint.tvals ~= currpoint.tvals;
            if any(testchanges)
                
                % Every crossing that is required occurs
                % Every crossing that is not required does not occur
                S_true  = double(cds.S(:,cds.ActTest)' == 0);  % Required crossings matrix
                S_false = double(cds.S(:,cds.ActTest)' == 1);  % Required noncrossings matrix
                S_change= double(cds.S(:,cds.ActTest)' == 2);  % DV: Required to change value (not sign)
                
                all_sings_detected = ...
                    (  testchanges * S_true  == sum(S_true)  ) & ...
                    ( ~testchanges * S_false == sum(S_false) ) & ...
                    ( testchanges2 * S_change== sum(S_change)); % DV
                singsdetected(cds.ActSing) = all_sings_detected(cds.ActSing);
                
                if sum(singsdetected) > 1
                    print_diag(3, 'More than one singularity detected: ')
%                     print_diag(0, 'contL: More than one singularity detected:\n')
                    reduce_stepsize = 1;
                elseif sum(singsdetected) == 1
                    if special_step
                        print_diag(3, 'Singularity detected at special step: ')
                        reduce_stepsize = 1;
                    else
                        si = find(singsdetected==1);  % singularity detected!
                        NeedToLocate = 1;
                    end
                end
            end
        end
        
        % User Function Evaluation
        if ~reduce_stepsize && Userfunctions
            [trialpoint.uvals,failed] = feval(cds.curve_userf, 0, UserInfo, trialpoint.x);
            failed = failed || isempty(trialpoint.uvals);
            if failed
                print_diag(1, 'contL: Unable to evaluate Test Functions at Point %d: ',cds.i)
                reduce_stepsize = 1;
            end
        end
        
        % Detect userfunctions
        if ~reduce_stepsize && Userfunctions
            % WM: use sign function and compare instead of multiply (for speed).
            userchanges = sign(trialpoint.uvals) ~= sign(currpoint.uvals);
            if special_step && any(userchanges)
                print_diag(3, 'Userfunction occus in Special Step: \n');
                reduce_stepsize = 1;
            end
        end
        
        % Reduce stepsize
        if reduce_stepsize
            % reduce stepsize
            print_diag(3, 'Reducing Stepsize\n')
            if cds.h > cds.h_min
                cds.h = max(cds.h_min, cds.h*cds.h_dec_fac);
                corrections = corrections+1;
                continue;
            else
                print_diag(0,'Current Stepsize is too small: %d.\n',cds.h);
                cds.lastpointfound = 1;
                break;
            end
        else
            % accept step
            break;
        end
    end
    
    %% Location of singularities
    if NeedToLocate
        % DV: There is always only one singularity detected
        %print_diag(1,'\n%s detected at step %d\n',cds.SingLables(si,:),cds.i);
        print_diag(1,'contL: \n%s detected at step %d\n',cds.SingLables(si,:),cds.i);
        
        locatorAvailable = ismember(si, find(UseLocators==1));  % do we have a locator?
        
        if locatorAvailable
            
            p1 = currpoint;
            p2 = trialpoint;
            
            % Use locator
            print_diag(3,'using locator\n');
            singpoint = feval(cds.curve_locate, si, p1, p2);
            
            if isempty(singpoint) % Singularity not located
                % MP print_diag(3,'Unable to locate %s with locator. \nTry increasing Locator_MaxIters, Loc_MaxIters, Loc_FunTolerance and Loc_VarTolerance\n', cds.SingLables(si,:));
                print_diag(3,'contL: Unable to locate %s with locator. \nTry increasing contL_Loc_MaxCorrIters, Loc_MaxIters, contL_Loc_FunTolerance and contL_Loc_VarTolerance\n', cds.SingLables(si,:));              
                singpoint = LocateSingularity(si, p1, p2);          % DV: Try location with testfunctions
            end
        else
            singpoint = LocateSingularity(si, currpoint, trialpoint);
        end
        
        if ~isempty(singpoint)
            % process singularity
            singpoint.CISdata = feval(cds.curve_CIS_step, singpoint.x, currpoint.CISdata);
            if isempty(singpoint.CISdata); print_diag(0, 'CIS algorithm failed at singular point'); end
            singpoint.R = normU(feval(cds.curve_func, singpoint.x, singpoint.CISdata));
            
            cds.i = cds.i + 1; s = [];
            s.label = cds.SingLables(si,:);
            
            [singpoint.tvals,failed] = EvalTestFunc(0, singpoint);
            [failed,s] = feval(cds.curve_process, si, singpoint, s );

            if ~failed
                singpoint.h = norm(currpoint.x - singpoint.x);
                if Userfunctions
                    [point.uvals,failed] = feval(cds.curve_userf, 0, UserInfo, singpoint.x);
                end
                singpoint = DefaultProcessor(singpoint, s);
            end
        else
            print_diag(0, 'Unable to locate %s\n', cds.SingLables(si,:));
        end
    end
    
    %% D. User Functions
    if Userfunctions
        if any(userchanges)
            useridx = find(userchanges);
            for ti=1:size(useridx,2)
                id = useridx(ti);
                print_diag(1,'\nStep %d: User Function %s detected ... \n', ...
                    cds.i+1,UserInfo{id}.label);
                [userpoint] = LocateUserFunction(id, UserInfo, currpoint, trialpoint);
                
                if isempty(userpoint)
                    % print_diag(1,'Unable to locate User Function %s \n', UserInfo{id}.label); % MP
                    print_diag(0,'Unable to locate User Function %s \n', UserInfo{id}.label);  % MP
                else % Point found
                    cds.i = cds.i+1; s = [];
                    s.label = UserInfo{id}.label;
                    
                    userpoint.h = norm(currpoint.x - userpoint.x);
                    
                    [userpoint.uvals,failed] = feval(cds.curve_userf, 0,UserInfo,userpoint.x);
                    s.data.userfunctions = userpoint.uvals;
                    
                    if Singularities
                        userpoint.CISdata = feval(cds.curve_CIS_step, userpoint.x, currpoint.CISdata);
                        [userpoint.tvals,failed] = EvalTestFunc(0, userpoint);
                        s.data.testfunctions = userpoint.tvals;
                    end
                    
                    s.msg  = sprintf('%s',UserInfo{id}.name);
                    
                    userpoint = DefaultProcessor(userpoint, s);
                end
            end
        end
    end
    
    %%      E. Process
    
    % DV: check if testfunctions should be updated
    tfUpdate = isfield(trialpoint, 'CISdata') && ~isempty(trialpoint.CISdata) && isstruct(trialpoint.CISdata) && ...
        ((contopts.CIS_AdaptOnOverlap && trialpoint.CISdata.overlap) || (Singularities && sum(singsdetected) > 0));
    
    if (mod(cds.i,AdaptSteps)==0) || tfUpdate
            
        [res,x2,v2,trialpoint.CISdata] ...
            = feval(cds.curve_adapt, trialpoint.x, trialpoint.v, trialpoint.CISdata, tfUpdate);
        trialpoint.x = x2;
        trialpoint.v = v2;
        
        if res == 1 && Singularities
            
            % recompute testvals            
            [trialpoint.tvals,failed] = EvalTestFunc(0,trialpoint);
        end
    end
    
    if trialpoint.v'*currpoint.v < 0,  trialpoint.v = -trialpoint.v; end
    currpoint = trialpoint;
    
    cds.i = cds.i + 1;
    currpoint  = DefaultProcessor(currpoint);

    % stepsize control
    if cds.h < cds.h_max && corrections==1
        cds.h = min(cds.h*cds.h_inc_fac, cds.h_max);
    end
    
    % closed curve check
    if CheckClosed>0 && cds.i>= CheckClosed && norm(trialpoint.x-x0) < cds.h
        cds.i=cds.i+1;
        cds.lastpointfound = 1;
        currpoint = firstpoint;
        print_diag(1,'\nClosed curve detected at step %d\n', cds.i);
        currpoint = DefaultProcessor(currpoint);
        break;
    end
end
%% III. Finalization
sout = cds.sout;

fclose(cds.dataFID);
if(cds.logFID~=1)
    fclose(cds.logFID);
end

%--< END OF CONTINUER >--



function v0 = find_initial_tangent_vector(x0, v0, CISdata0)
if isempty(v0)
    A = contjac(x0, CISdata0);
    if isempty(A); v0 = []; return; end
    v = zeros(length(x0),1);
    v(end) = 1;
    B = [A; v'];
    C = zeros(length(x0),1);
    C(end) = -1;
    v0 = bordCIS1(B,C,1);
    v0 = v0/norm(v0);
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
    [out,failed] = feval(cds.curve_testf, cds.ActTest, point.x, point.v, point.CISdata);
elseif isempty(id)  % DV
    out = [];
    failed = 0;
else
    [out,failed] = feval(cds.curve_testf, id, point.x, point.v, point.CISdata);
    %out = out(id);
end

print_diag(1,'Test Functions: [')
print_diag(1,' %+.3e',out)
print_diag(1,']\n')

%--< END OF TESTF EVAL >--


%------------------------------------------------------------
%
% Locate singularity xs between x1 and x2
% First locating zeros of testfunctions, then glue them
%
%------------------------------------------------------------

function singpoint = LocateSingularity(si, p1, p2)
global cds contopts

% zero locations of testf i is stored in testzero(:,i)

% if 1 zero/tf check if nonzero/tf _is_ nonzero at that point
% if more zero/tf then glue, nonzero/tf is kind of a problem because can always be found
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
        print_diag(3,'Testfunction %d for Singularity %s failed to converge\n', idx(ii), SingLables(si,:));
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
        
        if max(nm)-min(nm) < contopts.Loc_Testf_FunTolerance
            xs = mean(xs,2);
            vs = mean(vs,2);
            Rs = max(Rs);
        else
            print_diag(3, 'LocateSingularity: Testfunctions for %s vanish at different points \n', SingLables(si,:));
            singpoint = [];
            return;
        end
end

% check non-zero testfunctions
tval1 = EvalTestFunc(nzs, p1);
tval2 = EvalTestFunc(nzs, p2);
if any(sign(tval1(nzs)) ~= sign(tval2(nzs)))
    print_diag(3, 'LocateSingularity: Nonzero testfunction for %s vanishes \n', SingLables(si,:));
    singpoint = [];
end

% check for required changes
tval1 = EvalTestFunc(nconst, p1);
tval2 = EvalTestFunc(nconst, p2);
if any(sign(tval1(nconst)) == sign(tval2(nconst)))
    print_diag(3, 'LocateSingularity: Testfunction for %s does not change its value \n', SingLables(si,:));
    singpoint = [];
end

singpoint.x = xs;
singpoint.v = vs;
singpoint.R = Rs;
singpoint.tvals = [];
singpoint.uvals = [];

%----------------------------
%
% DefaultProcessor
%
%----------------------------

function point = DefaultProcessor(varargin)
global cds

point = feval(cds.curve_defaultprocessor, varargin{:});

%--< END OF DEFAULTPROCESSOR >--

%----------------------------------------------
%
%LocateTestFunction(id,x1,v1,x2,v2)
%
%----------------------------------------------
%function [x,v,R,i,p1,p2] = LocateTestFunction(id, p1, p2, MaxIters, FunTolerance, VarTolerance)
function [pout,p1,p2] = LocateTestFunction(id, p1, p2)
tic
global contopts cds

pout = [];
MaxIters  = contopts.MaxTestIters;            % DV MP: Use new options

FunTolerance  = contopts.contL_Testf_FunTolerance;   % MP:
VarTolerance  = contopts.contL_Testf_VarTolerance;   % MP:

% default locator: bisection
%print_diag(3,'Locating by test function %d\n', id);

[t1, failed1]   = EvalTestFunc(id, p1);
[t2, failed2]   = EvalTestFunc(id, p2);

if ((~isempty(failed1)) || (~isempty(failed2))) && (failed1 || failed2)
    print_diag(3, 'Evaluation of testfunctions failed in bisection');
    x = []; v = []; i = 0;
    return;
end

tmax = 10 * max(abs(t1), abs(t2));
p    = 1;
failed2 = 1;
R = [];
for i = 1:MaxIters
    if tmax < Inf
        % WM: make educated guess of where the zero point might be
        r = abs(t1/(t1-t2))^p;
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
    p3 = newtcorrL(x3,v3, CISdata3);
    if isempty(p3)
        return
    end
    
    p3.CISdata = feval(cds.curve_CIS_step, p3.x, p1.CISdata);
    [tval, failed] = EvalTestFunc(id,p3);
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
        failed2 = 0;
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
end

if failed2  % DV: When MaxIters is reached without meeting tolerances
    pout = [];
end

print_diag(1,'Time spent in bisection: %f\n', toc);

%--< END OF locatetestfunction>--
%---------------------------------------------
%----------------------------------------------
%
%LocateUserFunction(id,x1,v1,x2,v2)
%
%----------------------------------------------
function pout = LocateUserFunction(id, userinf, p1, p2)
% default locator: bisection
global cds contopts
%print_diag(3,'locating by user function %d\n', id);
pout = [];

i = 1;
t1 = feval(cds.curve_userf, id, userinf, p1.x);
t2 = feval(cds.curve_userf, id, userinf, p2.x);
tmax = 10*max(abs(t1),abs(t2));
p = 1;

failed2  = 1;

MaxTestIters  = contopts.contL_Userf_MaxIters;            % DV: Use new options
FunTolerance  = contopts.contL_Userf_FunTolerance;
VarTolerance  = contopts.contL_Userf_VarTolerance;

while i<=MaxTestIters
    if tmax < Inf
        % WM: make educated guess of where the zero point might be
        r = abs(t1/(t1-t2))^p;
    else
        r=0.5;
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
    p3 = newtcorrL(x3,v3, CISdata3);
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
        failed2 = 0;
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
end

if failed2
    x=[];
    v=[];
end
%--< END OF locateuserfunction>--