function out = equilibriumL
%
% Equilibrium curve definition file for a problem in odefile
%
out{1}  = @curve_func;
out{2}  = @defaultprocessor;
out{3}  = @options;
out{4}  = @jacobian;
out{5}  = [];%@hessians;
out{6}  = @testf;
out{7}  = @userf;
out{8}  = @process;
out{9}  = @singmat;
out{10} = @locate;
out{11} = @init;
out{12} = @done;
out{13} = @adapt;

out{14} = @first_point;
out{15} = @candidate_step;
out{16} = @accept_step;

%% -------------------------------------------------------
function func = curve_func(arg)
global cds
%print_diag(5,'In equilibriumL/curvefunc\n');
[x,p] = rearr(arg); p = num2cell(p);
func = feval(cds.func,0, x, p{:});

%% ---------------------------------------------------------------
function jac = jacobian(x)
%print_diag(5,'In equilibriumL/jacobian\n');
[x0,p] = rearr(x); p = num2cell(p);

jac = [ejac(x0,p) ejacp(x0,p)];

%% ---------------------------------------------------------------
function failed = defaultprocessor(varargin)
global cds
%print_diag(5,'In equilibriumL/defaultprocessor\n');

if nargin > 1 && cds.i ~= 1
    cds.tfUpdate = 1;     % Need to update test functions
end
failed = savePoint(varargin{:});

%% -------------------------------------------------------------
function options                % DV
global contopts
% DV:  change options in the global variable if necessary
% DV: Do not use input and output arguments anymore

%% ---------------------------------------------------------------
function [out, failed] = testf(id, x, v)
global cds
global A1 Q1 T1 evl1_r evl1_l
%print_diag(5,'In equilibriumL/testf\n');

ndim = cds.ndim;
J = jacobian(x);
A3 = J(:,1:ndim-1);

[T,pout] = contCIS_step(cds.h, A1, Q1, T1, A3, ...
    cds.h_min, evl1_r, evl1_l);

if isempty(T)
    print_diag(3,'EquilibriumL: contCIS_step failed to find Asub');
    out = [];
    failed = 1;
    return;
end

evl_r = pout.evl_r;
NSub = cds.NSub;

out(4) = 0;
failed = 0;

for i = id
    lastwarn('');
    switch i
        case 1 % BP
            % Jacobian extended with bordering vectors v and w
            B = [J; v'];
            A = B(1:ndim-1,1:ndim-1);
            b = B(1:ndim-1,ndim);
            c = B(ndim,1:ndim-1)';
            d = B(ndim,ndim);
            P      = colamd(A);
            [L, U] = lu(A(:,P));            % A(:,P) = L*U, A(:,P)' = U'*L'
            w = L'\(U'\c(P,:));             % Solve A(:,P)'*w = c(P,:)
            deltas = d - w'*b;
            evl_r_min = min(abs((evl_r)));
            m_u = sum((real(evl_r)>0));
            out(1) = sign(deltas)*(-1)^m_u*evl_r_min;
        case 2 % H1
            k=1;
            bialt= zeros(NSub*(NSub-1)/2,1);
            for ii = 1:(NSub-1)
                for j = (ii+1):NSub
                    bialt(k) = evl_r(ii)+evl_r(j);
                    k = k + 1;
                end
            end
            mu_min = min(abs(bialt));
            M_u = sum((real(bialt)>0));
            out(2) = mu_min*(-1)^M_u;
        case 3 % H2
            out(3) = sum(real(pout.evl_r >= 0));  % DV
        case 4 % LP
            evl_r_min = min(abs((evl_r)));
            m_u = sum((real(evl_r)>0));
            out(4) = evl_r_min*(-1)^m_u;
        otherwise
            error('No such testfunction');
    end
    
    if ~isempty(lastwarn)
        failed = 1;
    end
    
end

%% ---------------------------------------------------------------
function [out, failed] = userf(id, UserInfo, x)
global cds

%print_diag(5,'In equilibriumL/userf\n');
failed = 0;

if id == 0 % evaluate all
    id = 1:cds.nUserf;
end

for ii = id
    lastwarn('');
    [x0,p] = rearr(x); p = num2cell(p);
    
    p{cds.ActiveParams} = x(end); %JH
    
    if (UserInfo{ii}.state==1)
        out(ii) = feval(cds.userf{ii},0,x0,p{:});
    else
        out(ii)=0;
    end
    
    if ~isempty(lastwarn)
        failed = 1;
        return
    end
end

%% ---------------------------------------------------------------
function [failed,s] = process(id, point, s)
global cds contopts
global A1 Q1 T1

%print_diag(5,'In equilibriumL/process\n');
x = point.x;
s.data.eigenvals = [];
switch id
    case 1 % BP
        s.msg     = 'Branch point';
        if isfield(cds.s, 'v2')
            s.data.v1 = cds.s.v1;
            s.data.v2 = cds.s.v2;
            if contopts.EQ_BranchingMethod == 2
                try
                    s.data.xnext = cds.s.xnext;
                    s.data.vnext = cds.s.vnext;
                catch
                    s.data.v1 = [];
                    s.data.v2 = [];
                end
            end
            s.data.eigenvals = cds.s.eigenvals;
        end
    case 2 % H
        if ~isempty(cds.s)
            borders.v = cds.s.borders.v;
            borders.w = cds.s.borders.w;
            A0        = cds.s.A0;
            Q0        = cds.s.Q0;
            T0        = cds.s.T0;
            kapa      = cds.s.kapa;
            cds.s = [];
        else
            A0 = contjac(x);
            A0 = A0(:,1:end-1);
            [T0, pout] = contCIS_step(cds.h, A1, Q1, T1, A0, cds.h_min);
            Q0 = pout.Q;
            
            [V,D1] = eig(T0);
            d      = diag(D1);
            % find a complex eigenvalue with smallest absolute real part and set omega
            % = its imaginary part
            I     = find(imag(d) ~= 0);
            if isempty(I)
                print_diag(1, 'Equilibrium process: All eigenvalues in start point are real \n');
                failed = 1; s = []; 
                return
            end
            Rmin  = min(abs(real(d(I))));
            I1    = find(abs(real(d)) == Rmin);
            lamda = real(d(I1(1)));
            omega = imag(d(I1(1)));				% get omega of Asub
            
            kapa  = lamda*lamda + omega*omega;
            
            [W,D2]= eig(T0');
            d     = diag(D2);
            I     = find(imag(d) ~= 0);
            Rmin  = min(abs(real(d(I))));
            I2    = find(abs(real(d)) == Rmin);
            
            % Input
            % make the bordered system nonsigular
            [Q,R,E] = qr([real(V(:,I1(1))) imag(V(:,I1(1)))]);
            borders.v =  [real(V(:,I1(1))) imag(V(:,I1(1)))];
            
            [Q,R,E] = qr([real(W(:,I2(1))) imag(W(:,I2(1)))]);
            borders.w =  [real(W(:,I2(1))) imag(W(:,I2(1)))];
        end
        
        % Normal form coefficient
        [x0,p] = rearr(x); p = num2cell(p);
        [s.data.l1, s.data.kapa] = nf_H_L(A0, x0, p, kapa);
        
        s.data.borders = borders;
        s.data.A0 = A0;
        s.data.T0 = T0;
        s.data.Q0 = Q0;
        
        % Eigenvalues
        NumEvals = cds.NSub + contopts.CIS_NExtra;
        s.data.eigenvals = eigs(A0, NumEvals, 'SM');
        
        s.data.NSub = cds.NSub;
        s.data.NUnstable = cds.NUnstable;
        
        print_diag(1,'Lyapunov Coefficient = %e\n', s.data.l1);
        s.msg = 'Hopf point';
        
    case 3 % LP
        
        if ~isempty(cds.s)
            s.data.borders.v = cds.s.borders.v;
            s.data.borders.w = cds.s.borders.w;
            s.data.A0        = cds.s.A0;
            s.data.Q0        = cds.s.Q0;
            s.data.T0        = cds.s.T0;
            s.data.kapa      = cds.s.kapa;
            cds.s = [];
            
            A0 = s.data.A0;
        else
            A0 = contjac(x);
            A0 = A0(:,1:end-1);
        end
        
        % Normal form coefficient
        [x0,p] = rearr(x); p = num2cell(p);
        s.data.a = nf_LP_L(A0, x0, p);
        
        % Eigenvalues
        NumEvals = cds.NSub + contopts.CIS_NExtra;
        s.data.eigenvals = eigs(A0, NumEvals, 'SM');
        
        s.data.NSub = cds.NSub;
        s.data.NUnstable = cds.NUnstable;
        
        print_diag(1,'Quadratic coefficient = %e\n', s.data.a);
        s.msg = 'Limit point';
end

NumEvals = cds.NSub+contopts.CIS_NExtra;
ndim = cds.ndim;

J = contjac(x);
if ~issparse(J)
    [d,evl1]= eig(J(:,1:ndim-1)');
else
    opt.disp= 0;
    [d,evl1]= eigs(J(:,1:ndim-1),NumEvals, 'SM', opt);
end

if isempty(s.data.eigenvals)    %%MF 9/3/11
    evals = diag(evl1);
else
    evals = s.data.eigenvals;
end                              %%MF 9/3/11

[eval_re,indx]=sort(-real(evals));
evals = evals(indx);

eval_r = evals(1:cds.NSub);
eval_l = evals(cds.NSub+1:end);

print_diag(2,'\nEigenvalues: (Subspace)\n')
for j=1:length(eval_r)
    print_diag(2,'%+e',real(eval_r(j)))
    if imag(eval_r(j))
        print_diag(2,' %+e i',imag(eval_r(j)));
    end
    print_diag(2,'\n')
end
print_diag(2,'Eigenvalues: (Reference)\n')
for j=1:min(length(eval_l),contopts.CIS_NExtra)
    print_diag(2,'%+e',real(eval_l(j)))
    if imag(eval_l(j))
        print_diag(2,' %+e i',imag(eval_l(j)));
    end
    print_diag(2,'\n')
end

cds.tfUpdate=1;
failed = 0;

%% ------------------------------------------------------------
function [S,L] = singmat

% 0: testfunction must vanish
% 1: testfunction must not vanish
% 2: testfunction must change   DV
% everything else: ignore this testfunction

S = [ 0 8 8 8      % DV
    8 0 2 8
    1 8 8 0 ];

L = [ 'BP'; 'H '; 'LP' ]; % rows of S

%% ------------------------------------------------------------
function [x,v] = locate(id, x1, v1, x2, v2)
global A1 Q1 T1

%print_diag(5,'In equilibriumL/locate\n');
tic
switch id
    case 1
        [x,v]     = locateBP(id, x1, x2, v1, v2);
    case 2
        [x,v]     = locateH3(x1, x2, v1, v2, A1, Q1, T1);
    case 3
        [x,v]     = locateLP(x1, x2, v1, v2, A1, Q1, T1);
    otherwise
        error('No locator defined for singularity %d', id);
end
print_diag(1,'Time spent in locator function: %f\n', toc);
%% ------------------------------------------------------------
function [x0,v0,options,failed] = init(varargin)
% DV: Currently not used, but would be good to implement in the future
global cds contopts
if nargin > 4
    if strcmp(varargin{5}.label,'BP')
        try
            [x0,v0,options] = init_BP_EPL(varargin{:});
            failed = 0;
        catch
            x0=[]; v0=[]; options=[];
            failed = 1;
        end
    else
        try
            [x0,v0,options] = init_EP_EPL(varargin{:});
            failed = 0;
        catch
            x0=[]; v0=[]; options=[];
            failed = 1;
        end
    end
else
    try
        [x0,v0,options] = init_EP_EPL(varargin{1:end-1});
        failed = 0;
        return;
    catch
        x0=[]; v0=[]; options=[];
        failed = 1;
    end
end
cds.num_sings = 1;
feval(cds.curve_options);

%% ------------------------------------------------------------
function varargout = done

% Placeholder

%% ----------------------------------------------------------
function [res,x,v] = adapt(x,v)
global cds
global A1 Q1 T1 evl1_r evl1_l
global A2 Q2 T2 evl2_r evl2_l
if cds.tfUpdate
    
    res = 1;                                   % updates tfs
    print_diag(3,'Adapt: updating update test functions\n');               % cds test
    
    [x0,p] = rearr(x); p = num2cell(p);
    A1 = ejac(x0,p);
    [T1, Q1, evl1_r, evl1_l] = contCIS_init(A1, 1);
    if isempty(T1)
        print_diag(3,'equilibriumL: contCIS_init failed to initialize Asub\n');
        res = 1;
        return;
    end
    
    A2 = A1;
    Q2 = Q1;
    T2 = T1;
    evl2_r = evl1_r;
    evl2_l = evl1_l;
    
    cds.tfUpdate = 0;
    
else
    res = 0; %no re-evaluations needed
end

% ---------------------------------------------------------------

%% ---------------------------------------------------------------
function [x,v] = locateBP(id, x1, x2, v1, v2)
%(MF) reliability: combine Newton with bisection
global cds contopts

%print_diag(5,'In equilibriumL/locateBP\n');

x1init = x1;
x2init = x2;

cds.s  = [];
NumEvals = cds.NSub+contopts.CIS_NExtra;
ndim   = cds.ndim;
Incr   = contopts.Cont_IncrFinDiff;
BordEr = contopts.EQ_SingLinSystTol;
MaxCorrIters = contopts.Locator_MaxIters;           % DV: Use new options
VarTolerance = contopts.Locator_VarTolerance;
FunTolerance = contopts.Locator_FunTolerance;

normdx = inf;
normf  = inf;

x  = 0.5*(x1+x2);
vx = v1 + v2;
vx = vx/norm(vx);

J = contjac(x);
if ~issparse(J)
    [V,evl] = eig([J;vx']);
    [d,evl1]= eig(J(:,1:ndim-1)');
    evl_vect = diag(evl);
    index = find(abs(evl_vect) == min(abs(evl_vect)));
    v = V(:,index(1));
else
    
    shift = 1e-2;
    B = [J;vx'];
    lB = size(B);
    R = rand(lB(1),1);
    R = R/norm(R);
    B = B - shift*speye(lB);
    v = bordCIS2(B,R,1);% v = B\R
    v = v/norm(v);
    
    if v ~= real(v)
        v = real(v);
    end
    opt.disp = 0;
    [d,evl1] = eigs(J(:,1:ndim-1), min(6,ndim-1), 'SM', opt);  % dsb %Question here about "6" eigenvalues, why this magic number?
    
end

bord.v = vx/norm(vx);
bord.w = v;

[y2,j] = min(abs(diag(evl1)));
bord.d = d(:,j);
i   = 0;

f1 = [zeros(ndim-1,2); eye(2)];
f3 = [zeros(ndim,1); 1];
vp  = cell(2);

bord.v0 = zeros(size(bord.v));
bord.w0 = zeros(size(bord.w));
i_improve = 0;
firstnewtonstep = 1;
while i <= MaxCorrIters
    
    %Step 1
    bord.w  = bord.w - (bord.w'*bord.v)*bord.v;
    bord.w  = bord.w/norm(bord.w);
    bord.W  = [bord.v bord.w];
    er_w    = norm(bord.w - bord.w0);
    er_v    = norm(bord.v - bord.v0);
    bord.v0 = bord.v;
    bord.w0 = bord.w;
    
    B = [J       bord.d
        bord.W' zeros(2,1)];
    
    if (i_improve < 1) && (er_w > BordEr || er_v > BordEr)
        % improve bord.d, bord.W
        i_improve = i_improve + 1;
        %Step 2
        [vg,fail] = bordCIS2(B,f1,2);  % vg = B\f1
        if fail
            % Newton seems to fail, use bisection
            print_diag(3,'Locate BP: Newton Fail: bordCIS2(B,f1,2) failed.\n')
            break
        end
        vp{1} = vg(1:ndim,1);
        vp{2} = vg(1:ndim,2);
        
        %Step 4
        [psig,fail] = bordCIS2(B',f3,2);  % psig = B'\f3
        if fail
            print_diag(3,'Locate BP: Newton Fail: bordCIS2(B(trans),f3,2) failed.\n')
            break
        end
        psi  = psig(1:ndim-1);
        
        %Step 8
        bord.d = psi/norm(psi);
        
        %Step 9
        bord.v = vp{1}/norm(vp{1});
        bord.w = vp{2};
        
    else
        %Step 2, Step 3
        f2  = -[feval(cds.curve_func, x); zeros(2,1)];
        f12 = [f1 f2];
        [vgh,fail] = bordCIS2(B,f12,2);  % vgh = B\f12
        if fail
            print_diag(3,'Locate BP: Newton Fail: bordCIS2(B,f12,2) failed.\n')
            break
        end
        
        %JH trying to avoid degeneracies
        nulldim = size(vgh);
        if nulldim(2) ~= 3
            print_diag(3,'Locate BP: Null Dimension ~= 3, Degenerate Case\n')
            x = []; v = [];
            return
        end
        
        vp{1} = vgh(1:ndim,1);
        vp{2} = vgh(1:ndim,2);
        xh   = vgh(1:ndim,3);
        muh  = vgh(ndim+1,3);
        
        %Step 4
        %tic%JH Testing NF Times
        [psig,fail] = bordCIS2(B',f3,2);  % psig = B'\f3
        if fail
            print_diag(3,'Locate BP: Newton Fail: bordCIS2(B(trans),f3,2) failed.\n')
            break;
        end
        psi  = psig(1:ndim-1);
        g    = psig(ndim:ndim+1);
        
        %Step 5
        dv{1} = Incr*vp{1}/norm(vp{1});
        dv{2} = Incr*vp{2}/norm(vp{2});
        dxh   = Incr*xh/norm(xh);
        
        for j = 1:2
            if norm(xh,inf) < eps
                eta(j,1) = 0;
            else
                eta(j,1) = norm(vp{j}) * norm(xh) / (Incr*Incr) * psi' * ...
                    (feval(cds.curve_func, x+dv{j}+dxh) - ...
                    feval(cds.curve_func, x+dv{j}) -     ...
                    feval(cds.curve_func, x+dxh) +       ...
                    feval(cds.curve_func, x));
            end
            
            M(j,j) = norm(vp{j})^2/(Incr*Incr)*psi'*...
                (  feval(cds.curve_func, x+dv{j}) - ...
                2*feval(cds.curve_func, x      ) + ...
                feval(cds.curve_func, x-dv{j}));
            
        end
        % new
        %print_diag(0,'NF time: ')%JH Testing NF Times
        
        M(1,2) = norm(vp{1})*norm(vp{2})/(Incr*Incr)*psi'*...
            (feval(cds.curve_func, x+dv{1}+dv{2}) - ...
            feval(cds.curve_func, x+dv{1}) - ...
            feval(cds.curve_func, x+dv{2}) + ...
            feval(cds.curve_func, x));
        M(2,1) = M(1,2);
        %toc%JH Testing NF Times
        %Step 6
        ksi = M\(g - eta);
        
        %Step 7
        dx = xh + ksi(1)*vp{1} + ksi(2)*vp{2};
        
        %Step 8 (modified for line search -- DSB)
        bord.d = psi/norm(psi);
        
        normdx_old = normdx;
        normf_old  = normf;
        
        x = x + dx;
        
        %         if norm(x-x1) > cds.h || norm(x-x2) > cds.h
        %             print_diag(3,'Locate BP: Newton Produced Divergent Solution\n');
        %             x=[];  v=[];
        %             return
        %         end
        
        f = feval(cds.curve_func, x) + muh*bord.d;
        
        if ~isempty(cds.usernorm)
            normf  = normU(f,2);
            normdx = normU(dx,1);
        else
            normdx = norm(dx);
            normf  = norm(f);
        end
        
        print_diag(2,'Locate BP: Newton Step: norm(dx)=%6.9f norm(f1)=%6.9f\n', normdx, normf);
        
        %Original Stopping Condition
        
        %% normdx                                                %%%
        %%VarTolerance
        %%normf
        %%FunTolerance
        %%%pause                                              %%%
        
        
        if normdx < VarTolerance && normf < FunTolerance
            v = v1 + v2;
            v = v/norm(v);
            
            v1 = v;
            v2 = vp{2} - v1*(vp{2}'*v1);
            v2 = v2/norm(v2);
            
            switch contopts.EQ_BranchingMethod
                case 0  %Normal Form Coefficients
                    cds.s.v1 = v1;
                    cds.s.v2 = v2 - 0.5*(M(2,2)/M(1,2))*v1;
                    cds.s.v2 = cds.s.v2/norm(cds.s.v2);
                case 1  %Perpendicular Guess
                    cds.s.v1 = v1;
                    cds.s.v2 = v2;
                case 2  %Bisection Guess
                    [cds.s.v1,cds.s.v2] = BPswitchBisect(x,v1,v2,x1init,x2init);
            end
            if ~isempty(cds.s.v2)
                print_diag(2,'Branch Angle = %6.6f * pi\n', innerangle(cds.s.v1,cds.s.v2)/pi())
            end
            J = contjac(x);
            if ~issparse(J)
                [d,evl1]= eig(J(:,1:ndim-1)');
            else
                opt.disp = 0;
                [d,evl1] = eigs(J(:,1:ndim-1), NumEvals, 'SM', opt);  % dsb %Question here about "6" eigenvalues, why this magic number?
            end
            cds.s.eigenvals = diag(evl1);
            %             print_diag(1,'Eigenvalues at BP: \n')
            %             print_diag(1,'%+e \n',cds.s.eigenvals)
            
            return;
        end
        
        %JH: Newton converging slowly, exit and reduce stepsize
        if ((normf > normf_old/2) || isnan(normf) ) && ~firstnewtonstep
            print_diag(3,'Locate BP: Residual Converging Too Slowly\n');
            x=[]; v=[];
            return
        end
        
        firstnewtonstep = 0;
        %Step 9
        bord.v = vp{1}/norm(vp{1});
        bord.w = vp{2};
        
        J = contjac(x);
    end
    
    if i == MaxCorrIters
        print_diag(3,'Locate BP: Failed to find x in %d Newton iterations\n',MaxCorrIters);
        x=[];v=[];
        return;
    end
    
    i = i+1;
end

%% ------------------Hopf Locator----------------------------------------
function [x,v] = locateH3(x1, x2, v1, v2, A1, Q1, T1)
global cds contopts
%print_diag(5,'In equilibriumL/locateH3\n');
cds.s = [];

%JH: Bisection 9/1/06----------
%v = v1 + v2;
%v = v/norm(v);
normdx = inf;
normf  = inf;
%NBisec = contopts.NBisectionIters;
X1 = x1; V1 = v1; X2 = x2; V2 = v2;
distX1X2 = norm(X1 - X2);
%------------------------------

ndim = cds.ndim;
NSub = cds.NSub;
Incr = contopts.Cont_IncrFinDiff;
MaxCorrIters = contopts.Locator_MaxIters;           % DV: Use new options
VarTolerance = contopts.Locator_VarTolerance;
FunTolerance = contopts.Locator_FunTolerance;

JacMatrix = contjac(x1);
A2 = JacMatrix(1:ndim-1, 1:ndim-1);
% DV: T1 is computed at currpoint and not at x1, so recompute
[T1, pout] = contCIS_step(cds.h, A1, Q1, T1, A2, cds.h_min);
Q1 = pout.Q;
Q1_partial = Q1(:, 1:NSub);
A1 = A2;

Asub   = T1(1:cds.NSub,1:cds.NSub);
[V,D1] = eig(Asub);
d      = diag(D1);

% find a complex eigenvalue with smallest absolute real part and set omega
% = its imaginary part
I     = find(imag(d) ~= 0);
if isempty(I)
    print_diag(1, 'Locate Hopf: All eigenvalues in start point are real \n');
    x = []; v = [];
    return
end
Rmin  = min(abs(real(d(I))));
I1    = find(abs(real(d)) == Rmin);
lamda = real(d(I1(1)));
omega = imag(d(I1(1)));				% get omega of Asub

kapa  = lamda*lamda + omega*omega;

[W,D2]= eig(Asub');
d     = diag(D2);
I     = find(imag(d) ~= 0);
Rmin  = min(abs(real(d(I))));
I2    = find(abs(real(d)) == Rmin);

% Input
% make the bordered system nonsigular
[Q,R,E] = qr([real(V(:,I1(1))) imag(V(:,I1(1)))]);
borders.v =  [real(V(:,I1(1))) imag(V(:,I1(1)))];

[Q,R,E] = qr([real(W(:,I2(1))) imag(W(:,I2(1)))]);
borders.w =  [real(W(:,I2(1))) imag(W(:,I2(1)))];

x = [x1;kapa];   % the variables are (u,alpha,kapa)

%data structure for constructing minimally augmented system dynamically
two_index = zeros(6,2);

two_index(1,:) = [1 2];
two_index(2,:) = [1 3];
two_index(3,:) = [1 4];
two_index(4,:) = [2 3];
two_index(5,:) = [2 4];
two_index(6,:) = [3 4];

transferForm = [1,1;1,2;2,1;2,2];

for i = 1:MaxCorrIters
    
    RED = Asub*Asub+kapa*eye(NSub);
    %RED
    
    %Step 1
    Bord  = [RED borders.w; borders.v' zeros(2)];
    bunit = [zeros(NSub,2);eye(2)];
    
    %Step 2
    vext = Bord\bunit;
    v1   = vext(1:NSub,1);
    v2   = vext(1:NSub,2);
    
    %Step 3
    wext = Bord'\bunit;
    w1   = wext(1:NSub,1);
    w2   = wext(1:NSub,2);
    
    % DV: Compute true left eigenvector of A^2 + kapa*I
    RED2 = A1*A1+kapa*speye(length(A1));
    Bord2 = [RED2' Q1_partial*borders.v; borders.w'*Q1_partial' zeros(2)];
    bunit2 = [zeros(length(A1), 2); eye(2)];
    wext2 = bordCIS1(Bord2, bunit2, 2);
    Qw1   = wext2(1:end-2,1);
    Qw2   = wext2(1:end-2,2);
    
    %Step 4
    borders.v = [vext(1:NSub,1)/norm(vext(1:NSub,1)),vext(1:NSub,2)/norm(vext(1:NSub,2))];
    borders.w = [wext(1:NSub,1)/norm(wext(1:NSub,1)),wext(1:NSub,2)/norm(wext(1:NSub,2))];
    
    %Step 5    
    Qv1  = Q1_partial * v1;
    Qv2  = Q1_partial * v2;
    QCv1 = Q1_partial * Asub * v1;
    QCv2 = Q1_partial * Asub * v2;
    
    x1 = x(1:ndim);
    JacMatrix_x = JacMatrix(:,1:end-1);
    
    contCISjacDirDerSquare_x1_Q1_Asubv1 = DirDerJacobian(x1,Q1_partial*Asub*v1,JacMatrix_x);
    contCISjacDirDerSquare_x1_Q1_Asubv2 = DirDerJacobian(x1,Q1_partial*Asub*v2,JacMatrix_x);
    contCISjacDirDerSquare_x1_Q1_v1     = DirDerJacobian(x1,Q1_partial*v1     ,JacMatrix_x);
    contCISjacDirDerSquare_x1_Q1_v2     = DirDerJacobian(x1,Q1_partial*v2     ,JacMatrix_x);
    
    gx(1,:) = -Qw1'*contCISjacDirDerSquare_x1_Q1_Asubv1 - Qw1'*A1*contCISjacDirDerSquare_x1_Q1_v1;
    gx(2,:) = -Qw1'*contCISjacDirDerSquare_x1_Q1_Asubv2 - Qw1'*A1*contCISjacDirDerSquare_x1_Q1_v2;
    gx(3,:) = -Qw2'*contCISjacDirDerSquare_x1_Q1_Asubv1 - Qw2'*A1*contCISjacDirDerSquare_x1_Q1_v1;
    gx(4,:) = -Qw2'*contCISjacDirDerSquare_x1_Q1_Asubv2 - Qw2'*A1*contCISjacDirDerSquare_x1_Q1_v2;
    
    f_x1 = JacMatrix_x;
    x2   = x(1:ndim);
    x2(end) = x2(end) + Incr;
    f_x2 = contjac(x2);
    f_x2 = f_x2(1:end,1:ndim-1);
    
    gp(1,1) = -Qw1'*(f_x2*QCv1 - f_x1*QCv1)/(Incr) - w1'*Asub*Q1_partial'*(f_x2*Qv1 - f_x1*Qv1)/(Incr);
    gp(2,1) = -Qw1'*(f_x2*QCv2 - f_x1*QCv2)/(Incr) - w1'*Asub*Q1_partial'*(f_x2*Qv2 - f_x1*Qv2)/(Incr);
    gp(3,1) = -Qw2'*(f_x2*QCv1 - f_x1*QCv1)/(Incr) - w2'*Asub*Q1_partial'*(f_x2*Qv1 - f_x1*Qv1)/(Incr);
    gp(4,1) = -Qw2'*(f_x2*QCv2 - f_x1*QCv2)/(Incr) - w2'*Asub*Q1_partial'*(f_x2*Qv2 - f_x1*Qv2)/(Incr);
    
    gk(1,1) = -w1'*v1;
    gk(2,1) = -w1'*v2;
    gk(3,1) = -w2'*v1;
    gk(4,1) = -w2'*v2;
    
    %Step 6
    x1  = x(1:ndim);
    f_x = JacMatrix;
    
    if i == 1
        A = [f_x,zeros(ndim-1,1)];
        
        testMatrix = [A;gx gp gk];
        [Index1, Index2] = hopf_locator_construction(testMatrix, two_index);
    end;
    
    B   = [f_x,zeros(ndim-1,1); gx(Index1,:),gp(Index1,:),gk(Index1,:); gx(Index2,:),gp(Index2,:),gk(Index2,:)];
    f1  = [feval(cds.curve_func, x1); vext(NSub+transferForm(Index1,1),transferForm(Index1,2)); vext(NSub+transferForm(Index2,1),transferForm(Index2,2))];
    
    %Step 7
    dx = bordCIS1(B,f1,2);  % dx = B\f1
    
    x = x - dx;
    
    % DV: Seems we do not need this either
    %      if (norm(x(1:ndim) - X1) > distX1X2) & (norm(x(1:ndim) - X2) > distX1X2)
    %          print_diag(3,'Locate Hopf, Newton Produced Divergent Solution\n');
    %          x=[];        v=[];
    %          return
    %      end
    
    kapa = x(end);	% kapa is considered as variable now
    
    % user defined norm
    normdx_old = normdx;
    normf_old  = normf;
    if ~isempty(cds.usernorm)
        normdx = sqrt(normU(dx(1:end-1),1)^2+dx(end)^2);
        normf  = normU(f1,2);
    else
        normdx = norm(dx);
        normf  = norm(f1);
    end
    print_diag(2,'Locate Hopf: Newton Step: norm(dx)=%6.9f norm(f1)=%6.9f\n', normdx, normf);
    
    if normdx < VarTolerance && normf < FunTolerance
        % End Locator Function, accept point.
        break;
    end
    
    %JH: Newton converging slowly, exit and reduce stepsize
    %DV: Exit and reduce stepsize does note work like that anymore
    %DV: Try just once and do not worry about monotonous convergence
    %     if (normf > 0.9*normf_old) || isnan(normf)
    %         print_diag(3,'Locate Hopf: Residual Converging Too Slowly\n');
    %         x=[];        v=[];
    %         return
    %     end
    
    if i == MaxCorrIters
        %when fail to locate Hopf point, report message and return(FQ, 08/01/2006)
        x = [];     v = [];
        print_diag(3,'locateHopf: Failed to find x in %d Newton iterations\n',MaxCorrIters);
        return;
    end
    
    %Step 8 Update f_x, Asub
    x_temp  = x(1:ndim);
    f_x     = contjac(x_temp);
    A2      = f_x(:,1:ndim-1);
    JacMatrix = f_x;
    [T,pout] = contCIS_step(cds.h, A1, Q1, T1, A2, cds.h_min);
    if isempty(T)
        print_diag(3,'Locate Hopf: contCIS_step failed to find Asub');
    else
        A1   = A2;
        Q1   = pout.Q;
        Q1_partial = Q1(:, 1:NSub);
        T1   = T;
        Asub = T(1:NSub,1:NSub);
    end
end

kapa = x(end);
x = x(1:ndim);
v = V1+V2;
v = v/norm(v);

cds.s.borders.v = borders.v;
cds.s.borders.w = borders.w;
cds.s.A0 = A1;
cds.s.Q0 = Q1;
cds.s.T0 = T1;
cds.s.kapa = kapa;


%% ---------------Limit Point Locator------------------------------

function [x,v] = locateLP(x1, x2, v1, v2, A1, Q1, T1)
global cds contopts

%print_diag(5,'In equilibriumL/locateLP\n');
cds.s = [];

ndim = cds.ndim;
NSub = cds.NSub;
Incr = contopts.Cont_IncrFinDiff;
MaxCorrIters = contopts.Locator_MaxIters;  
VarTolerance = contopts.Locator_VarTolerance;
FunTolerance = contopts.Locator_FunTolerance;
sparse_flag = contopts.CIS_SparseSolvers;            %MP 4-2016   

f_x = contjac(x1);
A0  = f_x(1:ndim-1, 1:ndim-1);
% DV: T1 is computed at currpoint and not at x1, so recompute
if ~isequal(A0, A1)
    [T0, pout] = contCIS_step(cds.h, A1, Q1, T1, A0, cds.h_min);
    Q0 = pout.Q;
else
    T0 = T1;
    Q0 = Q1;
end

Q1_partial = Q0(:,1:NSub);
Asub = T0(1:NSub,1:NSub);
d    = eig(Asub);

% find a real eigenvalue with smallest absolute value and call it Rmin
I    = find(imag(d) == 0);
Rmin = min(abs((d(I))));
I1   = find(abs((d)) == Rmin);
RED  = Asub - d(I1(1))*eye(NSub);

% Input
borders.v = null(RED);
borders.w = null(RED');

x = x1;

normdx = inf;
normf  = inf;

X1 = x1; V1 = v1; X2 = x2; V2 = v2;
distX1X2 = norm(X1 - X2);

%------------------------------
for i = 1:MaxCorrIters
    %Step 1
    Bord  = [Asub borders.w; borders.v' 0];
    bunit = [zeros(NSub,1);1];
    
    %Step 2
    vext  = Bord\bunit;
    
    %Step 3
    wext  = Bord'\bunit;
    Q1_partial = Q1(:,1:NSub); %MP 4-2016 
    
    % DV: To compute g_x we need the left eigenvector, Q*w is not good
    
       % Bord2 = [A0' Q1*borders.v; borders.w'*Q1' 0]; orig %MP 4-2016 
    Bord2 = [A0' Q1_partial*borders.v; borders.w'*Q1_partial' 0];  %MP 4-2016 
    bunit2 = [zeros(length(A0), 1); 1];
    wext2 = bordCIS1(Bord2, bunit2, 1);
    Qw = wext2(1:end-1);
    
    %Step 4
    borders.v = vext(1:NSub)/norm(vext(1:NSub));
    borders.w = wext(1:NSub)/norm(wext(1:NSub));
    
    %Step 5
    Qv = Q1_partial * vext(1:NSub);
    
    x1 = x;
    x1(1:ndim-1) = x1(1:ndim-1) - Incr/norm(Qv)*Qv;
    f_x1 = contjac(x1);
    x2 = x;
    x2(1:ndim-1) = x2(1:ndim-1) + Incr/norm(Qv)*Qv;
    f_x2 = contjac(x2);
    
    g_x = -(Qw'*f_x2 - Qw'*f_x1) * norm(Qv)/(2*Incr); % DV: minus sign was missing
    
    %Step 6
    B   = [f_x; g_x];
        
   %f1  = [feval(cds.curve_func, x); vext(end)]; %orig %MP 4-2016
   f1_0 = feval(cds.curve_func, x);              %MP 4-2016
   if ~sparse_flag                               %MP 4-2016
      f1 = [f1_0 vext(end)];                     %MP 4-2016
      f1=f1';                                    %MP 4-2016  
   elseif sparse_flag                            %MP 4-2016
      f1 = [f1_0; vext(end)];                    %MP 4-2016
   end                                           %MP 4-2016
    
    %Step 7
    dx = bordCIS2(B,f1,1);  % dx = B\f1
    
    x = x - dx;
    
    if (norm(x(1:ndim) - X1) > distX1X2) & (norm(x(1:ndim) - X2) > distX1X2)
        print_diag(3,'Locate Hopf, Newton Produced Divergent Solution\n');
        x=[];        v=[];
        return
    end
    
    % user defined norm
    normdx_old = normdx;
    normf_old  = normf;
    if ~isempty(cds.usernorm)
        normdx = normU(dx,1);
        normf  = normU(f1,2);
    else
        normdx = norm(dx);
        normf  = norm(f1);
    end
    
    print_diag(2,'Locate LP, Newton Step: norm(dx)=%6.9f norm(f1)=%6.9f\n', normdx, normf);
    if normdx < VarTolerance && normf < FunTolerance
        break;
    end
    
    %JH: Newton converging slowly, exit and reduce stepsize
    if ((normf > normf_old) || isnan(normf))
        print_diag(3,'Locate LP: Residual Converging Too Slowly\n');
        x=[];       v = [];
        return
    end
    if i == MaxCorrIters
        print_diag(3,'Locate LP: Failed to find x with %d Newton iterations\n',MaxCorrIters);
        x = [];     v = [];
        return;
    end
    
    %Step 8                 Update f_x, Asub
    f_x = contjac(x);
    A2 = f_x(:,1:ndim-1);
    [T,pout] = contCIS_step(cds.h, A0, Q0, T0, A2, cds.h_min);
    if isempty(T)                                      % CIS failed
        print_diag(3,'Locate LP: contCIS_step failed to find Asub\n');
    else
        Q1   = pout.Q;
        Q1_partial = Q1(:,1:NSub);
        Asub = T(1:NSub,1:NSub);
    end
    
end

v = v1 + v2;
v = v/norm(v);


%% ---------------------------------------------------------
function WorkspaceInit(x,v)

% Placeholder

%% ------CIS First Point-------------------------------------------
function [failed] = first_point(x,v)
%%global cds                %MP 7-2017 orig
global cds contopts         %MP 7-2017
global A1 T1 Q1 evl1_r evl1_l
sparse_flag = contopts.CIS_SparseSolvers;            %MP 7-2017  

%print_diag(5,'In equilibriumL/first_point\n');

[x0,p] = rearr(x); p = num2cell(p);
% A1 = ejac(x0,p);                     %MP 7-2017
[x0,p] = rearr(x); p = num2cell(p);
   if ~sparse_flag                               %MP 7-2017
      A1 = full(ejac(x0,p));                     %MP 7-2017  
   elseif sparse_flag                            %MP 7-2017
      A1 = ejac(x0,p);                           %MP 7-2017
   end                                           %MP 7-2017

[T1, Q1, evl1_r, evl1_l] = contCIS_init(A1, 0);
if isempty(T1)
    print_diag(3,'Equilibrium: contCIS_init failed to initialize Asub\n');
    failed = 1;
    return;
end

cds.s = [];

failed = 0;

%% ---------------------------------------------------------
function [failed,special_step] = candidate_step(x,v)
global cds
global A1 Q1 T1 evl1_r evl1_l
global A2 Q2 T2 evl2_r evl2_l

%print_diag(5,'In equilibriumL/candidate_step\n');

[x0,p] = rearr(x); p = num2cell(p);
A2 = ejac(x0,p);

[T,pout] = contCIS_step(cds.h, A1, Q1, T1, A2, ...
    cds.h_min, evl1_r, evl1_l);

if isempty(T)
    failed = 1;
    special_step = 0;
    print_diag(3,'EquilibriumL: contCIS_step failed to find Asub candidate_step\n');
    return;
end


evl2_r = pout.evl_r;
evl2_l = pout.evl_l;
Q2     = pout.Q;
T2     = T;

if pout.overlap
    cds.tfUpdate = 1;
end

% all done succesfully
failed = 0;
special_step = pout.overlap;

%% Accept Step
% ---------------------------------------------------------
function failed = accept_step(point)
global A1 Q1 T1 evl1_r evl1_l
global A2 Q2 T2 evl2_r evl2_l
global contopts

failed = 0;
%print_diag(5,'In equilibriumL/accept_step\n');
NExtra = contopts.CIS_NExtra;
space_angle = subspace(Q1,Q2)/pi();
print_diag(2,'Angle Between Subspaces: %+e * pi\n',space_angle);
A1 = A2; Q1 = Q2; T1 = T2;
evl1_r = evl2_r; evl1_l = evl2_l;
print_diag(2,'\nEigenvalues: (Subspace)\n')
for j=1:length(evl1_r)
    print_diag(2,'%+e',real(evl1_r(j)))
    if imag(evl1_r(j))
        print_diag(2,' %+e i',imag(evl1_r(j)));
    end
    print_diag(2,'\n')
end
print_diag(2,'Eigenvalues: (Reference)\n')
for j=1:min(length(evl1_l),NExtra)
    print_diag(2,'%+e',real(evl1_l(j)))
    if imag(evl1_l(j))
        print_diag(2,' %+e i',imag(evl1_l(j)));
    end
    print_diag(2,'\n')
end

%% Branch Angle Bisection
% ---------------------------------------------------------
function [V1,V2] = BPswitchBisect(x,v1,v2,x1init,x2init)

global cds
h0 = cds.h;
%P = [v1 v2]*[v1'; v2'];

X1 = x1init;
X2 = x2init;

% calculate alpha1, alpha2 and define thetaMin = max(alpha1,alpha2)
alpha1   = innerangle(v1, [v1 v2]*([v1'; v2']*(x - X1)));
alpha2   = innerangle(v1, [v1 v2]*([v1'; v2']*(x - X2)));
thetaMin = max(alpha1,alpha2)+1e-6;
if thetaMin < (pi()/2 - 1e-6) % Tolerance
    gamma = atan(thetaMin);
else
    print_diag(3,'BPswitchBisect: Cannot locate bracketing points x1, x2 on initial curve\n')
    V1 = [];   V2 = [];
    return;
end

vright = (gamma*v2 + v1)/(sqrt(gamma^2 +1));
vleft  = (gamma*v2 - v1)/(sqrt(gamma^2 +1));
vguess = v2;

%JH:    Determines the current direction of the tangent vector
if abs(vguess(end)) < 1e-6
    Vdir = sign(sum(vguess(1:end-1)));
else
    Vdir = sign(vguess(end));
end
%JH: Ensures the tangent vector is calculated in the 'positive' direction.
if Vdir < 0
    vguess = -vguess;
end

h1 = h0;
%Tangent vector bisection calculation
while 1
    if h1 < cds.h_min/2
        print_diag(3,'BPswitchBisect: Stepsize too small\n')
        V1 = [];   V2 = [];
        return;
    end
    xguess = x + h1*vguess;
    %%%tol = 1e-4; %%% test BPswitvh #2
    %%%if abs(x(end)-xguess(end)) >= tol h1 = h1/2; continue; end;%%%
    %% x_end = x(end) %%%
    %%  xguess_end = xguess(end) %%%
    [xcorr,vcorr] = newtcorrL(xguess,vguess);
    if   isempty(xcorr)   h1 = h1/2; continue; end;
    
    theta = innerangle(v1,[v1 v2]*([v1'; v2']*(x-xcorr)));
    if theta < thetaMin && norm(vleft - vright) > 1e-6
        if v1'*[v1 v2]*([v1'; v2']*(x-xcorr)) > 0
            vright = vguess;
            print_diag(2,'Right Angle = %d\n', thetaMin - theta)%test
        else
            vleft  = vguess;
            print_diag(2,'Left  Angle = %d\n', thetaMin - theta)%test
        end
        
        vguess = vleft + vright;
        vguess = vguess/norm(vguess);
    else
        break;
    end
end
%Next point prediction.
h1 = h0;
while 1
    if h1 < cds.h_min/2
        print_diag(3,'BPswitchBisect: Cannot locate nearby point\n')
        V1 = [];   V2 = [];
        return;
    end
    xguess = x + h1*vcorr;
    [xnext,vnext] = newtcorrL(xguess,vcorr);
    if   isempty(xnext)
        h1 = h1/2;
    else
        theta = innerangle(v1,[v1 v2]*([v1'; v2']*(x-xnext)));
        if theta > thetaMin     %check for fallback or curves too close
            cds.s.xnext = xnext;
            cds.s.vnext = vnext;
            V1 = v1;  V2 = vguess; V2 = V2/norm(V2);
            return;
        else
            print_diag(3,'BPswitchBisect: Unresolved Fallback or Branches too close, Unable to switch branches.\n')
            V1 = v1;   V2 = [];
            return;
        end
        
    end
end