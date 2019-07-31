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

out{14} = @CIS_first_point;
out{15} = @CIS_step;

%% -------------------------------------------------------
function func = curve_func(X, ~) % unused argument is CISdata

global cds
[x,p] = rearr(X); p = num2cell(p);
func = feval(cds.func,0, x, p{:});

%% ---------------------------------------------------------------
function jac = jacobian(X, ~) % unused argument is CISdata

[x0,p] = rearr(X); p = num2cell(p);
jac = [ejac(x0,p) ejacp(x0,p)];

%% ---------------------------------------------------------------
function point = defaultprocessor(varargin)

global cds
point = varargin{1};
if nargin > 1 && cds.i ~= 1
    cds.tfUpdate = 1;     % Need to update test functions
end
savePoint(varargin{:});

%% -------------------------------------------------------------
function options                % DV
% DV:  change options in the global variable if necessary
% DV: Do not use input and output arguments anymore

%% ---------------------------------------------------------------
function [out, failed] = testf(id, X, v, CISdata)
print_diag(5,'In equilibriumL/testf\n');

evl_r  = CISdata.evl_r;
NSub   = CISdata.NSub;

out(4) = 0;
failed = 0;

for i = id
    lastwarn('');
    switch i
        case 1 % BP
            [x0, p] = rearr(X); p = num2cell(p);
            A = ejac(x0, p);
            b = ejacp(x0, p);
            c = v(1:end-1);
            d = v(end);
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
            if NSub == 1
                out(2) = 1;
            else    
                out(2) = mu_min*(-1)^M_u;     
            end   
        case 3 % H2
            out(3) = sum(real(evl_r >= 0));  % DV
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
print_diag(5,'In equilibriumL/userf\n');

global cds

failed = 0;

if id == 0 % evaluate all
    id = 1:cds.nUserf;
end

out(cds.nUserf) = 0;
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
print_diag(5,'In equilibriumL/process\n');

if isfield(point, 's'); s.data = point.s; else; s.data = []; end

s.data.CISdata = point.CISdata;

[x0,p] = rearr(point.x); p = num2cell(p);

switch id
    case 1 % BP
        s.msg     = 'Branch point';
    case 2 % H
        % Normal form coefficient
        if isfield(s.data, 'kapa') && isfield(s.data, 'borders')
            [l1, kapa] = nf_H_L(point.CISdata, x0, p, point.s.kapa, point.s.borders);
        else
            [l1, kapa] = nf_H_L(point.CISdata, x0, p);
        end
        s.data.l1 = l1;
        s.data.kapa = kapa;
        
        print_diag(1,'Lyapunov Coefficient = %e\n', s.data.l1);
        s.msg = 'Hopf point';
        
    case 3 % LP
        % Normal form coefficient
        if isfield(s.data, 'borders')
            a = nf_LP_L(point.CISdata, x0, p, s.data.borders);
        else
            a = nf_LP_L(point.CISdata, x0, p);
        end
        s.data.a = a;
        
        print_diag(1,'Quadratic coefficient = %e\n', s.data.a);
        s.msg = 'Limit point';
end

% Eigenvalues
s.data.eigenvals = [point.CISdata.evl_r; point.CISdata.evl_l];

eval_r = point.CISdata.evl_r;
eval_l = point.CISdata.evl_l;

in_log_file(eval_r,eval_l);

failed = 0;

%% ------------------------------------------------------------
function [S,L] = singmat
%print_diag(5,'In equilibriumL/singmat\n');

% 0: testfunction must vanish
% 1: testfunction must not vanish
% 2: testfunction must change   DV
% everything else: ignore this testfunction

S = [ 0 8 8 8      % DV
    8 0 2 8
    1 8 8 0 ];

L = [ 'BP'; 'H '; 'LP' ]; % rows of S

%% ------------------------------------------------------------
function [ps] = locate(id, p1, p2)
print_diag(5,'In equilibriumL/locate\n');
    
tic
switch id
    case 1
        ps     = EP_locateBP(p1, p2);
    case 2
        ps     = EP_locateH3(p1, p2);
    case 3
        ps     = EP_locateLP(p1, p2);
    otherwise
        error('No locator defined for singularity %d', id);
end
print_diag(1,'Time spent in locator function: %f\n', toc);

%% ------------------------------------------------------------
function varargout = init(~,~) % unused arguments are x and v

% Placeholder
varargout{1} = 0;

%% ------------------------------------------------------------
function varargout = done

% Placeholder
varargout{1} = 0;

%% ----------------------------------------------------------
function [res,X,V,CISdata] = adapt(X,V,CISdata,tfUpdate)

if tfUpdate
    res = 1;                                   % updates tfs
    print_diag(3,'Adapt: updating test functions\n');    % cds test
    
    NSub = CISdata.NSub;
    NUnstable = CISdata.NUnstable;
    [x,p] = rearr(X); p = num2cell(p); A = ejac(x,p);    
    CISdata = contCIS_init(A, 1, NSub, NUnstable); % recompute subspace
    
    if isempty(CISdata)
        print_diag(3,'equilibriumL: contCIS_init failed to initialize Asub\n');
        res = 1;
        return;
    end
    
else
    res = 0; %no re-evaluations needed
end

%% ------CIS First Point-------------------------------------------
function CISdata = CIS_first_point(X)

global contopts
NSub = contopts.CIS_NSub;            
NUnstable = contopts.CIS_NUnstable;
[x,p] = rearr(X); p = num2cell(p); A = ejac(x, p);
   
CISdata = contCIS_init(A, 0, NSub, NUnstable);

% Eigenvalues                            % MP 9/2018
s.data.eigenvals = [CISdata.evl_r; CISdata.evl_l];  % MP

eval_r = CISdata.evl_r;     % MP
eval_l = CISdata.evl_l;     % MP  
in_log_file(eval_r,eval_l); % MP

%% ---------------------------------------------------------
function CISdata = CIS_step(X, CISdata1)

[x,p] = rearr(X); p = num2cell(p); A = ejac(x, p);
CISdata = contCIS_step(A, CISdata1);

% Eigenvalues                            % MP 9/2018
if ~ isempty(CISdata)
  s.data.eigenvals = [CISdata.evl_r; CISdata.evl_l];

  eval_r = CISdata.evl_r;
  eval_l = CISdata.evl_l;  
  in_log_file(eval_r,eval_l); % MP
end
