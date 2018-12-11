function out = bogdanovtakensL
%
% Equilibrium curve definition file for a problem in odefile
%
out{1}  = @curve_func;
out{2}  = @defaultprocessor;
out{3}  = @options;
out{4}  = @jacobian;
out{5}  = @hessians;
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

return

%----------------------------------------------------
function func = curve_func(X, CISdata)
%print_diag(5,'In bogdanovtakensL/curve_func\n');

global cds

[x,p] = rearr(X); p     = num2cell(p);

NSub = CISdata.NSub;
T    = CISdata.T;
Asub = T(1:NSub,1:NSub);

% g
Bord    = [Asub cds.borders.w; cds.borders.v' 0];
bunit   = [zeros(NSub,1); 1];
vext    = Bord\bunit;

% g_lambda     
b       = vext;     % DV
b(end)  = 0;        % DV
vext2   = Bord\b;   % DV

func    = [feval(cds.func, 0, x, p{:}) ; vext(end); vext2(end)];

%---------------------------------------------------
function jac = jacobian(X, CISdata)

%print_diag(5,'In bogdanovtakensL/jacobian\n');

global cds contopts

Incr = contopts.Increment;

[x,p] = rearr(X); p = num2cell(p);

f_x = [CISdata.A lpjacp(x,p)];

NSub  = CISdata.NSub;
T     = CISdata.T;
A     = CISdata.A;
Asub  = T(1:NSub,1:NSub);
Q     = CISdata.Q(:, 1:NSub);

%Step 1
Bord  = [Asub cds.borders.w; cds.borders.v' 0];
bunit = [zeros(NSub,1);1];

%Step 2
vext = Bord\bunit;
bv       = vext;     % DV
bv(end)  = 0;        % DV
vext2   = Bord\bv;   % DV

%Step 3
wext    = Bord'\bunit;

% DV: Use true left eigenvectors
Bord2 = [A' Q*cds.borders.v; cds.borders.w'*Q' 0];
bunit2 = [zeros(length(A), 1); 1];
wext2 = Bord2 \ bunit2; % bordCIS1(Bord2, bunit2, 1);
Qw1 = wext2(1:end-1);

bunit3 = [Qw1; 0];
wext3 = Bord2 \ bunit3; % bordCIS1(Bord2, bunit3, 1);
Qw2 = wext3(1:end-1);

%Step 4
cds.borders.v = vext(1:NSub)/norm(vext(1:NSub));
cds.borders.w = wext(1:NSub)/norm(wext(1:NSub));

%Step 5
Qv = Q * vext(1:NSub);
% Qw = Q3 * wext(1:NSub);

x1 = X;
x1(1:cds.ncoo) = x1(1:cds.ncoo) - Incr/norm(Qv)*Qv;
[xt,p] = rearr(x1); p = num2cell(p);
f_x1 = [lpjac(xt,p) lpjacp(xt,p)];

x2 = X;
x2(1:cds.ncoo) = x2(1:cds.ncoo) + Incr/norm(Qv)*Qv;
[xt,p] = rearr(x2); p = num2cell(p);
f_x2 = [lpjac(xt,p) lpjacp(xt,p)];

g_x = -Qw1'* (f_x2 - f_x1) * norm(Qv)/(2*Incr);

%step 6
Qv2 = Q * vext2(1:NSub);                               % DV

x1 = X;                                                % DV
x1(1:cds.ncoo) = x1(1:cds.ncoo) - Incr/norm(Qv2)*Qv2;   % DV
[xt,p] = rearr(x1); p = num2cell(p);                    % DV
f_x3 = [lpjac(xt,p) lpjacp(xt,p)];                      % DV

x2 = X;                                                % DV
x2(1:cds.ncoo) = x2(1:cds.ncoo) + Incr/norm(Qv2)*Qv2;   % DV
[xt,p] = rearr(x2); p = num2cell(p);                    % DV
f_x4 = [lpjac(xt,p) lpjacp(xt,p)];                      % DV

g2_x = -Qw2'* (f_x2 - f_x1) * norm(Qv)/(2*Incr) -Qw1'* (f_x4 - f_x3) * norm(Qv2)/(2*Incr);        % DV

%Step 7
jac = [f_x; g_x; g2_x];

%---------------------------------------------------
function hess = hessians(varargin)

%print_diag(5,'In bogdanovtakensL/hessians\n');
hess =[];
%---------------------------------------------------
function point = defaultprocessor(varargin)
global cds
point = varargin{1};
if nargin > 1 && cds.i ~= 1
    cds.tfUpdate = 1;     % Need to update test functions
end
failed = savePoint(varargin{:});

%----------------------------------------------------
function options

%print_diag(5,'In bogdanovtakensL/options\n');
global contopts
contopts = contset(contopts, 'Locators', [0 0 0]);

% -------------------------------------------------------
function [out, failed] = testf(id, x, v, CISdata)

%print_diag(5,'In bogdanovtakensL/testf\n');

global cds              

[x0,p] = rearr(x); p = num2cell(p);

NSub   = CISdata.NSub;
Asub   = CISdata.T(1:NSub,1:NSub);
evl3_r = CISdata.evl_r;
evl3_l = CISdata.evl_l;

failed = 0;

if any(ismember(id, [1, 2]))
    c = nf_BT_L(CISdata, x0, p);
end

if any(ismember(id, [3, 4, 5]))
    [~, I]   = sort(abs(evl3_r));
    E = evl3_r(I);
    [~, ind] = sort(abs(E));
    E = E(ind);
    E = E(3:end);
    n_unstable = sum(real(E) > 0);
end

for i=id
    lastwarn('');
    switch i
        case 1 % a2 normal form coefficient (BT-CP)
            out(1) = c(1);
        case 2 % b2 normal form coefficient
            out(2) = c(2);
        case 3 % triple zero
            lambda_min = min(abs(E));
            out(3) = lambda_min * (-1)^n_unstable;
        case 4
            NSub2 = NSub - 2;
            bialt_eigs  = zeros(NSub2*(NSub2-1)/2,1);
            kk=1;
            for ii = 1:(NSub2-1)
                for jj = (ii+1):NSub2
                    bialt_eigs(kk) = E(ii)+E(jj);
                    kk = kk + 1;
                end
            end
            mu_min = min(abs(bialt_eigs));
            M_u  = sum(real(bialt_eigs)>0);
            out(4) = mu_min*(-1)^M_u;
        case 5 
            out(5) = n_unstable;
        otherwise % BP
            error('Internal error: Requested testfunction not available')
    end
    if ~isempty(lastwarn)
        msg = sprintf('Could not evaluate tf %d\n', i);
        failed = 1;
    end
    
end
%%%%%fprintf('%+e %+e %+e \n',out)%DEBUG
%------------------------------------------------------
function [out, failed] = userf(id, UserInfo, x)
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
%---------------------------------------------------------
function [failed,s] = process(id, point, s)

%print_diag(5,'In bogdanovtakensL/process\n');
global cds

if isfield(point, 's'); s.data = point.s; else; s.data = []; end

s.data.CISdata = point.CISdata;

[x0,p] = rearr(point.x); p = num2cell(p);

switch id
    case 1 % BTa2
        [~, b4] = nf_BT4_L(point.CISdata, x0, p);
        print_diag(1,'BTCP normal form coefficients (a2, b4) = (%e, %e)\n', point.tvals(1), b4);
        s.msg  = 'Bogdanov-Takens with vanishing normal form coefficient a2';
    case 2 % BTb2
        [a3, b3] = nf_BT3_L(point.CISdata, x0, p);
        [a4, ~ ] = nf_BT4_L(point.CISdata, x0, p);
        print_diag(1,'BTGH normal form coefficients (b2, a3, b3, a4) = (%e, %e, %e, %e)\n', point.tvals(2), a3, b3, a4);
        s.msg  = 'Bogdanov-Takens with vanising normal form coefficient b2';
    case 3 % BT3
        print_diag(1,'Triple zero. No normal form coefficients available');
        s.msg = sprintf('Triple zero');
    case 4 % BTH
        print_diag(1,'BT + Hopf. No normal form coefficients available');
        s.msg = sprintf('Bogdanov-Takens with Hopf');
    otherwise
        error('singularity not supported')
end

failed = 0;

%--------------------------------------------------------
function  [S,L] = singmat

%print_diag(5,'In bogdanovtakensL/singmat\n');

% 0: testfunction must vanish
% 1: testfunction must not vanish
% everything else: ignore this testfunction

S = [  0 8 1 1 8  
       8 0 1 1 8
       8 8 0 8 8
       8 8 8 0 2 ];

L = [ 'BTa2'; 'BTb2'; 'BT3 '; 'BTH '];        % rows of S

%------------------------------------------------------
function [ps] = locate(id, p1, p2)
% DV: There are no locators available at the moment
error('No locator defined for singularity %d', id);

%------------------------------------------------------
function [x0,v0,options,failed] = init(varargin)
% DV: Never called

%--------------------------------------------------------
function varargout = done

%---------------------------------------------------------
function [res,x,v,CISdata] = adapt(x,v,CISdata,tfUpdate)

%print_diag(5,'In bogdanovtakensL/adapt\n');

global cds

if tfUpdate && ~isempty(CISdata)                                  
    res = 1;                                   % updates tfs
    print_diag(3,'Adapt: updating test functions\n');    % cds test
    
    NSub = CISdata.NSub;
    NUnstable = CISdata.NUnstable;
    [x0,p] = rearr(x); p = num2cell(p); A = ejac(x0,p);    
    CISdata = contCIS_init(A, 1, NSub, NUnstable); % recompute subspace
    
    if isempty(CISdata)
        print_diag(3,'limitpointL: contCIS_init failed to initialize Asub\n');
        res = 1;
        return;
    end
    
    cds.borders = find_newborders_LP(CISdata);
    
else
    res = 0; % no re-evaluations needed
end

% ------------------------------------------------------

function WorkspaceDone

% -------------------------------------------------------


% ---------------------------------------------------------
function CISdata = CIS_first_point(X)

global contopts cds
NSub = contopts.CIS_NSub;            
NUnstable = contopts.CIS_NUnstable;
[x,p] = rearr(X); p = num2cell(p); A = lpjac(x, p);
   
CISdata = contCIS_init(A, 0, NSub, NUnstable);

if ~isempty(CISdata)
    cds.borders = find_newborders_LP(CISdata);
end

%% ---------------------------------------------------------
function CISdata = CIS_step(X, CISdata1)

[x,p] = rearr(X); p = num2cell(p); A = lpjac(x, p);
CISdata = contCIS_step(A, CISdata1);