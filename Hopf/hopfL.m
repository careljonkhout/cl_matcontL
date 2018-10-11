function out = hopfL
%
% Hopf curve definition file for a problem in odefile

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

%-------------------------------------------------------
function func = curve_func(X, CISdata)
%print_diag(5,'in hopfL/curve_func\n');

[x,p,k] = rearr(X); p = num2cell(p);

global cds

NSub = CISdata.NSub;
Asub = CISdata.T(1:NSub,1:NSub);

Bord = [Asub*Asub+k*eye(NSub) cds.borders.w;cds.borders.v' zeros(2)];
bunit= [zeros(NSub,2); eye(2)];

vext  = Bord\bunit;
vext1 = vext(NSub+cds.index1(1),cds.index1(2));
vext2 = vext(NSub+cds.index2(1),cds.index2(2));
func  = [feval(cds.func, 0, x, p{:}); vext1;vext2];

%--------------------------------------------------------
function jac = jacobian(X, CISdata)
%print_diag(5,'in jacobian of hopfL\n');

global cds contopts

nap = cds.nap;
Incr = contopts.Increment;
[x,p_temp,k] = rearr(X); p = num2cell(p_temp);

NSub = CISdata.NSub;
A2   = CISdata.A;
Asub = CISdata.T(1:NSub,1:NSub);
Q3   = CISdata.Q(:, 1:NSub);

RED   = Asub*Asub+k*eye(NSub);
Bord  = [RED cds.borders.w;cds.borders.v' zeros(2)];
bunit = [zeros(NSub,2);eye(2)];
vext  = Bord\bunit;
wext  = Bord'\bunit;
J     = hjac(x,p);
JP    = hjacp(x,p);

% DV: Compute true left eigenvector of A^2 + kapa*I
RED2 = A2*A2+k*speye(length(A2));
Bord2 = [RED2' Q3*cds.borders.v; cds.borders.w'*Q3' zeros(2)];
bunit2 = [zeros(length(A2), 2); eye(2)];
wext2 = bordCIS1(Bord2, bunit2, 2);

w1   = wext(1:NSub,cds.index1(1));
v1   = vext(1:NSub,cds.index1(2));
Qw1  = wext2(1:end-2,cds.index1(1));
Qv1  = Q3*v1;
QCv1 = Q3*Asub*v1;
w2   = wext(1:NSub,cds.index2(1));
v2   = vext(1:NSub,cds.index2(2));
Qw2  = wext2(1:end-2,cds.index2(1));
Qv2  = Q3*v2;
QCv2 = Q3*Asub*v2;

gx1 = -Qw1'*hjacDirDerSquare(x,p,Q3,Asub*v1) - Qw1'*A2*hjacDirDerSquare(x,p,Q3,v1);
gk1 = -w1'*v1;

gx2 = -Qw2'*hjacDirDerSquare(x,p,Q3,Asub*v2) - Qw2'*A2*hjacDirDerSquare(x,p,Q3,v2);
gk2 = -w2'*v2;

ap = cds.ActiveParams;
gp = zeros(2, nap);
for i=1:nap
    p1        = p_temp;
    p1(ap(i)) = p1(ap(i)) - Incr;
    p1        = num2cell(p1);
    f_x1      = hjac(x,p1);
    p2        = p_temp;
    p2(ap(i)) = p2(ap(i)) + Incr;
    p2        = num2cell(p2);
    f_x2      = hjac(x,p2);
    
    gp(1,i) = -Qw1'*(f_x2*QCv1 - f_x1*QCv1)/(2*Incr) - Qw1'*A2*(f_x2*Qv1 - f_x1*Qv1)/(2*Incr);
    gp(2,i) = -Qw2'*(f_x2*QCv2 - f_x1*QCv2)/(2*Incr) - Qw2'*A2*(f_x2*Qv2 - f_x1*Qv2)/(2*Incr);
end

jac = [J JP zeros(cds.ncoo,1); gx1,gp(1,:),gk1; gx2,gp(2,:),gk2];

%debug('exit jacobian\n');

%------------------------------------------------------
function hess = hessians(varargin)
hess =[];

%------------------------------------------------------
function point = defaultprocessor(varargin)
global cds
point = varargin{1};
if nargin > 1 && cds.i ~= 1
    cds.tfUpdate = 1;     % Need to update test functions
end
failed = savePoint(varargin{:});

%-------------------------------------------------------
function options
global contopts
contopts = contset(contopts, 'Locators', [0 0 0 0]);

% ---------------------------------------------------------------

function [out, failed] = testf(id, x, v, CISdata)
print_diag(5,'In HopfpointL/testf\n');

global cds

[x0,p,kapa] = rearr(x); p1 = num2cell(p);

NSub  = CISdata.NSub;
evl_r = CISdata.evl_r;

for i=id
    lastwarn('');
    
    switch i
        case 1 % BT
            out(1) = kapa;
        case 2 % ZH
            evl_r_min = min(abs((evl_r)));
            m_u       = sum((real(evl_r)>0));
            out(2)    = evl_r_min*(-1)^m_u;
        case 3 % DH1
            NSub2    = NSub-2;
            [~, I]   = sort(abs(real(evl_r)));
            E = evl_r(I);
            if isreal(E(1)) && ~isreal(E(2)) && ~isreal(E(3))
                E(3) = E(1);
            end
            E(1:2) = [];
            kk=1;
            bialt_eigs  = zeros(NSub2*(NSub2-1)/2,1);
            for ii = 1:(NSub2-1)
                for jj = (ii+1):NSub2
                    bialt_eigs(kk) = real(E(ii)+E(jj));
                    kk = kk + 1;
                end
            end
            M_u  = sum(bialt_eigs>0);
            mu_min = min(abs(bialt_eigs));
            out(3) = mu_min*(-1)^M_u;
        case 4  %GH
            if kapa > 0
                l1 = nf_H_L(CISdata, x0, p1, kapa, cds.borders);
                if isempty(l1) % failed
                    out = [];
                    return
                end
                out(4) = l1;
            else
                out(4) = 500;
            end
        case 5  %DH2
            [~, I]   = sort(abs(real(evl_r)));
            E = evl_r(I);
            if isreal(E(1)) && ~isreal(E(2)) && ~isreal(E(3))
                E(3) = E(1);
            end
            E(1:2) = [];
            out(5) = sum(real(E) > 0);  % DV: More robust test function
        otherwise
            error('No such testfunction');
    end
end
failed = 0;
%-------------------------------------------------------
function [out, failed] = userf(id, x, v)
print_diag(5,'In HopfpointL/userf\n');

global cds

failed = 0;

if id == 0 % evaluate all
    id = 1:cds.nUserf;
end

out(cds.nUserf) = 0;
for ii = id
    lastwarn('');
    [x0,p] = rearr(x); p = num2cell(p);
    
    p{cds.ActiveParams} = x(end);
    
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
% ---------------------------------------------------------------

function [failed,s] = process(id, point, s)
print_diag(5,'In HopfpointL/process\n');

global cds

if isfield(point, 's'); s.data = point.s; else; s.data = []; end

s.data.CISdata = point.CISdata;

[x0,p,kapa] = rearr(point.x); p = num2cell(p);

switch id
    case 1 % BT
        s.data.c = nf_BT_L(point.CISdata, x0, p);
        print_diag(1,'BT quadratic coefficients (a2, b2) = (%e, %e)\n', s.data.c(1), s.data.c(2));
        s.msg  = 'Bogdanov-Takens point';
    case 2 % ZH
        s.data.c = nf_ZH_L(point.CISdata, x0, p);
        print_diag(1,'ZH coefficients (B0, C0, E0) = (%e, %e, %e)\n', s.data.c(1), s.data.c(2), s.data.c(3));
        s.msg  = 'Zero-Hopf point';
    case 3 % DH
        s.data.c = nf_HH_L(point.CISdata, x0, p);
        print_diag(1,'HH coefficients (G2100, G1011, H1110, H0021) = (%e, %e, %e, %e)\n', ...
            s.data.c(1,1), s.data.c(1,2), s.data.c(2,1), s.data.c(2,2));
        s.msg  = 'Double Hopf point';
    case 4 % GH
        s.data.c = nf_GH_L(point.CISdata, x0, p, kapa, cds.borders);
        print_diag(1,'Second Lyapunov coefficient l2 = %e\n', s.data.c);
        s.msg  = 'Generalized Hopf point';
end

failed = 0;
%-------------------------------------------------------
function [S,L] = singmat
%print_diag(5,'In HopfpointL/singmat\n');

%elseif strcmp(arg, 'singmat')
% 0: testfunction must vanish
% 1: testfunction must not vanish
% everything else: ignore this testfunction

S = [  0 8 8 8 8
    1 0 8 8 8
    8 8 0 8 2
    8 1 8 0 8 ]; % DV:TEST

L = [ 'BT'; 'ZH'; 'DH';'GH'];

%-------------------------------------------------------
function ps = locate(id, p1, p2)
% DV: There are no locators available at the moment
error('No locator defined for singularity %d', id);

%-------------------------------------------------------
function [x0,v0,options,failed] = init(varargin)
% DV: Never called

%-------------------------------------------------------
function varargout = done

%------------------------------------------------------
function [res,x,v,CISdata] = adapt(x,v,CISdata, tfUpdate)
print_diag(5,'In HopfpointL/adapt\n');

res = 0;
global cds

if ~isempty(CISdata)
    if tfUpdate
        res = 1;                                   % updates tfs
        print_diag(3,'Adapt: updating test functions\n');    % cds test
        
        NSub = CISdata.NSub;
        NUnstable = CISdata.NUnstable;
        [x0,p] = rearr(x); p = num2cell(p); A = ejac(x0,p);
        CISdata = contCIS_init(A, 1, NSub, NUnstable); % recompute subspace
        
        if isempty(CISdata)
            print_diag(3,'hopfL: contCIS_init failed to initialize Asub\n');
            res = 1;
            return;
        end
        
        [cds.borders, cds.index1, cds.index2] = find_new_borders_H(CISdata,x0,p);
        
    else
        
        [x0,p,kapa] = rearr(x); p = num2cell(p);
        [cds.borders, cds.index1, cds.index2] = find_new_borders_H(CISdata, x0, p, kapa);
        
    end
end

% DV: We cannot continue a branch of neutral saddles when we are using CIS
kapa = x(end);
if kapa < 0
    cds.lastpointfound = 1;
    print_diag(0,'Hopf curve becomes neutral saddle: Continuation stopped.\n');
end


%% ------CIS First Point-------------------------------------------
function CISdata = CIS_first_point(X)

global contopts cds
NSub = contopts.CIS_NSub;
NUnstable = contopts.CIS_NUnstable;
[x,p] = rearr(X); p = num2cell(p); A = ejac(x, p);

CISdata = contCIS_init(A, 0, NSub, NUnstable);

[cds.borders, cds.index1, cds.index2] = find_new_borders_H(CISdata,x,p);

%% ---------------------------------------------------------
function CISdata = CIS_step(X, CISdata1)

[x,p] = rearr(X); p = num2cell(p); A = ejac(x, p);
CISdata = contCIS_step(A, CISdata1);