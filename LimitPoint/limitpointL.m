function out = limitpointL
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

print_diag(5,'In limitpointL/curve_func\n');

global cds

NSub   = CISdata.NSub;
Asub   = CISdata.T(1:NSub,1:NSub);

Bord    = [Asub cds.borders.w;cds.borders.v' 0];
bunit   = [zeros(NSub,1); 1];
vext    = Bord\bunit;

[x, p] = rearr(X); p = num2cell(p);
func    = [feval(cds.func, 0, x, p{:}) ; vext(end)];

%---------------------------------------------------
function jac = jacobian(X, CISdata)

print_diag(5,'In limitpointL/jacobian\n');

global cds contopts
Incr = contopts.Increment;

[x,p] = rearr(X); p = num2cell(p);
f_x = [CISdata.A lpjacp(x,p)];

NSub   = CISdata.NSub;
T      = CISdata.T;
Asub   = T(1:NSub,1:NSub);
Q3     = CISdata.Q(:, 1:NSub);

%Step 1
Bord  = [Asub cds.borders.w; cds.borders.v' 0];
bunit = [zeros(NSub,1);1];

%Step 2
vext = Bord\bunit;

%Step 3
wext = Bord'\bunit;

% DV: Use true left eigenvector
Bord2 = [CISdata.A' Q3*cds.borders.v; cds.borders.w'*Q3' 0];
bunit2 = [zeros(length(CISdata.A), 1); 1];
wext2 = bordCIS1(Bord2, bunit2, 1);

% %Step 4
cds.borders.v = vext(1:NSub)/norm(vext(1:NSub));
cds.borders.w = wext(1:NSub)/norm(wext(1:NSub));

%Step 5
Qv = Q3 * vext(1:NSub);
Qw = wext2(1:end-1);
% Qw = Q3 * wext(1:NSub);

x1 = X;
x1(1:cds.ncoo) = x1(1:cds.ncoo) - Incr/norm(Qv)*Qv;
[xt,p] = rearr(x1); p = num2cell(p);
f_x1 = [lpjac(xt,p) lpjacp(xt,p)];

x2 = X;
x2(1:cds.ncoo) = x2(1:cds.ncoo) + Incr/norm(Qv)*Qv;
[xt,p] = rearr(x2); p = num2cell(p);
f_x2 = [lpjac(xt,p) lpjacp(xt,p)];

g_x = -Qw'* (f_x2 - f_x1) * norm(Qv)/(2*Incr);

%Step 6
jac = [f_x; g_x];

%---------------------------------------------------
function hess = hessians(varargin)

%print_diag(5,'In limitpointL/hessians\n');
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

%print_diag(5,'In limitpointL/options\n');
global contopts
contopts = contset(contopts, 'Locators', [0 0 0]);

% -------------------------------------------------------
function [out, failed] = testf(id, x, v, CISdata)
print_diag(5,'In limitpointL/testf\n');

global cds           

[x0, p] = rearr(x); p = num2cell(p);

NSub   = CISdata.NSub;
Asub   = CISdata.T(1:NSub,1:NSub);
evl3_r = CISdata.evl_r;

Bord   = [Asub cds.borders.w; cds.borders.v' 0];
bunit  = [zeros(NSub,1);1];
vext   = Bord\bunit;
wext   = Bord'\bunit;
cds.borders.v = vext(1:NSub)/norm(vext(1:NSub));
cds.borders.w = wext(1:NSub)/norm(wext(1:NSub));

[~, I]   = sort(abs(evl3_r));
E = evl3_r(I);
E(1) = [];

failed = 0;


for i=id
    lastwarn('');
    switch i
        case 1 % BT
            out(1) = wext(1:NSub)'*vext(1:NSub);
        case 2 % ZH1
            NSub2    = NSub-1;
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
            out(2) = mu_min*(-1)^M_u;
        case 3 % CP
            out(3) = nf_LP_L(CISdata, x0, p, cds.borders);
        case 4 % ZH2
            out(4) = sum(real(E) > 0);   % DV: New test function
        otherwise % BP
            bp = size(cds.BranchParams,2);
            out(3+(1:bp))=wext'*lpjacbr(x0,p);
    end
    if ~isempty(lastwarn)
        msg = sprintf('Could not evaluate tf %d\n', i);
        failed = 1;
    end
    
end
%------------------------------------------------------
function [out, failed] = userf(id, UserInfo, x)
%print_diag(5,'In limitpointL/userf\n');

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
%---------------------------------------------------------
function [failed,s] = process(id, point, s)

print_diag(5,'In limitpointL/process\n');

global cds

if isfield(point, 's'); s.data = point.s; else; s.data = []; end

s.data.CISdata = point.CISdata;

[x0,p] = rearr(point.x); p = num2cell(p);

switch id
    case 1 % BT
        s.data.c = nf_BT_L(point.CISdata, x0, p); % MP
        print_diag(1,'BT quadratic coefficients (a2, b2) = (%e, %e)\n', s.data.c(1), s.data.c(2));
        s.msg  = 'Bogdanov-Takens point';
    case 2 % ZH
        s.data.c = nf_ZH_L(point.CISdata, x0, p); % MP
        print_diag(1,'ZH coefficients (B0, C0, E0) = (%e, %e, %e)\n', s.data.c(1), s.data.c(2), s.data.c(3));
        s.msg  = 'Zero-Hopf point';
    case 3 % CP
        s.data.c = nf_CP_L(point.CISdata, x0, p);
        print_diag(1,'CP cubic coefficient c = %e\n', s.data.c);
        s.msg = sprintf('Cusp point');
    otherwise
        s.msg = sprintf('Branch point (parameter %d)',cds.BranchParams(id-3));
end
s.data.eigenvals = [point.CISdata.evl_r; point.CISdata.evl_l];

failed = 0;

%--------------------------------------------------------
function  [S,L] = singmat
%print_diag(5,'In limitpointL/singmat\n');

global cds
% 0: testfunction must vanish
% 1: testfunction must not vanish
% everything else: ignore this testfunction

% S = [  0 0 8 8
%        1 0 8 0 
%        8 8 0 8 ];

S = [  0 8 8 8     % DV: Changed similar to CL_Matcont5p4
       1 0 8 2     % DV: Changed for new testfunction ZH2
       1 8 0 8 ]; 
bp = size(cds.BranchParams,2);
for i = 1:bp
    S(3+i,end+1)   = 0;
    S(3+i,1:end-1) = 8;
    S(1:end-1,3+i) = 8;
end
L = [ 'BT'; 'ZH'; 'CP'];        % rows of S
for i = 1:bp
    L(3+i,:)= strcat('BP',num2str(cds.BranchParams(i)));
end

%------------------------------------------------------
function [ps] = locate(id, p1, p2)
% DV: There are no locators available at the moment
error('No locator defined for singularity %d', id);

%------------------------------------------------------
function [varargin] = init(varargin)
% Placeholder
varargout{1} = 0;

%--------------------------------------------------------
function varargout = done
% Placeholder
varargout{1} = 0;
%---------------------------------------------------------
function [res,x,v,CISdata] = adapt(x,v,CISdata,tfUpdate)
print_diag(5,'In limitpointL/adapt\n');

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

% -------------------------------------------------------

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