function [x,v] = init_LC_LC_L(odefile, x, v, s, par, ap, ntst, ncol)
% 
% [x0,v0] = init_LC_LC(odefile, x, v, s, ap)
%
% Initializes a limit cycle continuation from a cycle calculated
% in a previous run.
%
global lds cds
% Make sure size of x and v vectors is right for calling init_LC_LC
% x = x(1:lds.ncoords+2,:);
% v = v(1:lds.ncoords+2,:);

% check input
n_par = size(ap,2);
if n_par~=1 && n_par~= 2
    error('One active parameter and the period or 2 active parameters are needed for limt cycle continuation');
end

% initialize lds XXX
oldlds = lds;
lds = [];


odefile_handles = feval(odefile);
symord = 0; 
symordp = 0;

if     ~isempty(odefile_handles{9}),   symord = 5; 
elseif ~isempty(odefile_handles{8}),   symord = 4; 
elseif ~isempty(odefile_handles{7}),   symord = 3; 
elseif ~isempty(odefile_handles{5}),   symord = 2; 
elseif ~isempty(odefile_handles{3}),   symord = 1; 
end
if     ~isempty(odefile_handles{6}),   symordp = 2; 
elseif ~isempty(odefile_handles{4}),   symordp = 1; 
end
if ~isfield(cds,'options')
    cds.options = contset();
end
cds.options = contset(cds.options, 'SymDerivative', symord);
cds.options = contset(cds.options, 'SymDerivativeP', symordp);
cds.symjac = 1;
cds.symhess = 0;
cds.probfile = odefile;
cds.P0 = par;
cds.ncoo = length(x) - 1;
cds.nap = length(ap);
cds.ActiveParams = ap;
cds.usernorm = odefile_handles{10};
cds.userfunc = [];
for i = 11:length(odefile_handles)
    cds.userfunc{i-10} = odefile_handles{i};
end

lds.odefile = odefile;
lds.func = odefile_handles{2};
lds.Jacobian  = odefile_handles{3};
lds.JacobianP = odefile_handles{4};
lds.Hessians  = odefile_handles{5};
lds.HessiansP = odefile_handles{6};
lds.Der3 = odefile_handles{7};
lds.Der4 = odefile_handles{8};
lds.Der5 = odefile_handles{9};
siz = size(odefile_handles,2);
if siz > 10 % DV: was 9
    j=1;
    for i=10:siz
        lds.user{j}= odefile_handles{i};
        j=j+1;
    end
else
    lds.user=[];
end
lds.nphase = round((size(x,1)-2)/(s.data.ntst*s.data.ncol+1));
lds.ActiveParams = ap;
lds.P0 = par;
set_ntst_ncol(s.data.ntst,s.data.ncol,s.data.timemesh);
if isfield(s.data,'T')
    lds.T = s.data.T;
else
    lds.T = oldlds.T;
end
% get x and v
x = x(:,s.index);
if ~isempty(v)
    v = v(:,s.index);
end
if n_par==2
    x(end-n_par+1:end)=lds.P0(ap);
else
    x(end-n_par+1) = lds.T;
    x(end) = lds.P0(ap);
end
% generate a new mesh and interpolate
[x,v]=new_mesh(x,v,ntst,ncol);
%v=[];

 %lds.ups = []; % DV: this apparently also important to initialize
 %lds.vps = [];
 %lds.BranchParams=[]; 
