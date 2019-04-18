function [x0,v0] = initOrbLC_old(odefile, t, y, p, ap, ntst, ncol,~)
% unused argument is tolerance

%
%
% Initializes a limit cycle continuation from a computed orbit stored in t
% (mumber of points x 1) and y (number of points x 3) 
% 
%
global cds lds 

% check input
n_par = size(ap,2);
if n_par~=1
    error('One active parameter and the period are needed for this limit cycle continuation');
end
lds = [];
% initialize lds
if isempty(cds) || ~isfield(cds,'options')
    cds.options = contset();
end
init_EP_EP_L(odefile,y(1,:)',p,ap); 
init_lds(odefile,y,p,ap,ntst,ncol);


func_handles = feval(odefile);
symord = 0; 
symordp = 0;

if     ~isempty(func_handles{9}),   symord = 5; 
elseif ~isempty(func_handles{8}),   symord = 4; 
elseif ~isempty(func_handles{7}),   symord = 3; 
elseif ~isempty(func_handles{5}),   symord = 2; 
elseif ~isempty(func_handles{3}),   symord = 1; 
end
if     ~isempty(func_handles{6}),   symordp = 2; 
elseif ~isempty(func_handles{4}),   symordp = 1; 
end
cds.options = contset(cds.options, 'SymDerivative', symord);
cds.options = contset(cds.options, 'SymDerivativeP', symordp);
cds.symjac = 1;
cds.symhess = 0;
cds.probfile = odefile;
cds.nap = length(ap);
cds.curve_CIS_step = [];
cds.integration = @ode15s; % needed for multipliers_variational

x=y';
norms       = vecnorm(x-x(:,1));
differences = diff(norms); % i.e. differences = [ norms(2) - norms(1) ...

minimum_norm       = Inf;%tolerance;
return_point_index = -1;

for i=2:size(x,2)
  % If the difference is negative we are returning to the starting point
  % and hence may find the point where the cycle closes.
  if differences(i-1) < 0 && norms(i) < minimum_norm
    minimum_norm = norms(i);
    return_point_index = i;
  end
end

if return_point_index == -1
  error('No cycle can be found');
end

x=x(:,1:return_point_index);
t=t(1:return_point_index)';
t_scaled_to_unit_interval=(t-t(1))/(t(end)-t(1));

[lds.msh,~,~] = ...
  newmeshcycle(x,t_scaled_to_unit_interval,size(x,2)-1,1,ntst,ncol);



lds.finemsh = zeros(1,ntst*ncol+1);
fine_mesh_index = 2;

for i=1:ntst
  mesh_interval_width = lds.msh(i+1) - lds.msh(i);
  fine_mesh_dt        = mesh_interval_width / ncol;
  for j=1:ncol
    lds.finemsh(fine_mesh_index) = lds.msh(i) + j * fine_mesh_dt;
    fine_mesh_index              = fine_mesh_index + 1;
  end
end
x = interp1(t_scaled_to_unit_interval,x',lds.finemsh);
x = reshape(x',size(x,2)*size(x,1),1);
x(end+1) = t(end)-t(1);

lds.ntst = ntst;
lds.ncol = ncol;
lds.T = x(end);
x(end+1) = p(ap);
x0=x;
v0=[];
%-----------------------------------------------------------------
function init_lds(odefile,y,p,ap,ntst,ncol)
global lds cds
lds.odefile = odefile;
func_handles = feval(lds.odefile);
lds.func = func_handles{2};
lds.Jacobian  = func_handles{3};
lds.JacobianP = func_handles{4};
lds.Hessians  = func_handles{5};
lds.HessiansP = func_handles{6};
lds.Der3=func_handles{7};
lds.Der4=func_handles{8};
lds.Der5=func_handles{9};
symord = 0; 
if ~isempty(func_handles{9}),   symord = 5; 
elseif ~isempty(func_handles{8}),   symord = 4; 
elseif ~isempty(func_handles{7}),   symord = 3; 
elseif ~isempty(func_handles{5}),   symord = 2; 
elseif ~isempty(func_handles{3}),   symord = 1; 
end
cds.options = contset(cds.options, 'SymDerivative', symord);
siz = size(func_handles,2);
if siz > 9
    j=1;
    for k=10:siz
        lds.user{j}= func_handles{k};
        j=j+1;
    end
else; lds.user=[];
end
lds.nphase = size(y,2);
lds.ActiveParams = ap;
lds.P0 = p;
set_ntst_ncol(ntst,ncol,(0:ntst)/ntst);
lds.T = [];
lds.cols_p1 = 1:(lds.ncol+1);
lds.cols_p1_coords = 1:(lds.ncol+1)*lds.nphase;
lds.ncol_coord = lds.ncol*lds.nphase;
lds.col_coords = 1:lds.ncol*lds.nphase;
lds.pars = lds.ncoords+(1:3);
lds.phases = 1:lds.nphase;
lds.ntstcol = lds.ntst*lds.ncol;
lds.wp = kron(lds.wpvec',eye(lds.nphase));
lds.pwwt = kron(lds.wt',eye(lds.nphase));
lds.pwi = lds.wi(ones(1,lds.nphase),:);

lds.PD_psi = [];
lds.PD_phi = [];
lds.PD_new_phi = [];
lds.PD_new_psi = [];
lds.PD_switch = 0;

lds.BP_psi = [];
lds.BP_phi = [];
lds.BP_psi1 = [];
lds.BP_phi1 = [];
lds.BP_new_phi = [];
lds.BP_new_psi = [];
lds.BP_new_psi1 = [];
lds.BP_new_phi1 = [];
lds.BP_switch = 0;

lds.BPC_switch = 0;
lds.BPC_psi = [];
lds.BPC_phi1 = [];
lds.BPC_phi2 = [];

lds.LPC_phi = [];
lds.LPC_psi = [];
lds.LPC_new_phi = [];
lds.LPC_new_psi = [];
lds.LPC_switch = 0;

lds.NS_psi0 = [];
lds.NS_psi1 = [];
lds.NS_phi0 = [];
lds.NS_phi1 = [];
lds.NS1_new_phi = [];
lds.NS2_new_phi = [];
lds.NS1_new_psi = [];
lds.NS2_new_psi = [];
lds.NS_new_phi = [];
lds.NS_new_psi = [];
lds.NS_switch = 0;
lds.NS1_switch = 0;
lds.NS2_switch = 0;


lds.bialt_M1 = [];
lds.bialt_M2 = [];
lds.bialt_M3 = [];
lds.bialt_M4 = [];
lds.multipliers = nan;
lds.monodromy = [];
lds.multi_r1 = [];
lds.multi_r2 = [];
lds.BranchParam = lds.ActiveParams;
lds.ups = [];
lds.vps = [];
lds.BranchParams=[];
lds.tsts = 1:lds.ntst;
lds.cols = 1:lds.ncol;

