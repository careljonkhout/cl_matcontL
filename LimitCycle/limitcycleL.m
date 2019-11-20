function out = limitcycleL
%
% Limit cycle curve definition file for a problem in an odefile, using
% discretization by orthogonal collocation. (a.k.a. spline collocation)
% 

  out{1}  = @curve_func;
  out{2}  = @defaultprocessor;
  out{3}  = @options;
  out{4}  = @jacobian;
  out{5}  = @hessians;
  out{6}  = @testfunctions;
  out{7}  = @userf;
  out{8}  = @process_singularity;
  out{9}  = @cycle_singularity_matrix;
  out{10} = @locate;
  out{11} = @init;
  out{12} = @done;
  out{13} = @adapt;

  out{14} = @CIS_first_point;
  out{15} = @CIS_step;
  return
end
%-------------------------------------------------------------------------------
% Computes the curve function of the curve of limit cycles. If curve_func
% evaluates to a zero vector, then arg corresponds exactly to a limitcycle.
%
% arg is the continuation state vector. The continuation state
% vector contains a sequence of points on the limit cycle, the period and the
% value of the active parameter. The sequence of points on the cycle contain one
% point for every point in lds.finemsh.
function val = curve_func(arg, ~)  % unused argument is CISdata
  [x,p,T] = rearr(arg);
  val = BVP('BVP_LC_f','BVP_LC_bc','BVP_LC_ic',x,p,T);
end
%-------------------------------------------------------------------------------
% Computes the Jacobian matrix of the curve function w.r.t. all points on the
% limit cycle the parameter and the period. This Jacobian matrix is a N by N+1
% matrix.
function varargout = jacobian(varargin)
  [x, parametervalues, period] = rearr(varargin{1});
  varargout{1} = BVP_jac('BVP_LC_jac', x, parametervalues, period, 2, 2);
end
%-------------------------------------------------------------------------------
function varargout = hessians(varargin)
  varargout = {};
end
%-------------------------------------------------------------------------------
function point = defaultprocessor(varargin) 
  global lds
  point = varargin{1};
  [~, parametervalues, period] = rearr(point.x);
  update_upoldp(point.x, point.v);
  
  % calculate multipliers if requested
  if lds.CalcMultipliers 
    update_multipliers_if_needed(point.x);
  end
  if lds.CalcPRC || lds.CalcdPRC
    [lds.PRCdata, lds.dPRCdata] = calcPRC(point.x,lds.PRCInput,[0 0]);
  end
  if ~ lds.CalcPRC
    lds.PRCdata = [];
  end
  if ~ lds.CalcdPRC
    lds.dPRCdata = [];
  end
  point.multipliers     = lds.multipliers;
  point.timemesh        = lds.msh;
  point.ntst            = lds.ntst;
  point.ncol            = lds.ncol;
  point.parametervalues = parametervalues;
  point.T               = period;
  point.phi             = lds.PD_phi(lds.coords);
  point.PRCdata         = lds.PRCdata;
  point.dPRCdata        = lds.dPRCdata;
    
  if lds.CalcMultipliers == 0
    lds.multipliers = [];
  end

  savePoint(point, varargin{2:end});
end
%-------------------------------------------------------------------------------
function options
end
%-------------------------------------------------------------------------------
% Test functions are used for detecting AND locating singularities by bisection.
% When detecting ids_testf_requested will be cds.ActTest, and when locating
% ids_testf_requested will contain only those ids of the testfunctions relevant
% to the bifurcation that is being located.
function [out, failed] = testfunctions(ids_testf_requested, x0, v, ~) 
  % unused argument is CISdata
  
  global lds;
  
  failed = false;
  
  const = Constants;
  
  if any(ismember([const.BPC_id const.PD_id const.NS_id], ids_testf_requested))
    update_multipliers_if_needed(x0)
  end
  
  out = cycle_testfunctions(ids_testf_requested, lds.multipliers, v);
end
%-------------------------------------------------------------------------------
function update_multipliers_if_needed(x)
  global cds lds
  if ~ isfield(lds,'multipliersX') || isempty(lds.multipliersX) || ...
      ~ all(lds.multipliersX == x)
    lds.multipliersX = x;
    try
      lds.multipliers = multipliers(...
          cjac(cds.curve_func,cds.curve_jacobian,x,[]));     
    catch error
      % One could try to use multipliers_variational.
      % lds.multipliers = multipliers_variational(x);
      print_diag(0, 'LimitcycleL: Failed to compute multipliers\n');
      print_diag(1,  getReport(error));
    end
  end
end
%-------------------------------------------------------------------------------
function [out, failed] = userf( userinf, id, x, ~) % unused argument is v
  global  lds
  dim =size(id,2);
  failed = [];
  out(dim) = 0;
  for i=1:dim
    lastwarn('');
    [x0,p] = rearr(x); p = num2cell(p);
    if (userinf(i).state==1)
        out(i)=feval(lds.user{id(i)},0,x0,p{:});
    else
        out(i)=0;
    end
    if ~isempty(lastwarn)
      failed = [failed i]; %#ok<AGROW>
      % the array failed will be small, will not cause performance issues
    end
  end
end
%-------------------------------------------------------------------------------
function [failed, s] = process_singularity(id, point, s)
  x = point.x;
  global cds lds contopts
  switch id
  case Constants.BPC_id
    format_string = 'Branch Point cycle (period = %e, parameter = %e)\n'; 
    print_diag(0, format_string, x(end-1), x(end));
    s.msg  = sprintf('Branch Point cycle'); 
  case Constants.PD_id
    if contopts.enable_nf_pd
      [~,p,T] = rearr(x); % unused argument is x0
      J = BVP_jac('BVP_PD_jac',x,p,T,1,1);
      [LJ,UJ] = lu(J);
      b = []; b(lds.ncoords+1)=1; b=b';
      ss = UJ\(LJ\b);
      lds.PD_phi = reshape(ss(lds.coords),lds.nphase,lds.tps);
      s.data.phi = lds.PD_phi(lds.coords);
      s.data.pdcoefficient = nf_PD(x);
      format_string = 'Period Doubling (period = %e, parameter = %e)\n';
      print_diag(0, format_string, x(end-1), x(end));
      s.msg  = 'Period Doubling';
      print_diag(0, 'Normal form coefficient = %d\n', s.data.pdcoefficient);
    else
      s.data.phi = NaN;
      s.data.pdcoefficient = NaN;
      format_string = 'Period Doubling (period = %e, parameter = %e)\n';
      print_diag(0, format_string, x(end-1), x(end));
      s.msg  = 'Period Doubling';
    end
  case Constants.LPC_id
    s.msg = 'Limit point cycle';
    if contopts.enable_nf_lpc
      s.data.lpccoefficient = nf_LPC(x);
      format_string = 'Limit point cycle (period = %e, parameter = %e)\n';
      print_diag(0, format_string, x(end-1), x(end));
      print_diag(0,' Normal form coefficient = %d\n', s.data.lpccoefficient);
    else
      s.data.lpccoefficient = NaN;
      format_string = 'Limit point cycle (period = %e, parameter = %e)\n';
      print_diag(0, format_string, x(end-1), x(end));
    end
  case Constants.NS_id

    if contopts.enable_nf_ns
      % todo: rewrite to de-duplicate for reporting neutral saddle
      s.data.nscoefficient = nf_NS(x);
      if strcmp(s.data.nscoefficient, 'Neutral saddle')
        s.msg = 'Neutral saddle cycle';
        format_string = 'Neutral Saddle Cycle (period = %e, parameter = %e)\n';
        % A neutral saddle is not really a bifurcation, therefore we use
        % priority 1 instead of 0, so that it is only logged if
        % contopts.contL_DiagnosticsLevel is set higher than the default value
        % which is zero.
        print_diag(1, format_string, x(end-1), x(end));
      else
        s.msg = 'Neimark Sacker';
        format_string = 'Neimark-Sacker (period = %e, parameter = %e)\n';
        print_diag(0, format_string, x(end-1) ,x(end));
        print_diag(0, 'Normal form coefficient = %d\n', s.data.nscoefficient);
      end
    else
      d = lds.multipliers;
      smallest_sum = Inf;
      for jk=1:lds.nphase-1
        [val,idx] = min(abs(d(jk+1:lds.nphase)*d(jk)-1));
        if val < smallest_sum
          idx2 = jk+idx;
          smallest_sum = val;
        end
      end
      singularity_is_neutral_saddle = ...
        abs(imag(d(idx2))) < cds.deviation_of_trivial_multiplier;
      if singularity_is_neutral_saddle
        s.msg = 'Neutral saddle cycle';
        format_string = 'Neutral Saddle Cycle (period = %e, parameter = %e)\n';
        % A neutral saddle is not really a bifurcation, therefore we use
        % priority 1 instead of 0, so that it is only logged if
        % contopts.contL_DiagnosticsLevel is set higher than the default value
        % which is zero.
        print_diag(1, format_string, x(end-1), x(end));
      else
        s.msg = 'Neimark Sacker';
        format_string = 'Neimark-Sacker (period = %e, parameter = %e)\n';
        print_diag(0, format_string, x(end-1) ,x(end));
        %print_diag(0, 'Normal form coefficient = %d\n', s.data.nscoefficient);
      end
    end
  end
  failed = 0;
end
%-------------------------------------------------------------------------------
function p_out = locate(id, p1,p2) %#ok<INUSD,STOUT>
  switch id   
    otherwise
      error('No locator defined for singularity %d', id);
  end
end
%-------------------------------------------------------------------------------
function varargout = init(x, v, varargin)
  WorkspaceInit(x, v);
  varargout{1} = 0;
end
%-------------------------------------------------------------------------------
function varargout = done
  varargout = {};
end
%-------------------------------------------------------------------------------
function [res,x,v,CISdata] = adapt(x,v,CISdata,~) 
  % unused argument (~) is tfUpdate
  global lds cds

  % calculate phi and psi for next point

  cds.adapted = 1;

  if lds.BP_switch == 0
    lds.BP_phi = lds.BP_new_phi;
    lds.BP_phi =  lds.BP_phi/norm(lds.BP_phi(lds.coords));
    lds.BP_phi1 = lds.BP_new_phi1;
    lds.BP_phi1 = (lds.BP_phi1*lds.BP_phi')*lds.BP_phi-(lds.BP_phi*lds.BP_phi')*lds.BP_phi1;
    lds.BP_phi1 = lds.BP_phi1/norm(lds.BP_phi1(lds.coords));  
  else
    lds.BP_psi = lds.BP_new_psi;
    lds.BP_psi =  lds.BP_psi/norm(lds.BP_psi);
    lds.BP_psi1 = lds.BP_new_psi1;
    lds.BP_psi1 = (lds.BP_psi1*lds.BP_psi')*lds.BP_psi-(lds.BP_psi*lds.BP_psi')*lds.BP_psi1;
    lds.BP_psi1 = lds.BP_psi1/norm(lds.BP_psi1);
  end
  lds.BP_switch = 1-lds.BP_switch;

  [x,v] = adapt_mesh(x,v);
  res = 1;
end
%-------------------------------------------------------------------------------
%{
function [x,v] = locateBPC(~, x1, v1, ~, v2)
  % unused arguments are id and x2
  global  cds

  ndim = cds.ndim;

  initpq(x1);
  b = 0;
  x = x1;
  i = 0;

  v = 0.5*(v1+v2);
  u = [x; b];
  [A,f]=locjac(x,b);
  while i < 4
    du = A\f;
    u = u - du;

    x = u(1:ndim);
    b = u(ndim+1);

    [A,f]=locjac(x,b);
    % WM: VarTol and FunTol were switched
    print_diag(3,"locateBPC norm(du):%.5f norm(f): %.5f\n",norm(du),norm(f))
    if norm(du) < cds.options.VarTolerance && norm(f) < cds.options.FunTolerance 
        return; 
    end

    i = i+1;

  end
  x = [];
end
%-------------------------------------------------------------------------------
function [A, f] = locjac(x0, b)
  % A = jac of system
  % f = system evaluated at (x,b)
  global cds lds
  [x,p,T] = rearr(x0);
  % append g
  J = BVP_BPCjac('BVP_BPC_jacCC',x,p,T,1,1);

  b1=[zeros(lds.ncoords+1,2);eye(2)];
  sn = J\b1;
  f = [feval(cds.curve_func, x0) + b*lds.BPC_psi'; sn(end,:)'];

  A = BVP_jac('BVP_LC_jac',x,p,T,2,2);
  j = size(A,1)+1;
  A(:,j+1) = lds.BPC_psi';
  A(j,:)   = 0;
  A(j+1,:) = 0;

  b1 = []; b1(lds.ncoords+3)=1;

  st = A'\b1';
  if any(st==Inf) || any(isNaN(st))
    st = lsqminnorm(A',b1');
  end
  v11 = sn(1:lds.ncoords,1)';
  v21 = sn(1:lds.ncoords,2)';
  v12 = sn(lds.ncoords+1,1);
  v22 = sn(lds.ncoords+1,2);
  v13 = sn(lds.ncoords+2,1);
  v23 = sn(lds.ncoords+2,2);
  w1 = st(1:lds.ncoords-lds.nphase)';

  % calculate g'
  ups = reshape(x,lds.nphase,lds.tps);
  p = num2cell(p);

  range0 = lds.cols_p1;
  range1 = lds.col_coords;
  range2 = lds.cols_p1_coords;

  t = lds.nphase:((lds.ncol+2)*lds.nphase-1);
  kr1 = fix(t/lds.nphase);
  kr2 = rem(t,lds.nphase)+1;

  for tstpt = lds.tsts
      xp  = ups(:,range0)*lds.wt;
      cv1 = v11(range2)';
      cv2 = v21(range2)';
      cw1 = w1(range1);

      range = lds.phases;
      for c=lds.cols
          xt = xp(:,c);
          sysj   = cjac(lds.func,lds.Jacobian,xt,p,lds.ActiveParams);
          sysp   = cjacp(lds.func,lds.JacobianP,xt,p,lds.ActiveParams);
          sysh   = chess(lds.func,lds.Jacobian,lds.Hessians,xt,p,lds.ActiveParams);
          syshp  = chessp(lds.func,lds.Jacobian,lds.HessiansP,xt,p,lds.ActiveParams);
          sysbr  = cjacbr(lds.func,lds.JacobianP,xt,p,lds.ActiveParams,lds.BranchParam);
          syshbr = chessbr(lds.func,lds.Jacobian,lds.HessiansP,xt,p,lds.ActiveParams,lds.BranchParam);
          syshbrp = chesspbr(lds.func,lds.JacobianP,lds.HessiansP,xt,p,lds.ActiveParams,lds.BranchParam);
          wtk = lds.wt(kr1,c(ones(1,lds.nphase)))';
          for d=lds.phases
              sh1(:,d) = (wtk.*sysh(:,kr2,d))*cv1;
              sh2(:,d) = (wtk.*sysh(:,kr2,d))*cv2;
          end    
          t11 = T* wtk.*sh1(:,kr2) + wtk.*sysj(:,kr2)*v12 + T*wtk.*syshbr(:,kr2,1)*v13;
          t12 = T* wtk.*sh2(:,kr2) + wtk.*sysj(:,kr2)*v22 + T*wtk.*syshbr(:,kr2,1)*v23;
          t21 = (wtk.*sysj(:,kr2))*cv1 + sysbr*v13;
          t22 = (wtk.*sysj(:,kr2))*cv2 + sysbr*v23;
          t31 = T* wtk.*syshp(:,kr2,1)* cv1 + sysp(:,1)*v12 + T*syshbrp(:,:,1)*v13;
          t32 = T* wtk.*syshp(:,kr2,1)* cv2 + sysp(:,1)*v22 + T*syshbrp(:,:,1)*v23;
          syshess1(range,:) = [t11 t21 t31 zeros(size(t11,1),1)];      
          syshess2(range,:) = [t12 t22 t32 zeros(size(t12,1),1)];
          range = range + lds.nphase;
      end

      A(j,[range2 lds.ncoords+(1:3)])   = A(j,[range2 lds.ncoords+(1:3)])   + cw1*syshess1;
      A(j+1,[range2 lds.ncoords+(1:3)]) = A(j+1,[range2 lds.ncoords+(1:3)]) + cw1*syshess2;
      range0 = range0 + lds.ncol;
      range1 = range1 + lds.ncol_coord;
      range2 = range2 + lds.ncol_coord;
  end
end
%-------------------------------------------------------------------------------
function  initpq(x0)
  % A = jac of system
  % f = system evaluated at (x,b,p)
  global lds

  [x,p,T] = rearr(x0);

  % append g
  ups = reshape(x,lds.nphase,lds.tps);
  p = num2cell(p);
  pars1 = lds.ncoords+1;
  pars2 = lds.ncoords+2;

  jac = spalloc(lds.ncoords+1,lds.ncoords+2,(lds.ncol+4)*lds.nphase);
  % function
  range0 = lds.cols_p1;
  range1 = lds.col_coords;
  range2 = lds.cols_p1_coords;
  for j = lds.tsts
    xp = ups(:,range0)*lds.wt;
    jac(range1,[range2 pars1 pars2]) = bordBVP_BPC_f(lds.odefile,xp,p,T,j);
    range0 = range0 + lds.ncol;
    range1 = range1 + lds.ncol_coord;
    range2 = range2 + lds.ncol_coord;
  end
  % boundary conditions
  range  = (lds.tps-1)*lds.nphase+ (lds.phases);
  range1 = lds.ncoords-lds.nphase+lds.phases;
  jac(range,[lds.phases range1]) = bordBVP_LPC_bc1;
  % integral constraint
  ic = zeros(1,lds.ncoords);
  range1 = lds.cols_p1;
  range2 = lds.cols_p1_coords;
  for j=lds.tsts
    pt = lds.dt(j)*(ups(:,range1).*lds.pwi);
    ic(range2) = ic(range2)+pt(lds.cols_p1_coords);
    range1 = range1 + lds.ncol;
    range2 = range2 + lds.ncol_coord;
  end
  jac(range(end)+1,1:lds.ncoords)= ic;
  %compute borders
  [Q,R,E] = qr(jac);
  R(end,end) = 0;
  R(end,end-1) = 0;
  p = E*[R(1:end-1,1:end-2)\-R(1:end-1,end-1:end);eye(2)];
  p = full(p);
  p = p'/norm(p);
  lds.BPC_phi1=p(1,:);
  lds.BPC_phi2=p(2,:);
  p = Q(:,end);
  lds.BPC_psi = p';
end
%}
%-------------------------------------------------------------------------------
function [x,p,T] = rearr(x0)
  %
  % [x,p] = rearr(x0)
  %
  % Rearranges x0 into coordinates (x) and parameters (p)
  global lds

  p = lds.P0;
  if length(lds.ActiveParams) == 1
      p(lds.ActiveParams) = x0(lds.PeriodIdx+1);         
      x = x0(lds.coords);                                  
      T = x0(lds.PeriodIdx);
  else
      p(lds.ActiveParams) = x0(lds.PeriodIdx:lds.PeriodIdx+1);
      x = x0(lds.coords);
      T = lds.T;
  end
end
%-------------------------------------------------------------------------------
function f = BVP(BVP_f,BVP_bc,BVP_ic,x,p,T)
  global lds 

  % extract ups
  ups = reshape(x,lds.nphase,lds.tps);
  p = num2cell(p);

  % function
  range1 = lds.cols_p1;
  range2 = lds.phases;
  f = zeros((lds.ncol* lds.ntst + 1) * lds.nphase + 1, 1);
  for j=lds.tsts
    % value of polynomial on each collocation point
    xp = ups(:,range1)*lds.wt;

    % derivative of polynomial on each collocation point
    t  = ups(:,range1)*lds.wpvec/lds.dt(j);

    % evaluate function value on each collocation point
    for c=lds.cols
      f(range2) = feval(BVP_f,lds.func,t(:,c),xp(:,c),p,T);
      range2 = range2+lds.nphase;
    end

    range1 = range1+lds.ncol;
  end
  % boundary conditions
  f(range2) = feval(BVP_bc,ups(:,1),ups(:,lds.tps));

  % integral constraint
  f(end) = feval(BVP_ic,ups);
end
%-------------------------------------------------------------------------------
function jacx = BVP_jac(BVP_func,x,p,T,pars,nc)
  global lds
  
  print_diag(4,'running %s ... ',BVP_func);
  
  p2 = num2cell(p);
  tic_time = tic;
  jacx = feval(BVP_func,lds.func,x,p,T,pars,nc,lds,p2,lds.Jacobian, ...
                lds.ActiveParams,lds.JacobianP);
	time_elapsed = toc(tic_time);
  print_diag(4,'time elapsed %.4f\n', time_elapsed);
end
%-------------------------------------------------------------------------------
function WorkspaceInit(x,v)
  global lds contopts
  ntst               = lds.ntst;
  ncol               = lds.ncol;
  nphase             = lds.nphase;
  ncoords            = lds.ncoords; % == (ntst * ncol + 1) * nphase
  lds.cols_p1        = 1:(ncol+1);
  lds.cols_p1_coords = 1:((ncol + 1) * nphase);
  lds.ncol_coord     = ncol * nphase;
  lds.col_coords     = 1:(ncol * nphase);
  lds.coords         = 1:ncoords;
  lds.pars           = ncoords+(1:2);
  lds.tsts           = 1:ntst;
  lds.cols           = 1:ncol;
  lds.phases         = 1:nphase;
  lds.ntstcol        = ntst * ncol;
  % fix rounds towards zero
  lds.idxmat         = reshape(fix((1:((ncol+1)*ntst))/(1+1/ncol))+1,ncol+1,ntst);
  lds.dt             = lds.msh(lds.tsts+1)-lds.msh(lds.tsts);

  lds.wp             = kron(lds.wpvec',eye(nphase));
  lds.pwwt           = kron(lds.wt',eye(nphase));
  lds.pwi            = lds.wi(ones(1,nphase),:);

  lds.wi             = nc_weight(ncol)';


  lds.PD_psi = reshape(exp((1:ncoords)/ncoords),nphase,lds.tps);
  lds.PD_psi = lds.PD_psi/norm(lds.PD_psi(lds.coords));
  lds.PD_phi = reshape(ones(lds.ncoords,1),lds.nphase,lds.tps);
  lds.PD_phi = lds.PD_phi/sqrt(BVP_PD_jac_ic*lds.PD_phi(lds.coords)');

  lds.PD_new_phi = lds.PD_phi;
  lds.PD_new_psi = lds.PD_psi;
  lds.PD_switch = 0;

  lds.BP_psi = reshape(exp((1:ncoords)/ncoords),nphase,lds.tps);
  lds.BP_psi = lds.BP_psi/norm(lds.BP_psi(lds.coords));
  lds.BP_phi = reshape(ones(ncoords,1),nphase,lds.tps);
  lds.BP_phi = lds.BP_phi/sqrt(BVP_BPC_jac_ic*lds.BP_phi(lds.coords)');

  lds.BP_psi1 = reshape(ones(ncoords,1),nphase,lds.tps);
  lds.BP_psi1 = (lds.BP_psi1*lds.BP_psi')*lds.BP_psi-(lds.BP_psi*lds.BP_psi')*lds.BP_psi1;
  lds.BP_psi1 = lds.BP_psi1/norm(lds.BP_psi1(lds.coords));

  lds.BP_phi1 = reshape(exp((1:ncoords)/ncoords),nphase,lds.tps);
  lds.BP_phi1 = (lds.BP_phi1*lds.BP_phi')*lds.BP_phi-(lds.BP_phi*lds.BP_phi')*lds.BP_phi1;
  lds.BP_phi1 = lds.BP_phi1/norm(lds.BP_phi1(lds.coords));

  lds.BP_new_phi  = lds.BP_phi;
  lds.BP_new_psi  = lds.BP_psi;
  lds.BP_new_psi1 = lds.BP_psi1;
  lds.BP_new_phi1 = lds.BP_phi1;
  lds.BP_switch   = 0;

  lds.BPC_switch = 0;
  lds.BPC_psi    = reshape(ones(ncoords+1,1),1,ncoords+1);
  lds.BPC_phi1   = reshape(ones(ncoords+2,1),1,ncoords+2);
  lds.BPC_phi2   = reshape(ones(ncoords+2,1),1,ncoords+2);

  lds.LPC_phi=[];
  lds.LPC_psi=[];
  lds.LPC_new_phi = lds.LPC_phi;
  lds.LPC_new_psi = lds.LPC_psi;
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


  lds.CalcMultipliers = contopts.Multipliers || contopts.Singularities;
  lds.CalcPRC = contopts.PRC;
  lds.CalcdPRC = contopts.dPRC;
  lds.PRCInput = contopts.Input;
  lds.multipliersX = [];
  lds.multipliers = nan;

  if isempty(v)
    v = zeros(size(x));
  end
  update_upoldp(x,v);
end
%-------------------------------------------------------------------------------
function jac = BVP_BPCjac(BVP_func,x,p,T,pars,nc) %#ok<DEFNU> 
  global lds
  p2 = num2cell(p);
  print_diag(4,'running %s\n', BVP_func)
  jac = feval(BVP_func,lds.func,x,p,T,pars,nc,lds,p2,lds.Jacobian,...
              lds.ActiveParams,lds.JacobianP,lds.BranchParam); 
end
%-------------------------------------------------------------------------------
function CISdata = CIS_first_point(~) % unused argument is x
  CISdata = 1;
end
%-------------------------------------------------------------------------------
function CISdata = CIS_step(~, ~) % unused arguments are x and CISdata1
  CISdata = 1;
end
%-------------------------------------------------------------------------------
function update_upoldp(x, v)
  global lds           
  [x,p,T] = rearr(x); 
  v       = rearr(v);
  lds.ups = reshape(x,lds.nphase,lds.tps);
  lds.vps = reshape(v,lds.nphase,lds.tps);
  % update upoldp
  p1 = num2cell(p);
  for i=1:lds.tps
    lds.upoldp(:,i) = T * feval(lds.func, 0, lds.ups(:,i), p1{:});
  end
end
%-------------------------------------------------------------------------------