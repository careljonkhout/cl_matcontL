function out = archcurve
  
  out{1} = @init;
  out{2} = @f;
  out{3} = @jacobian;
  out{4} = []; %@jacobianp;
  out{5} = [];
  out{6} = [];
  out{7} = [];
  out{8} = [];
  out{9} = [];
  out{10} = [];

% --------------------------------------------------------------------------
function y0 = init(p,density)

  global arch_mesh_setup

  [nodex, ix, Iactive, cdof, vdof] = arch_mesh(density);
  nodeu = 0*nodex;
  E = 1e3;
  nu = 0.3;

  arch_mesh_setup = [];
  arch_mesh_setup.nodex   = nodex;
  arch_mesh_setup.ix      = ix;
  arch_mesh_setup.Iactive = Iactive;
  arch_mesh_setup.cdof    = cdof;
  arch_mesh_setup.vdof    = vdof;
  arch_mesh_setup.E       = E;
  arch_mesh_setup.nu      = nu;

  y0      = zeros(length(Iactive), 1);

% --------------------------------------------------------------------------
function dydt = f(t,x,p,density)

  global arch_mesh_setup
  nodex   = arch_mesh_setup.nodex;
  ix      = arch_mesh_setup.ix;
  Iactive = arch_mesh_setup.Iactive;
  cdof    = arch_mesh_setup.cdof;
  E       = arch_mesh_setup.E;
  nu      = arch_mesh_setup.nu;

  nodeu = 0*nodex;
  nodeu(Iactive) = x;
  nodeu(cdof)    = -p;
  
  R = eval_Kf(nodeu, nodex, ix, E, nu);
  dydt = R(Iactive);
  
% --------------------------------------------------------------------------
function K = jacobian(t,x,p,density)

  global arch_mesh_setup
  nodex   = arch_mesh_setup.nodex;
  ix      = arch_mesh_setup.ix;
  Iactive = arch_mesh_setup.Iactive;
  cdof    = arch_mesh_setup.cdof;
  E       = arch_mesh_setup.E;
  nu      = arch_mesh_setup.nu;

  nodeu = 0*nodex;
  nodeu(Iactive) = x;
  nodeu(cdof)    = -p;
  
  [R,K] = eval_Kf(nodeu, nodex, ix, E, nu);
  K = -K(Iactive,Iactive);

% --------------------------------------------------------------------------
function jacp = jacobianp(x,p,density,cdof)

  global arch_mesh_setup
  nodex   = arch_mesh_setup.nodex;
  ix      = arch_mesh_setup.ix;
  Iactive = arch_mesh_setup.Iactive;
  cdof    = arch_mesh_setup.cdof;
  E       = arch_mesh_setup.E;
  nu      = arch_mesh_setup.nu;

  nodeu = 0*nodex;
  nodeu(Iactive) = x;
  nodeu(cdof)    = -p;
  
  [R,K] = eval_Kf(nodeu, nodex, ix, E, nu);
  K = -K(Iactive,cdof);
  
% --------------------------
function [nodex, ix, Iactive, cdof, vdof] = arch_mesh(density)

  nx = 16*density;
  ny = 1*density;
  [ix, nodex] = block2d9(pi,9.75, 0,10.75, nx,ny);
  idf     = reshape(1:prod(size(nodex)), size(nodex));

  % Index computations
  mx = 2*nx+1;
  my = 2*ny+1;
  nbc1 = (my*0    + ny)+1;  % Left midpoint
  nbc2 = (my*2*nx + ny)+1;  % Right midpoint
  nmid = (my*nx   +  0)+1;  % Bottom center
  nobs = (my*nx/2 +  0)+1;  % Bottom quarter

  % Boundary conditions
  idf(:, nbc1) = 0;
  idf(:, nbc2) = 0;

  % Center node (for control) and second node (quarter way along)
  cdof = idf(2,nmid);
  vdof = idf(2,nobs);

  idf(cdof) = 0;
  Iactive = find(idf);
  nodex = [cos(nodex(1,:)).*nodex(2,:);
           sin(nodex(1,:)).*nodex(2,:)];
