function testbruss_Orb_LC()
% Test script for 1d brusselator
%% Options
opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',          1);
opt = contset(opt,'contL_DiagnosticsLevel', 3);
opt = contset(opt,'Backward',               0); %new
opt = contset(opt,'InitStepsize',         0.1);
opt = contset(opt,'MaxStepsize',          0.1);
opt = contset(opt,'MinStepsize',         1e-5);
opt = contset(opt,'Singularities',          0);
opt = contset(opt,'Userfunctions',          0);
opt = contset(opt,'MaxNumPoints',          10);
opt = contset(opt,'Locators',         [0 0 0]);
opt = contset(opt,'CIS_SparseSolvers',      1);
opt = contset(opt,'CIS_MaxUnstable',        5); %new
opt = contset(opt,'CIS_NStableRef',         6);
opt = contset(opt,'CIS_NExtra',             6);
opt = contset(opt,'CIS_Ric_Cayley_Shift',  10);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename', 'testbruss_HP0');

%% Start point
N = 4; L = 1.1; A = 1; B = 3; Dx = 0.008; Dy = 0.004;
parameters = {N; L; A; B; Dx; Dy};
handles = bruss_1d;

integration_opt = odeset(     'AbsTol', 1e-10);
integration_opt = odeset(integration_opt, 'RelTol', 1e-13);
integration_opt = odeset(integration_opt, 'Jacobian', @(t,y) feval(handles{3},t,y,parameters{:}));

x0 = [ones(N,1); ones(N,1)];
dydt = handles{2};
f =@(t, y) dydt(t, y, parameters{:});
[t1, x1] = ode15s(f, [0 30], x0, integration_opt);  %#ok<*ASGLU>
draw_plot = false || false;
if draw_plot
  plot(t1,x1)
  return
end
approximate_period = 7;
[t2, x2] = ode15s(f, [0 approximate_period],x1(end,:)' ,integration_opt); 
opt.TSearchOrder = 0;
%% Continue limit cycle from Hopf
tolerance = 1e-3;
ntst = 3;
ncol = 4;
odefile = @bruss_1d;
y = x2';
ap = 4;
[x0,v0] = initOrbLC(odefile, t2, y, cell2mat(parameters), ap, ntst, ncol,tolerance);
% to circumvent bialtaa error comment line 603 of limitcycleL.m
[~, datafile] = contL(@limitcycleL,x0,v0,opt);

