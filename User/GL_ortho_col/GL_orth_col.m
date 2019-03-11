
% continuation of cycles in brusselator
odefile = @cubic_quintic_Ginzburg_Landau ;
N = 2; L = 1.1; A = 1; B = 2.2; Dx = 0.008; Dy = 0.004;
parameters = {L; A; B; Dx; Dy};
handles = feval(odefile);

integration_opt = odeset(...
  'AbsTol',      1e-3,    ...
  'RelTol',      1e-6,    ...
  'BDF',         'off',   ...
  'MaxOrder',     5,      ...
  'NormControl',  'off',  ...
  'Refine',       1,      ...
  'Jacobian',     @(t,y) feval(handles{3},t,y,parameters{:}) ...
);


x0 = [ones(N,1); ones(N,1)];
dydt = handles{2};
f =@(t, y) dydt(t, y, parameters{:});
[t1, x1] = ode15s(f, [0 400], x0, integration_opt);

draw_plots = false || false;
if draw_plots
  figure(1)
  plot(t1,x1)
  title(sprintf( ...
    'Brusselator N:%d L:%.2f A:%.2f B:%.2f Dx:%.3f Dy:%.3f', ...
    N,L,A,B,Dx,Dy));
  xlabel('t');
  ylabel('x_1, ..., x_N, y_1, ..., y_N');
end

approximate_period = 8;
[t2,x2] = ode15s(f, 0:0.001:approximate_period, x1(end,:)',integration_opt); 

if draw_plots
  figure(2)
  plot(t2,x2)
  title(sprintf( ...
    'Brusselator N:%d L:%.2f A:%.2f B:%.2f Dx:%.3f Dy:%.3f', ...
    N,L,A,B,Dx,Dy));
  xlabel('t')
  ylabel('x_1, ..., x_N, y_1, ..., y_N')
end

if draw_plots
  figure(3)
  plot(x2(:,1),x2(:,N+1))
  title(sprintf( ...
    'Brusselator N:%d L:%.2f A:%.2f B:%.2f Dx:%.3f Dy:%.3f', ...
    N,L,A,B,Dx,Dy));
  xlabel('x_1')
  ylabel('y_1')

  drawnow
end

%% Continue limit cycle from orbit
t = t2;
y = x2;
p = cell2mat(parameters);
ntst = 20;
ncol = 4;
tolerance = 1e-2;

opt = contset();
opt = contset(opt,'contL_LogFile',             1); 
opt = contset(opt,'contL_DiagnosticsLevel',    3);  
opt = contset(opt,'Backward',                  1);
opt = contset(opt,'InitStepsize',              0.01);
opt = contset(opt,'MinStepsize',            1e-5);
opt = contset(opt,'MaxStepsize',               0.1);
opt = contset(opt,'MaxCorrIters',             12);
opt = contset(opt,'MaxNewtonIters',            5);
opt = contset(opt,'FunTolerance',           1e-6);
opt = contset(opt,'VarTolerance',           1e-5);
opt = contset(opt,'contL_SmoothingAngle',   pi/30); 
opt = contset(opt,'Singularities',             1);
opt = contset(opt,'Userfunctions',             0);
opt = contset(opt,'MaxNumPoints',              50);

opt = contset(opt,'MaxTestIters',             10);
%opt = contset(opt,'contL_Testf_VarTolerance', 1e-4); 
%opt = contset(opt,'contL_Testf_FunTolerance', 1e-5); 
opt = contset(opt,'CIS_SparseSolvers',         1);
opt = contset(opt,'CIS_NStableRef',            4);
opt = contset(opt,'CIS_MaxUnstable',           5); %new
opt = contset(opt,'CIS_Ric_Cayley_Shift',      1); 
opt = contset(opt,'contL_EQ_BranchingMethod',  2); 
opt = contset(opt, 'CIS_UsingCIS',   false);
    % disable smoothing by angle:
  opt = contset(opt, 'contL_SmoothingAngle', pi/2);
opt = contset(opt,'Locators',            [0 0 0]);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename',     'testbruss_Orb_LC');
opt = contset(opt, 'Multipliers',   1);
opt = contset(opt, 'contL_DiagnosticsLevel', 5);
ap = 5;
opt.TestTolerance = 1e-5;





[x0,v0] = initOrbLC_L(odefile, t, y, p, ap, ntst, ncol,tolerance);

[sout,datafile] = contL(@limitcycleL,x0,v0,opt);
s = sout;
[xlc, vlc, ~,  mult] = loadPoint(datafile); % DV: load computed cycles
load('Data\testbruss_Orb_LC.mat')            % DV: load singular points
plotcycle(xlc,vlc,s,[size(xlc,1) 1 3]);
figure
hold on;
for i=1:2*N
  plot(xlc(end,:),mult(i, :))
end

