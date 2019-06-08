function testbruss_LP1()
% Test script for 1d brusselator (original MATCONT data)

% This example performs a continuation of limit points starting from a
% point located in testbruss_LP0

%% Options
opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',               1);
%opt = contset(opt,'contL_DiagnosticsLevel',     3);
opt = contset(opt,'contL_DiagnosticsLevel',      5);
opt = contset(opt,'Backward',                    0);
opt = contset(opt,'InitStepsize',                3);
opt = contset(opt,'MaxStepsize',                 3);
opt = contset(opt,'MinStepsize',              1e-8);
opt = contset(opt,'MaxCorrIters',               12);
opt = contset(opt,'MaxNewtonIters',              5);
opt = contset(opt,'FunTolerance',             1e-5);
opt = contset(opt,'VarTolerance',             1e-5);
opt = contset(opt,'contL_SmoothingAngle',    pi/30);
opt = contset(opt,'Singularities',               1);
opt = contset(opt,'Userfunctions',               0);
opt = contset(opt,'MaxNumPoints',              150);
opt = contset(opt,'Locators',              [1 1 1]);
opt = contset(opt,'contL_Testf_FunTolerance', 1e-4);
opt = contset(opt,'contL_Testf_VarTolerance', 1e-4);
opt = contset(opt,'MaxTestIters',               16);
opt = contset(opt,'CIS_SparseSolvers',           1);
opt = contset(opt,'CIS_NStableRef',              9);
opt = contset(opt,'CIS_MaxUnstable',             5); % MP new
opt = contset(opt,'CIS_NExtra',                  4);
opt = contset(opt,'CIS_Ric_Cayley_Shift',       10);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename',     'testbruss_LP1');

%% Continuation
path_to_this_script = get_path;
LP0_file = [path_to_this_script, 'Data/testbruss_LP0.mat'];
load(LP0_file, 's');
ap = [2, 3];
ID = 2;
if(strcmp(s(ID).label ,'LP'))
  data = s(ID).data;
  %
  [x0,v0]      = init_LP_LP_L(@brusselator_1d, [], [], ap, data);
  contL(@limitpointL,x0,v0,opt);

  %% Plot results

  datafile         = [path_to_this_script 'Data/testbruss_LP1.dat'];
  singularity_file = [path_to_this_script 'Data/testbruss_LP1.mat'];

  x = loadPoint(datafile);
  load(singularity_file, 's');

    
  plot(x(end-1, :), x(end, :))
  hold on
  for ii = 1:length(s)
      xii = s(ii).data.x;
      plot(xii(end-1), xii(end), 'r.')
      text(xii(end-1), xii(end), s(ii).label)
  end
  xlabel 'l'
  ylabel 'a'
end    