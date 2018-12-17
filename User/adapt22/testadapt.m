function testadapt()
%global x v s h f opt
%disp('>> init;');

%disp('>> [x0,v0]=init_EP_EP(@adaptx,[0;0;0],[-10;1],[1]);');
%disp('>> opt = contset; opt = contset(opt,''Singularities'',1);');
%opt=contset;opt=contset(opt,'Singularities',1);

opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',               1);
opt = contset(opt,'contL_DiagnosticsLevel',      3);
opt = contset(opt,'Backward',                    0);  % {0,1} = {Backward, Forward}
opt = contset(opt,'InitStepsize',              0.2);   
opt = contset(opt,'MaxStepsize',               0.5);
%opt = contset(opt,'MaxStepsize',               2.0); 
opt = contset(opt,'MinStepsize',              1e-5);   
opt = contset(opt,'Singularities',               1); 
opt = contset(opt,'MaxNumPoints',               40); 

opt = contset(opt,'CIS_SparseSolvers',           0);
opt = contset(opt,'CIS_NUnstable',              -1);
opt = contset(opt,'CIS_NExtra',                  0);
opt = contset(opt,'CIS_NSub',                    3);  % NUnstable + NStableRef
opt = contset(opt,'CIS_resizeable_flag',         0);  % resize subspace MP 2018
opt = contset(opt,'CIS_Ric_SubspaceSelect',  'eig');  % (default 'ric') MP 2018 
opt = contset(opt,'CIS_DetectOverlap',           0);

opt = contset(opt,'Locators',      [1 1 1]); % new
opt = contset(opt,'MaxTestIters',               30); % new
opt = contset(opt,'contL_Testf_FunTolerance', 1e-7); 
opt = contset(opt,'contL_Testf_VarTolerance', 1e-6); 

%  matcont defaults:
%    cds.options.FunTolerance      = contget(cds.options, 'FunTolerance', 1e-6);
%    cds.options.VarTolerance      = contget(cds.options, 'VarTolerance', 1e-6);
%    cds.options.TestTolerance     = contget(cds.options, 'TestTolerance', 1e-5);

opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt,'Filename','testadapt');



%% Continuation 
%[x0,v0]=init_EP_EP(@adapt22,[0;0;0],[-10;1],[1]);
%disp('>> [x,v,s,h,f]=cont(@equilibrium,x0,[],opt);');
%[x,v,s,h,f]=cont(@equilibrium,x0,[],opt);
%N = 3; 
par_alpha = -10; par_beta = 1;
p = [par_alpha, par_beta]; ap1 = 1;

[x0,v0]=init_EP_EP_L(@adapt22,[], p, ap1);
contL(@equilibriumL,x0,v0,opt); 

%% Plot results
x = loadPoint(fullfile('Data','testadapt.dat'));
load(fullfile('Data','testadapt.mat'),'s')
%N = s(1).data.P0(1);
N = 3;
plot(x(N+1, :), x(1, :));
hold on
for sii = s
    plot(x(N+1, sii.index), x(1, sii.index), 'r.');
    text(x(N+1, sii.index), x(1, sii.index), sii.label);
end
