function testadapt1()
%global x v s h f xlc vlc slc hlc flc opt
%testadapt;
opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',               1);
opt = contset(opt,'contL_DiagnosticsLevel',      3);
opt = contset(opt,'Backward',                    0);  % {0,1} = {Backward, Forward}
opt = contset(opt,'InitStepsize',              0.01);   
opt = contset(opt,'MaxStepsize',               0.1);
opt = contset(opt,'MinStepsize',              1e-5);    
opt = contset(opt,'Singularities',               0); 
opt = contset(opt,'MaxNumPoints',              200); 

opt = contset(opt,'CIS_UsingCIS',                 0);

%opt = contset(opt,'CIS_SparseSolvers',           1);  % new
%opt = contset(opt,'CIS_NUnstable',              -1);
%opt = contset(opt,'CIS_NExtra',                  0);
%opt = contset(opt,'CIS_NSub',                    3);  % NUnstable + NStableRef
%opt = contset(opt,'CIS_resizeable_flag',         0);  % resize subspace MP 2018
%opt = contset(opt,'CIS_Ric_SubspaceSelect',  'eig');  % (default 'ric') MP 2018 
%opt = contset(opt,'CIS_DetectOverlap',           0);

opt = contset(opt,'Locators',      [1 0 0 0]); % new
opt = contset(opt,'MaxTestIters',               30); % new
opt = contset(opt,'contL_Testf_FunTolerance', 1e-7); 
opt = contset(opt,'contL_Testf_VarTolerance', 1e-6); 

opt = contset(opt,'Multipliers',                 1);
opt = contset(opt,'Adapt',                       1);

opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt,'Filename','testadapt1');


%% Continuation
load('Data\testadapt')
ap = [1];
ID = 2;
opt = contset(opt,'Singularities',               1);  
if(strcmp(s(ID).label ,'H '))
    data = s(ID).data;
    x  = data.x;    
    x1 = x(1:end-1);
    P0 = data.P0;
    p  = [x(end);P0(2)];
    
    h = 1e-6;
    ntst = 20;
    ncol = 4;
    %
    [x0,v0] = init_H_LC_L(@adapt22, x1, p, ap, h, ntst, ncol);
    [~, datafile] = contL(@limitcycleL,x0,v0,opt);
end


disp('>> plotcycle(xlc,vlc,slc,[size(xlc,1) 1 2]);');
figure
axes
[xlc, vlc, ~] = loadPoint(datafile); % DV: load computed cycles
load('Data\testadapt1.mat')            % DV: load singular points
plotcycle(xlc,vlc,s,[size(xlc,1) 1 2]);



opt = contset(opt,'TestPath',[mfilename('fullpath'), '_run2']);
opt = contset(opt,'Filename','testadapt1_run2');
opt = contset(opt,'Singularities',               1);  
[x0,v0]=init_LC_LC_L(@adapt22,xlc,vlc,s(end),[xlc(end);1],[1 2],ntst,ncol);
[~, datafile]=contL(@limitcycleL,x0,v0,opt);
figure
axes
[xlc, vlc, ~] = loadPoint(datafile); % DV: load computed cycles
load('Data\testadapt1_run2.mat')       % DV: load singular points
plotcycle(xlc,vlc,s,[size(xlc,1) 1 2]);