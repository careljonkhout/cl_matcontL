
%cd \Users\Carel\Downloads\cl_matcontL_2018_8_27\CL_MATCONTL\
%init
%cd User\hopf_normal_form\
my_equilibrium = zeros(20,1);

[x0,v0]=init_EP_EP_L(@hopf_normal_L,my_equilibrium,[0],[1]);
options=contset;options=contset(options,'Singularities',1);
options = contset(options,'Filename','hopf_normal');
[sout,datafile]=contL(@equilibriumL,x0,[],options);


[x0,v0] = init_H_LC_L(@hopf_normal_L,my_equilibrium,[0],[1],1,20,4);
options=contset;
options.MaxNumPoints=10;
options.InitStepsize=0.5;
options.MinStepsize=0.5;
options.MaxStepsize=0.5;
[sout,datafile]=contL(@limitcycleL,x0,v0,options);

[xlc, vlc, slc] = loadPoint(datafile); % DV: load computed cycles
load('Data\hopf_normal.mat')            % DV: load singular points

% function [x,v] = init_LC_LC(odefile, x, v, s, par, ap, ntst, ncol)
last_parameter_value=xlc(end:end);
[x0,v0] = init_LC_LC_L(@hopf_normal_L,xlc,vlc,sout(end),...
                [last_parameter_value],[1]);
%[xlc,vlc,slc,hlc,flc]=cont(@limitcycleL,x0,v0,options);
plotcycle(xlc,vlc,s,[1 2]);
[sout,datafile]=contL(@limitcycleL,x0,v0,options);
%figure
%plotcycle(xlcc,vlcc,slcc,[1 2]);