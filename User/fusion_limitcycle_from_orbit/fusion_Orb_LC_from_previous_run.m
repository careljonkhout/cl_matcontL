run_init_if_needed
clc
format long
global lds cds contopts %#ok<NUSED>
clear('lds')
clear('cds')
clear('contopts')
% continuation of limit cycles in fusion system
N = 25;                     
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));

opt = contset();
opt = contset(opt, 'MaxNumPoints',   5);
opt = contset(opt, 'InitStepsize',   1e-1);
opt = contset(opt, 'MinStepsize',    1e-12);
opt = contset(opt, 'MaxStepsize',    1e-1);
opt = contset(opt, 'MaxNewtonIters', 1);
opt = contset(opt, 'MaxCorrIters',   10);
opt = contset(opt, 'VarTolerance',   1e-6);
opt = contset(opt, 'FunTolerance',   1e-8);
opt = contset(opt, 'Adapt',          3);
opt = contset(opt, 'MaxNumPoints',   1000);
opt = contset(opt, 'CheckClosed',    50);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       false);
opt = contset(opt, 'Singularities',  true);
opt = contset(opt, 'CIS_UsingCIS',   false);
opt = contset(opt, 'Locators'    , [0 0 0 0]);
    % disable smoothing by angle:
opt = contset(opt, 'contL_SmoothingAngle', pi/2);
opt = contset(opt, 'contL_DiagnosticsLevel', 5);
opt = contset(opt, 'enable_bpc'            , true);
opt = contset(opt, 'enable_nf_lpc'         , false);
opt = contset(opt, 'enable_nf_pd'          , false);
opt = contset(opt, 'console_output_level'  ,  5);
opt = contset(opt, 'contL_Testf_FunTolerance', 1e-5);
opt = contset(opt, 'contL_Testf_VarTolerance', 1e-4);
opt = contset(opt, 'bpc_tolerance', Inf);
opt = contset(opt, 'newtcorrL_use_max_norm', true);
opt = contset(opt, 'MaxTestIters',   30);
opt = contset(opt, 'SingularTestFunction',   true);
opt = contset(opt, 'MoorePenrose',          false);

% right before alleged bpc

filename = 'fusion_Orb_LC_from_previous_run_12-Feb-2019_14_52_13';
%filename = 'fusion_Orb_LC_with_bifurcations_12-Feb-2019_14_25_42';
mat_filename = fullfile('Data',[filename '.mat']);
dat_filename = fullfile('Data',[filename '.dat']);
load(mat_filename,'s');
s = s(end);
[x,v] = loadPoint(dat_filename);
par = s(end).data.parametervalues;
ap  = s(end).data.ap;
ntst = 80; % s(end).data.ntst;
ncol = s(end).data.ncol;

[x,v] = init_LC_LC_L(odefile, x, v, s, par, ap,ntst,ncol);

[sout, datafile] = contL(@limitcycleL,x,v,opt); %#ok<*ASGLU>


figure
hold on;
coordinate1 = 1;
coordinate2 = 2;
 title(sprintf( ...
    'fusion N:%d a:%.2f n:%.2f q_inf:%.2f', ...
     N,a,b,q_inf));
 xlabel('n_1')
 ylabel('U_1')
[x, v, h, mult] = loadPoint(datafile); % DV: load computed cycles



close all
figure
hold on
for i=1:10:opt.MaxNumPoints
  xx = x(1:end-1-length(lds.ActiveParams),i);
  coord1_vals = xx(coordinate1:lds.nphase:end);
  coord2_vals = xx(coordinate2:lds.nphase:end);
  plot(coord1_vals,coord2_vals,'b');
  title(sprintf( ...
    'fusion N:%d a:%.2f n:%.2f q_{inf}:varied', ...
     N,a,b));
  xlabel('n_1')
  ylabel('U_1')
end
figure
plot(x(end,:),x(end-1,:))
title(sprintf( ...
    'fusion N:%d a:%.2f n:%.2f q_{inf}:varied', ...
     N,a,b));
xlabel('q_{inf}')
ylabel('period')

figure
hold on;
nMults = size(mult,1);
for i=nMults-10:nMults
  plot(x(end,:),mult(i, :))
end

