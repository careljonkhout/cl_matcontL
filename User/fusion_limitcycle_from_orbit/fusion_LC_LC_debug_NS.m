run_init_if_needed
clc
format long
% continuation of limit cycles in fusion system
N = 25;                     
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
a = -1;
b = -0.3;
q_inf = -0.72;
poincare_tolerance = 1e-4;
handles = feval(odefile);

%% Continue limit cycle from orbit
p = cell2mat(parameters);
ntst = 20;
ncol = 4;
tolerance = 1e-2;

opt = contset();
opt = contset(opt, 'MaxNumPoints',   3);
opt = contset(opt, 'InitStepsize',   3);
opt = contset(opt, 'MinStepsize',    1e-12);
opt = contset(opt, 'MaxStepsize',    1);
opt = contset(opt, 'MaxNewtonIters', 1);
opt = contset(opt, 'MaxCorrIters',   10);
opt = contset(opt, 'VarTolerance',   1e-6);
opt = contset(opt, 'FunTolerance',   1e-8);
opt = contset(opt, 'Adapt',          3);
opt = contset(opt, 'MaxNumPoints',   400);
opt = contset(opt, 'CheckClosed',    50);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       true);
opt = contset(opt, 'Singularities',  true);
opt = contset(opt, 'CIS_UsingCIS',   false);
    % disable smoothing by angle:
opt = contset(opt, 'contL_SmoothingAngle', pi/2);
opt = contset(opt, 'contL_DiagnosticsLevel', 5);
opt = contset(opt, 'enable_bpc'            , true);
opt = contset(opt, 'console_output_level'  ,  5);
opt = contset(opt, 'contL_Testf_FunTolerance', 1e-5);
opt = contset(opt, 'contL_Testf_VarTolerance', 1e-4);
opt = contset(opt, 'bpc_tolerance', 5e-9);

opt = contset(opt, 'MaxTestIters',   20);


ap = 3;
filename = 'fusion_Orb_LC_with_bifurcations_06-Feb-2019_15_17_32';
datafile = fullfile('Data', [filename '.dat']);
matrix_file = fullfile('Data', filename);
load(matrix_file,'s');
singularities25 = s;
[x_from_file, v_from_file, ~, ~] = loadPoint(datafile);
singularity_index = 5;
singularity = singularities25(singularity_index);
i = singularity.index;
data = singularity.data;
par = data.parametervalues;
ap = data.ap; % active parameter


[x,v] = init_LC_LC_L(odefile, x_from_file, v_from_file, singularity, par, ap,ntst,ncol);

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


global lds
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

