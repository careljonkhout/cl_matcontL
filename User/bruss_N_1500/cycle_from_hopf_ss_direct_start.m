% Continuation of cycles in the Brusselator 
% using single shooting with Newton Picard

% todo: fix initializer, so that continuation does not return to hopf

format long
format compact
clear global

odefile = @brusselator_1d;
% parameters for the system of ODEs:
N = 1500;
A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
L = pi * sqrt((Dx + Dy) / (B - A^2-1));
ode_parameters = [N L A B Dx Dy];



x0=[ones(N,1)*A; ones(N,1)*B/A; L];
h = 1e-2;
subspace_size = 10;


[x0, v0] = init_single_shooting_from_hopf(...
    odefile, x0, ode_parameters, active_parameter, h, subspace_size);

% we specify the options for the continuation of cycles using single shooting:
opts_h_lc = contset();
opts_h_lc = contset(opts_h_lc, ...
  ...
  'MaxNumPoints',           2, ...
  'InitStepsize',           0.1, ...
  'MaxStepsize',            0.1, ...
  'contL_SmoothingAngle',   100, ...
  'newtcorrL_use_max_norm', true, ...
  'Singularities',          true, ...
  'contL_DiagnosticsLevel', 5, ...
  'console_output_level',   5, ...
  'integration_abs_tol',    1e-12, ...
  'integration_rel_tol',    1e-12, ...
  'FunTolerance',           1e-6, ...
  'VarTolerance',           1e-6, ...
  'NewtonPicard',           true);

% we open a plot window:
figure
hold on
title_format_string = 'Brusselator N:%d  A:%.0f  B:%.1f  Dx:%.3f  Dy:%.3f';
title_format_args = {N,A,B,Dx,Dy};
xlabel('L')
ylabel('period')
title(sprintf(title_format_string, title_format_args{:}));

% we run the cycle continuation:
singularities = contL(@single_shooting, x0, v0, opts_h_lc, 'callback', ...
                                          @plot_T_versus_param);
                                        
% after the continuation has finished, we plot the singularities:
for singularity = singularities
  plot_sing(singularity)
end

function plot_sing(s, varargin)
  L = s.data.x(end);
  T = s.data.x(end-1); % period of cycle
  plot(L,T,'r*')
  text(L,T,s.label,varargin{:})
end



