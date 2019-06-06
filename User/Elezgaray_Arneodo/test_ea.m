handles = feval(@elezgaray_arneodo);
dydt = handles{2};
global N;
N = 20;

y = [ -2 * ones(N,1) ; -4 * ones(N,1)];

D = 0.14;
eps = 0.01;
alpha = 0.01;



dy = feval(dydt,0,y,D,eps,alpha);

[t,y] = ode15s(@(t,y) dydt(t,y,D,eps,alpha),[0 100], y);


equilibrium = y(end,:)';

[x0,v0] = init_EP_EP_L(@elezgaray_arneodo, equilibrium, [D,eps,alpha], 1);

opt = contset( ...
  'MaxNumPoints', 50, ...
  'Backward', true, ...
  'MaxStepsize', 1, ...
  'Singularities', true, ...
  'CIS_UsingCIS', true, ...
  'contL_DiagnosticsLevel', 1, ...
  'console_output_level',   1 ...
);

  


% 2   H  :  +3.230753e-02     2.081301e+01     9.832189e-08     3.452756e-02
s=contL(@equilibriumL,x0,v0,opt);

hopf = s(2);

active_parameter = 1;

x = hopf.data.x;
ode_parameters = [D, eps, alpha];
ode_parameters(active_parameter) = hopf.data.x(end);

% h will be the amplitude of the initial cycle
h = 0.001;
% direction of the parameter change
dp = -1;

% ntst is the number of mesh intervals. That is, the cycle will be represented
% using a piecewise polynomial function of ntst pieces.
ntst = 20;
% ncol is the number of collocation points in each mesh interval
ncol = 4;

% We run the initializer for continuation of cycles by collocation from a Hopf
% point:
[x0, v0] = init_collocation_from_hopf(...
    @elezgaray_arneodo, x, ode_parameters, active_parameter, h, dp, ntst, ncol);

% We specify the option for the cycle continuation.
opts_h_lc = contset( ...
  ...
  'MaxNumPoints',          208, ...
  'InitStepsize',           0.05, ...
  'MaxStepsize',            0.05, ...
  'contL_SmoothingAngle',   100, ...
  'newtcorrL_use_max_norm', true, ...
  'Singularities',          true, ...
  'enable_nf_ns',           false, ...
  'enable_nf_lpc',          false, ...
  'enable_nf_pd',           false, ...
  'contL_DiagnosticsLevel', 0, ...
  'console_output_level',   0);
figure
hold on

xlabel('D')
ylabel('period')
title('Elezgaray-Arneodo');

% We run the cycle continuation:
singularities = contL(@limitcycleL, x0, v0, opts_h_lc, 'callback', ...
                                          @plot_T_versus_param);

% We plot the singularities from the cycle continuation
vertical_alignment = 'top';
for singularity = singularities
  plot_sing(singularity, 'VerticalAlignment', vertical_alignment)
  % switch vertical alignment after each singularity to reduce overlap of labels
  if strcmp(vertical_alignment, 'top')
    vertical_alignment = 'bottom';
  else
    vertical_alignment = 'top';
  end
end

function plot_sing(s, varargin)
  L = s.data.x(end);
  T = s.data.x(end-1);
  plot(L,T,'r*')
  text(L,T,s.label,varargin{:})
end

