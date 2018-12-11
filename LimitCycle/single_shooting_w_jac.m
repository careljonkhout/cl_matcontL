function out = single_shooting
%
% Curve file of cycle continuation with single shooting
%
    out{1}  = @curve_func;
    out{2}  = @defaultprocessor;
    out{3}  = @options;
    out{4}  = [];%@jacobian;
    out{5}  = [];%@hessians;
    out{6}  = [];%@testf;
    out{7}  = [];%@userf;
    out{8}  = [];%@process;
    out{9}  = [];%@singmat;
    out{10} = [];%@locate;
    out{11} = @init;
    out{12} = [];%@done;
    out{13} = @adapt;
    out{14} = @curve_CIS_first_point;
    out{15} = [];%@curve_CIS_step;
%---------------------------------------------------------  
function func = curve_func(varargin)
  global cds
  x = varargin{1};
  p0 = x(end);
  T  = x(end-1);
  x0 = x(1:end-2);
  L = 1.1; A = 1; B = 2.2; Dx = p0; Dy = 0.004;
  parameters = {L; A; B; Dx; Dy};
  trajectory = shoot(x0, T, parameters);
  func = [trajectory(end,:)'-x0;sum(sum(cds.x0_prime.*trajectory))]; 
%---------------------
function trajectory = shoot(x, T, parameters)
  global cds
  odefile = @brusselator_N_2;
  handles = feval(odefile);
  dydt = handles{2};
  f =@(t, y) dydt(t, y, parameters{:});
  integration_opt = odeset(...
    'AbsTol',      1e-10,    ...
    'RelTol',      1e-10,    ...
    'BDF',         'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1,      ...
    'Jacobian',     @(t,y) feval(handles{3},t,y,parameters{:}) ...
  );
  [~, trajectory] = ode15s(f, ...
    linspace(0, T, cds.nDiscretizationPoints), x, integration_opt);
  
  
  
function [value, isterminal, direction] = returnToPlane(t, x, x0, v0)
   % v0 should be a row vector
   % x and x0 should be column vectors
   value = v0*(x-x0);
   isterminal =  sum((x-x0).^2) < poincare_tolerance && t > 1;
   direction = 1;

%---------------------------------------------------------
function init(~,~)
%----------------------------------------------------------
function out=defaultprocessor(varargin)
  out = varargin{1};

function options

function out = curve_CIS_first_point(~)
  out = 1;

function adapt(varargin)
