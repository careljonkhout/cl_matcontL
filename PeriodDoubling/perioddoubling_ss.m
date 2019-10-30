function out = perioddoubling_ss
% period doubling with single shooting curve definition file
  out{1}  = @curve_function;
  out{2}  = @default_processor;
  out{3}  = @options;
  out{4}  = @jacobian;
  out{5}  = [];%@hessians;
  out{6}  = @testf;
  out{7}  = @userf;
  out{8}  = @process;
  out{9}  = @singmat;
  out{10} = @locate;
  out{11} = @init;
  out{12} = @done;
  out{13} = @adapt;
  out{14} = @CIS_first_point;
  out{15} = @CIS_step;
end

function f = curve_function(varargin)
  global cds
  [phases_0, v, period, parameters] = pd_ss_extract_data(varargin{1});
  shoot                        = @NewtonPicard.shoot;
  phases_end                   = shoot(phases_0, period, parameters);
  cds.phases_0                 = phases_0;
  M                            = @NewtonPicard.SingleShooting.monodromy_map;
  Mv                           = M(v, period, parameters);
  f = [phases_end - phases_0; 
      (phases_0 - cds.previous_phases)' * cds.previous_dydt_0;
      Mv + v;
      cds.l' * v - 1];
end

function J = jacobian(varargin)
  global contopts
  x = varargin{1};
  x1 = x;
  x2 = x;
  J = zeros(length(x)-1, length(x));
  for i = 1:length(x)
    x1(i) = x1(i) - contopts.Increment;
    x2(i) = x2(i) + contopts.Increment;
    J(:,i) = curve_function(x2) - curve_function(x1);
    x1(i) = x(i);
    x2(i) = x(i);
  end
  J = J/(2*contopts.Increment);
end

function point = default_processor(varargin)
  point = varargin{1};
  savePoint(point);
end

function options(varargin)
end

function [values, failed] = testf(varargin)
  values = 1;
  failed = false;
end

function userf(varargin)
end

function process(varargin)
  x = varargin{1};

end

function [S,L] = singmat(varargin)
  S = [];%[ 0 8 0 8;
        %8 0 8 8
        %1 1 0 8
        %1 1 8 0];
  L = [];%['R2  ';'LPPD';'GPD ';'PDNS'];
end

function locate(varargin)
end

function init(varargin)
end

function done(varargin)
end

function [has_changed, x_out, v_out, CISdata] = adapt(x, v, CISdata, ~)
  x_out = x;
  v_out = v;
  has_changed = false;
end

function out = CIS_first_point(varargin)
  out = 1;
end

function out = CIS_step(varargin)
  out = 1;
end
