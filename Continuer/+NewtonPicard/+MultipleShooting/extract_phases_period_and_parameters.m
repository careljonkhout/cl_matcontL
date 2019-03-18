% extracts 
% - y ( the current approximation of points on the cycle )
% - period
% - parameters ( of the ode system in which cycles are continued )
% from the continuation state vector x.
% The non-active parameters, i.e. the parameters that do not change during
% the continuation are extracted from the global struct cds
% (i.e.) curve description structure.
% The parameters are returned as a cell array, so that that can be passed to
% cds.dydt_ode in an syntactically elegant manner.

function [phases,period,parameters] = extract_phases_period_and_parameters(x)
  global cds
  phases                       = x(1:cds.nphases*cds.nShootingPoints);
  phases                       = reshape(phases,cds.nphases,cds.nShootingPoints);
  period                       = x(end-1);
  parameter_value              = x(end);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = parameter_value;
  parameters                   = num2cell(parameters);
end