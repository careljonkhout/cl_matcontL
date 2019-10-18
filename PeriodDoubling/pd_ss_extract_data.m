

function [phases_0, v, period, parameters] = pd_ss_extract_data(x)
  global cds
  phases_0                     = x(1:cds.nphases);
  v                            = x(cds.nphases+1:2* cds.nphases);
  period                       = x(end-2);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = x(end-1:end);
  parameters                   = num2cell(parameters);
  