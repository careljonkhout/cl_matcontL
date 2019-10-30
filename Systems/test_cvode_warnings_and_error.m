% Calls cvode with a IVP that goes backwards in time, which cvode cannot solve,
% to test if error messages from cvode are displayed correctly in the matlab
% command window.

N = 10;
brusselator_1d_N_10.recompile_cvode_mex
global contopts
contopts = contset('console_output_level', 5);
[t, y] = brusselator_1d_N_10.cvode( ...
  'initial_point',   - ones(2 * N,1), ...
  'ode_parameters',    ones(5,1), ...
  't_values',          [ 0 -1 ]);
  