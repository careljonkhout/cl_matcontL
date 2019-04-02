load('arguments_for_multipliers', 'J')

global contopts
contopts = contset();
contopts.console_output_level = 3;
contopts.contL_DiagnosticsLevel = 0;

tic
v = find_initial_tangent_vector(J);
toc

tic

jacobian_solve(J,v)



function v0 = find_initial_tangent_vector(J)
  v = zeros(size(J,2),1);
  v(end) = 1;
  B = [J; v'];
  C = zeros(size(J,2),1);
  C(end) = -1;
  v0 = bordCIS1(B,C,1);
  v0 = v0/norm(v0);
end
