load('arguments_for_multipliers', 'J')

global contopts
contopts = contset();
contopts.console_output_level = 3;
contopts.contL_DiagnosticsLevel = 0;

disp pqz
tic
multipliers_pqz_schur(J);
toc 




disp without_p_retrieve_J_first
tic
multipliers_retrieve_J_first(J);
toc


