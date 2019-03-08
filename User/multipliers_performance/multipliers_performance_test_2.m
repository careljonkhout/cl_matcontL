load 'arguments_for_multipliers'

global contopts
contopts = contset();
contopts.console_output_level = 0;
contopts.contL_DiagnosticsLevel = 0;

multipliers_pqz_schur_2(J);