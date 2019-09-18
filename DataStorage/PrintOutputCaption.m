function PrintOutputCaption
  global cds

  print_diag(0,'             :  ')
  for i = 1:cds.nap
      print_diag(0,'p(%d)              ', cds.ActiveParams(i))
  end

  if has_period(cds.curve)
    caption = ['period           point norm       curve func norm' ...
               '  step size       mult. accuracy  largest mult\n'];
    print_diag(0,caption)
  else
  % example of output for reference:
  %              :  p(2)            norm of point 	 curve function norm 	 step size 
  %   1   1   00 :  +5.000000e-01     1.882005e+01     0.000000e+00     1.000000e-02
  %   2          :  +5.100000e-01     1.882005e+01     0.000000e+00     1.000000e-02
  %   3   2   H  :  +5.128157e-01     1.882005e+01     0.000000e+00     2.815666e-03
    print_diag(0,'point norm       curve func norm  step size\n')
  end
end