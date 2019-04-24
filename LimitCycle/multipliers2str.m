function str = multipliers2str(multipliers)
  global contopts cds
  
  distance_to_one                     = abs(multipliers - 1);
  cds.deviation_of_trivial_multiplier = min(distance_to_one);
  print_diag(0,'deviation of trivial multiplier: %.2e\n', ... 
                          cds.deviation_of_trivial_multiplier);
  
  if contopts.NewtonPicard
    print_diag(0, 'number of multipliers: %d\n', length(multipliers));
  end
  
  % Printing every multiplier that is computed would produce too much output,
  % therefore we print only the multipliers of which the norm is above a
  % threshold that can be configured by the user.
  threshold = contopts.multiplier_print_threshold;
  multipliers = multipliers(abs(multipliers) > threshold);
  str = sprintf('multipliers with norm larger than %.3f:\n', threshold);
  i=1;
  while i <= length(multipliers)
    m = multipliers(i);
    if abs(imag(m)) > contopts.real_v_complex_threshold
      % This text will not be so long that it will negatively impact
      % performance. Therefore, the 'array grows on every loop
      % iteration'-warning is ignored
      str  = [ str sprintf('%.15f +/- %.15fi norm: %.15f\n', ...
                real(m), abs(imag(m)), abs(m))]; %#ok<AGROW>
	    i = i + 2;
    else
      str  = [ str sprintf('%.15f\n', real(m)) ]; %#ok<AGROW>
      i = i + 1;
    end
  end
end