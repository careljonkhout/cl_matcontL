function str = multipliers2str(multipliers)
  global contopts cds
  
  distance_to_one                     = abs(multipliers - 1);
  cds.deviation_of_trivial_multiplier = min(distance_to_one);
  log_10_of_deviation                 = log10(min(distance_to_one));
  str = sprintf('deviation of trivial multiplier: 10^%.2f\n', ... 
                                                           log_10_of_deviation);
  
  if contopts.NewtonPicard
    str = [str sprintf('number of multipliers: %d\n', length(multipliers))];
  end
  
  % Printing every multiplier that is computed would produce too much output,
  % therefore we print only the multipliers of which the norm is above a
  % threshold that can be configured by the user.
  threshold = contopts.multiplier_print_threshold;
  multipliers = multipliers(abs(multipliers) > threshold);
  str = [str sprintf('multipliers with norm larger than %.3f:\n', threshold)];
  i=1;
  while i <= length(multipliers)
    m = multipliers(i);
    if abs(imag(m)) > cds.deviation_of_trivial_multiplier
      % This text will not be so long that it will negatively impact
      % performance. Therefore, the 'array grows on every loop
      % iteration'-warning is ignored
      str  = [ str sprintf('%.6f +/- %.6fi norm: %.6f\n', ...
                real(m), abs(imag(m)), abs(m))]; %#ok<AGROW>
	    i = i + 2;
    else
      str  = [ str sprintf('%.6f\n', real(m)) ]; %#ok<AGROW>
      i = i + 1;
    end
  end
end