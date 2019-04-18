function str = multipliers2str(multipliers)
  global contopts
  
  distance_to_one = abs(multipliers - 1);
  accuracy        = min(distance_to_one);
  print_diag(0,'deviation of trivial multiplier: %.2e\n', accuracy);
  
  if contopts.NewtonPicard
    print_diag(0, 'number of multipliers: %d\n', length(multipliers));
  end
  
  % Printing every multiplier that is computed would produce too much output,
  % therefore we print only the multipliers of which the norm is above a
  % threshold that can be configured by the user.
  threshold = contopts.multiplier_print_threshold;
  multipliers = multipliers(abs(multipliers) > threshold);
  str = sprintf('multipliers with norm larger than %.3f:\n', threshold);
  
  for i=1:length(multipliers)
    m = multipliers(i);
    if abs(imag(m)) > contopts.real_v_complex_threshold
      if  sign(imag(m)) > 0
        sign_text = '+';
      else
        sign_text = '-';
      end
      % This text will not be so long that it will negatively impact
      % performance. Therefore, the 'array grows on every loop
      % iteration'-warning is ignored
      str  = [ str sprintf('%.15f %s %.15fi norm: %.15f\n', ...
                real(m), sign_text, abs(imag(m)), abs(m))]; %#ok<AGROW>
    else
      str  = [ str sprintf('%.15f\n', real(m)) ]; %#ok<AGROW>
    end
  end
end