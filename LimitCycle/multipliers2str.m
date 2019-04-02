function str = multipliers2str(multipliers)
  str = '';
  for i=1:length(multipliers)
    m = multipliers(i);
    if ~ isreal(m)
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