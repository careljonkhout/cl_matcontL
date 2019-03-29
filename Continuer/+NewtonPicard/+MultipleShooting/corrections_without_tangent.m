function x = corrections_without_tangent(x)
  
  global cds contopts
 
  curve_function_norm = max(abs(feval(cds.curve_func,x)));
  corrections = 0;
  while (~ (curve_function_norm < contopts.FunTolerance)) ...
      && corrections < contopts.MaxCorrIters
    period = x(end-1);
    fprintf('function_norm: %.8e period: %.8e corrections: %d\n', ...
      curve_function_norm, period, corrections);
    corrections = corrections + 1;
 
    x = NewtonPicard.MultipleShooting.one_correction_without_tangent(x);
    
    if period < 0
      fprintf('period less than zero, correction aborted\n');
      break
    end
    new_curve_function_norm = max(abs(feval(cds.curve_func,x)));
    if (new_curve_function_norm > 2 * curve_function_norm) 
      print_diag(0,'function_norm: %.8e period: %.8e corrections: %d\n', ...
        new_curve_function_norm, period, corrections);
      print_diag(0,['Curve function norm is strongly increasing.' ...
        ' Aborting Corrections.\n']);
      break;
    end
    curve_function_norm = new_curve_function_norm;
  end
end