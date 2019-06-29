function point = do_corrections(x0,v)
  
  global cds contopts
 
  curve_function_norm = max(abs(feval(cds.curve_func,x0)));
  lowest_curve_function_norm = curve_function_norm;
  x = x0;
  corrections = 0;
  done = false;
  while ~ done && corrections < contopts.MaxCorrIters
    period = x(end-1);
    if corrections >= 1
      print_diag(1, ...
        ['log_10 of function norm: %1.2f ' ...
         'log_10 of correction norm: %1.2f period: %.8e corrections: %d\n'], ...
          log10(curve_function_norm), ...
          log10(correction_norm), period, corrections);
    else
      print_diag(1, ...
       'log_10 of function norm: %1.2f period: %.8e corrections: %d\n', ...
        log10(curve_function_norm), period, corrections);
    end
    corrections = corrections + 1;
    old_x = x;
    if     isequal(cds.curve, @single_shooting)
      x = NewtonPicard.SingleShooting.do_one_correction(x0,x,v);
    elseif isequal(cds.curve, @multiple_shooting)
      x = NewtonPicard.MultipleShooting.do_one_correction(x0,x,v);
    else
      print_diag(0,'Newton_Picard_Corrections: wrong curvefile.\n');
      point = [];
      return
    end
   
    if isempty(x)
      break;
    end
    
    period = x(end-1);
    
    if period < 0
      print_diag(0, 'Period less than zero, abborting corrections\n')
      point = [];
      return;
    end
    
    correction_norm = max(abs(x - old_x));
    
    new_curve_function_norm = max(abs(feval(cds.curve_func,x)));
    
    if new_curve_function_norm > 20 * lowest_curve_function_norm
      print_diag(1, ...
        ['log_10 of function norm: %1.2f ' ...
         'log_10 of correction norm: %1.2f period: %.8e corrections: %d\n'], ...
          log10(curve_function_norm), ...
          log10(correction_norm), period, corrections);
      print_diag(0,['Curve function norm increased by a factor of 20.' ...
        ' Aborting Corrections.\n']);
      break;
    end
    lowest_curve_function_norm = min(lowest_curve_function_norm, ...
                                          new_curve_function_norm);
    curve_function_norm = new_curve_function_norm;
    
    done = curve_function_norm < contopts.FunTolerance && ...
               correction_norm < contopts.VarTolerance;
  end
  print_diag(1, ...
        ['log_10 of function norm: %1.2f ' ...
         'log_10 of correction norm: %1.2f period: %.8e corrections: %d\n'], ...
          log10(curve_function_norm), ...
          log10(correction_norm), period, corrections);
  if done
    point.R = curve_function_norm; 
    point.x = x;
    point.v = v;
    point.iters = corrections;
  else
    point = [];
    return
  end
end