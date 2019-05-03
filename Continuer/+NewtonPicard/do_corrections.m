function point = do_corrections(x0,v)
  
  global cds contopts
 
  curve_function_norm = max(abs(feval(cds.curve_func,x0)));

  if curve_function_norm > 10
    point = [];
    return;
  end
  
  x = x0;
  corrections = 0;
  done = false;
  while ~ done && corrections < contopts.MaxCorrIters
    period = x(end-1);
    fprintf('function_norm: %.8e period: %.8e corrections: %d\n', ...
                            curve_function_norm, period, corrections);
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
    correction_norm = max(abs(x - old_x));
    if isempty(x)
      break;
    end

    if period < 0
      fprintf('period less than zero, correction aborted\n');
      break
    end
    new_curve_function_norm = max(abs(feval(cds.curve_func,x)));
    if new_curve_function_norm > 2 * curve_function_norm
      print_diag(0,'function_norm: %.8e period: %.8e corrections: %d\n', ...
        new_curve_function_norm, period, corrections);
      print_diag(0,['Curve function norm is strongly increasing.' ...
        ' Aborting Corrections.\n']);
      break;
    end
    curve_function_norm = new_curve_function_norm;
    
    done = curve_function_norm < contopts.FunTolerance && ...
               correction_norm < contopts.VarTolerance;
  end
  if curve_function_norm < contopts.FunTolerance
    point.R = curve_function_norm; 
    point.x = x;
    point.v = v;
    point.iters = corrections;
  else
    point = [];
    return
  end
end