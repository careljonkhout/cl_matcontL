function point = NP_do_corrections(x0,v)
  
  global cds contopts
  curve_function_norm        = max(abs(feval(cds.curve_func,x0)));
  correction_norm            = 0;
  lowest_curve_function_norm = curve_function_norm;
  x                          = x0;
  corrections                = 0;
  done                       = false;
  while ~ done && corrections < contopts.MaxCorrIters
    print_stats(curve_function_norm, correction_norm, corrections)
    corrections = corrections + 1;
    old_x = x;
    switch func2str(cds.curve)
      case 'single_shooting'
        x = NP_SS_do_one_correction(x0,x,v);
      case 'multiple_shooting'
        x = NP_MS_do_one_correction(x0,x,v);
      case 'perioddoubling_ss'
        x = NP_SS_do_one_correction_pd(x0,x,v);
      otherwise
        print_diag(0,'Newton_Picard_Corrections: wrong curvefile.\n');
        point = [];
        return
    end
    if isempty(x)
      break;
    end

    correction_norm = max(abs(x - old_x));

    curve_function_norm = max(abs(feval(cds.curve_func,x)));

    if curve_function_norm > ...
               contopts.max_rel_funcnorm_increase * lowest_curve_function_norm
      print_diag(1,[ ...
        'Current curve function norm is now more than %d times ' ...
        'the lowest curve function norm that ' ...
        'was attained in this continuation step. ' ...
        'Aborting Corrections.\n'], contopts.max_rel_funcnorm_increase);
      break;
    end

    lowest_curve_function_norm = min(lowest_curve_function_norm, ...
                                            curve_function_norm);

    done = curve_function_norm < contopts.FunTolerance && ...
               correction_norm < contopts.VarTolerance;
  end
  print_stats(curve_function_norm, correction_norm, corrections);
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

function print_stats(curve_function_norm, correction_norm, corrections)
  if corrections >= 1
    print_diag(1, ...
      'function norm: 10^%1.2f correction norm: 10^%1.2f corrections: %d\n', ...
      log10(curve_function_norm), log10(correction_norm), corrections);
  else
    print_diag(1, ...
      'function norm: 10^%1.2f corrections: %d\n', ...
      log10(curve_function_norm), corrections);
  end
end