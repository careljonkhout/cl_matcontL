function p_out = locate_NS(p1, p2, testfunctions)

  tstart = tic;
  global contopts
  
  p_out = [];
  
  MaxIters      = contopts.MaxTestIters;          
  VarTolerance  = contopts.contL_Testf_VarTolerance;


  [t1, failed1] = testfunctions(Constants.NS_id, p1.x, p1.v, 0);
  t1 = t1(Constants.NS_id);
  % todo: prevent computation of things have already been computed
  [t2, failed2] = testfunctions(Constants.NS_id, p2.x, p2.v, 0);
  t2 = t2(Constants.NS_id);
  print_diag(3,'t1:%d t2:%d\n',t1,t2);
  if ((~isempty(failed1)) || (~isempty(failed2))) && (failed1 || failed2)
    print_diag(3, 'Evaluation of testfunctions failed in bisection');
    return;
  end

  for i = 1:MaxIters
    % todo: perhaps it is possible to locate using secant method, instead of
    % simple bisection.
    x3 = p1.x + (p2.x-p1.x)/2;
    v3 = p1.v + (p2.v-p1.v)/2;
    v3 = v3/norm(v3);
    if contopts.NewtonPicard
      p3 = NewtonPicard.do_corrections(x3,v3);
      if isempty(p3)
        print_diag(3, 'NewtonPicard algorithm failed during bisection')
        return
      end
    else
      p3 = newtcorrL(x3,v3,1);
      if isempty(p3)
        print_diag(3, 'newtcorrL algorithm failed during bisection')
        return
      end
    end
    
    
    [tval, failed] = testfunctions(Constants.NS_id, p3.x, p3.v, 0);
    tval = tval(Constants.NS_id);
    if failed
      print_diag(3, 'Evaluation of testfunctions failed in bisection');
      return;
    end
    
    %JH: Changed to make the check for relative difference. 9/25/06
    %dist1 = norm(x-x1);
    %dist2 = norm(x-x2);
    dist1 = 2 * norm(p3.x-p1.x) / (norm(p1.x)+norm(p3.x));
    dist2 = 2 * norm(p3.x-p2.x) / (norm(p2.x)+norm(p3.x));
    
    % 
    if min(dist1,dist2) < VarTolerance
      failed2 = 0;
      p_out = p3;
      break;
    elseif tval == t1
      p1 = p3;
      % in case of 2 NS bifurcations we locate the one that occurs first along
      % the curve
    elseif min(t1,t2) <= tval && tval <= max(t1,t2)
      p2 = p3;
      if tval ~= t2
        warning('Two NS bifurcations detected. The second NS bifurcation along the curve will be ignored')
      end
    else
      print_diag(1, ...
        ['Error in locating NS in function locateNS in ' mfilename '\n']);
      return;
    end
  end

  if failed2
    print_diag(1, ['Maximum number of iterations reached ' ...
      ' without meeting tolerance in locateNS in ' mfilename ]);
  end

  print_diag(1,'Time spent in bisection: %f\n', toc(tstart));

end
