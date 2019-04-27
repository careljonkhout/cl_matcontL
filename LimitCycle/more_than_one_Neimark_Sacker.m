function more_than_1_NS = more_than_one_Neimark_Sacker(tvals1,tvals2)
  global cds;
  % Neimark-Sacker bifurcation only occur in continuation of limit cycles
  % hence we test if the current curve is actualy a limit cycle
  if ~ ( isequal(cds.curve, @single_shooting) || ...
         isequal(cds.curve, @multiple_shooting) || ...
         isequal(cds.curve, @limitcycleL) )
    more_than_1_NS = false;
    return;
  else
    more_than_1_NS = abs(tvals1(Constants.NS_id) - tvals2(Constants.NS_id)) > 2;
  end
end