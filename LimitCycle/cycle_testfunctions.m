function out = cycle_testfunctions(ids_testf_requested, multipliers, v)
  global cds
  const = Constants;
  print_diag(6,'evaluating testfunctions with multipliers:\n')
  print_diag(6,'%s\n', multipliers2str(multipliers));
  
  
  if any(ismember(const.BPC_id, ids_testf_requested))
    distance_to_one = abs(multipliers - 1);
    [~,index]       = min(distance_to_one);
    multipliers(index) = 0;
    out(const.BPC_id) = real(prod(multipliers - ones(size(multipliers))));
  end
  if ismember(const.PD_id, ids_testf_requested)
    % real is needed to ensure that small complex parts induced by roundoff
    % errors do not make the result complex. A complex valued test function
    % would cause false positives when detecting bifurcations.
    out(const.PD_id) = real(prod(multipliers + ones(size(multipliers))));
  end
  if ismember(const.LPC_id, ids_testf_requested)
    out(const.LPC_id) = v(end);
  end
  if any(ismember(const.NS_id, ids_testf_requested))   
    threshold              = cds.deviation_of_trivial_multiplier;
    complex_mults          = multipliers(abs(imag(multipliers)) > threshold);
    moduli                 = abs(complex_mults(1:2:end));
    out(const.NS_id)       = prod(1-moduli);
  end
end

