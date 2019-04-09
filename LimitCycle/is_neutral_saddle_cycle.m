function is_neutral_saddle_cycle = is_neutral_saddle_cycle(singularity_index)
  global cds lds contopts
  if singularity_index ~= 4 
    % if the Neimark-Sacker testfunction did not change sign, there is no
    % neutral saddle cycle.
    is_neutral_saddle_cycle = false;
    return;
  end
  if isequal(cds.curve, @limitcycleL)
    mults = lds.multipliers;
    mults = mults(abs(mults) > 0.1);
    mults = mults(abs(imag(mults)) > contopts.real_v_complex_threshold);
    is_neutral_saddle_cycle = isempty(mults);
  elseif isequal(cds.curve, @single_shooting) || ...
         isequal(cds.curve, @multiple_shooting)
    mults = cds.multipliers;
    mults = mults(abs(mults) > 0.1);
    mults = mults(abs(imag(mults)) > contopts.real_v_complex_threshold );
    is_neutral_saddle_cycle = isempty(mults);
  else
    is_neutral_saddle_cycle = false;
  end
end