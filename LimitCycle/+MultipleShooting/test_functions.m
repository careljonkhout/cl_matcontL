% test functions are used for detecting AND location singularities by bisection
% when detecting ids will be cds.ActTest, and when locating ids will contain
% only those ids relevant to the bifurcation that is being located.
function [out, failed] = test_functions(ids, x0, v, ~) 
  % unused arguments are v and CISdata
  global cds
  
  failed = false;
  
  PD_id  = 6; % id for period doubling test function
  LPC_id = 7; % id for limit point of cycles test function
  NS_id  = 8; % id for Neimarck-Sacker test function
  
  if any(ismember([PD_id NS_id],ids))
    update_multipliers_if_needed(x0)
  end
  if any(ismember(1:5,ids))
    % detection of Branching points of cycles is currently not implemented
    out(1:5) = ones(5,1);
  end
  if ismember(PD_id, ids)
    % real is needed to ensure that small complex parts induced by roundoff
    % errors do not make the result complex.
    % A complex valued test function would cause false positives when detecting 
    % bifurcations.
    out(6) = real(prod(cds.multipliers + ones(size(cds.multipliers))));
  end
  if ismember(LPC_id, ids)
    out(7) = v(end);
  end
  if ismember(NS_id, ids)
    mults = cds.multipliers;
    psi_ns = 1;
    for i = 1 : (length(mults) - 1)
      for j = (i + 1) : length(mults)
        psi_ns = psi_ns * (real(mults(i) * mults(j)) - 1);
      end
    end
    out(8) = psi_ns;
  end
end