% After a singularity is detected and located, the contL calls this function.
% Usually a normal form coefficient (nfc) is computed, via a function call from
% this function. However, for continuation of cycles by single shooting or
% multiple shooting, the computation of normal form coefficients is not
% implemented yet.
function [failed,s] = shooting_process_singularity(id,point,s)
  global cds
  x = point.x;
  switch id
  case 1
    format_string = 'Branch Point cycle(period = %e, parameter = %e)\n'; 
    print_diag(0, format_string, x(end-1), x(end));
    s.msg  = sprintf('Branch Point cycle'); 
  case 2
    format_string = 'Period Doubling (period = %e, parameter = %e)\n';
    print_diag(0, format_string, x(end-1), x(end));
    s.msg  = 'Period Doubling';
  case 3
    s.msg = 'Limit point cycle';
    format_string = 'Limit point cycle (period = %e, parameter = %e)\n';
    print_diag(0, format_string, x(end-1), x(end));
  case 4
    d = cds.multipliers;
    smallest = Inf;
    % we find the pair of multipliers whose product is closest to one, and check
    % to see these the multipliers have a nonzero imaginary part.
    for jk=1:length(d)-1
      [val,idx] = min(abs(d(jk+1:length(d))*d(jk)-1));
      if val < smallest
        idx2 = jk+idx;
        smallest = val;
      end
    end
    threshold = cds.deviation_of_trivial_multiplier;
    singularity_is_neutral_saddle = abs(imag(d(idx2))) < threshold;
    if singularity_is_neutral_saddle
      s.msg = 'Neutral saddle cycle';
      format_string = 'Neutral Saddle Cycle (period = %e, parameter = %e)\n';
      % A neutral saddle is not really a bifurcation
      print_diag(0, format_string, x(end-1), x(end));
    else
      s.msg = 'Neimark Sacker';
      format_string = 'Neimark-Sacker (period = %e, parameter = %e)\n';
      print_diag(0, format_string, x(end-1) ,x(end));
      %print_diag(0, 'Normal form coefficient = %d\n', s.data.nscoefficient);
    end
  end
  failed = 0;
end  