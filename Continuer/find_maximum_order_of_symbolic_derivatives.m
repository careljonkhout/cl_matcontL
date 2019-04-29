function [max_order, max_order_params] = ...
                             find_maximum_order_of_symbolic_derivatives(odefile)
  max_order = 0; 
  max_order_params = 0;
  
  odefile_handles = feval(odefile);

  if     ~isempty(odefile_handles{9}),   max_order = 5; 
  elseif ~isempty(odefile_handles{8}),   max_order = 4; 
  elseif ~isempty(odefile_handles{7}),   max_order = 3; 
  elseif ~isempty(odefile_handles{5}),   max_order = 2; 
  elseif ~isempty(odefile_handles{3}),   max_order = 1; 
  end
  if     ~isempty(odefile_handles{6}),   max_order_params = 2; 
  elseif ~isempty(odefile_handles{4}),   max_order_params = 1; 
  end

end

