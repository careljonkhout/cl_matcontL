classdef Assert
  methods(Static)
    
    function greater_than(n, name, value)
      assert(value > n, [name ' must be greater than ' num2str(n) '.']);
    end
    
    function less_than(n, name, value)
      assert(value < n, [name ' must be less than ' num2str(n) '.']);
    end
    
    function positive(name, value)
      assert(value > 0 , [name ' must be strictly positive.']);
    end
    
    function scalar(name, value)
      assert(numel(value) == 1, [name ' must be a scalar.']);
    end
    
    function integer(name, value)
      assert(round(value) == value, [name ' must be an integer.']);      
    end
    
    function function_handle(name, value)
      assert(isa(value, 'function_handle'), [name ' must be a function handle.']);
    end
    
  end
end

