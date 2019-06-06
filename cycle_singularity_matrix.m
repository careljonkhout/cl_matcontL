% Defines which changes in testfunctions correspond to which singularity type
% for cycle continuations. The columns of matrix correspond to testfunctions and
% the rows correspond to singularities BPC, PD, LPC, and NS respectively.
function [matrix, labels] = cycle_singularity_matrix()

  sign_change    = Constants.sign_change;
  sign_constant  = Constants.sign_constant;
  value_change   = Constants.value_change;
  o              = Constants.ignore; 

  matrix = [ 
  % BPC_testfunc PD_testfunc  LPC_testfunc   NS_testfunc
    sign_change  o            sign_constant  o            % BPC    
    o            sign_change  o              o            % PD   
    o            o            sign_change    o            % LPC   
    o            o            o              value_change % NS
  ];
  labels = [ 'BPC';'PD '; 'LPC'; 'NS ' ];
  
  assert(size(matrix,1) == length(labels))
end


