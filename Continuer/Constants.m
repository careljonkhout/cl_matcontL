classdef Constants

  
  properties( Constant = true )
    % bifurcation ids for limitcycles:
    BPC_id = 1
    PD_id  = 2;
    LPC_id = 3;
    NS_id  = 4;
    
    % testfunction behaviours:
    sign_change   = 0;
    sign_constant = 1;
    value_change  = 2;
    ignore        = 8;
  end
 
end

