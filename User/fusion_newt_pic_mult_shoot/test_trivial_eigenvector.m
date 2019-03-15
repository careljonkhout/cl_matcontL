trivial_eigenvector = cds.dydt_ode(0,x(1:end-2),parameters{:}); 

Mv = d_phi_d_x(x(1:end-2),trivial_eigenvector,period,parameters);

MV./trivial_eigenvector