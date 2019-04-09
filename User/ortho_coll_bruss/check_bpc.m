load('s')
handles = limitcycleL;
jac_handle = handles{4};
jac = feval(jac_handle,cds.sout(3).data.x);

