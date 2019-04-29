function loadOdeFile(curve_struct, odefile)

  odefile_handles = feval(odefile);

  curve_struct.odefile   = odefile; % todo: carefully remove this duplicate data (the other copy is cds.probfile)
  curve_struct.func      = odefile_handles{2};
  curve_struct.Jacobian  = odefile_handles{3};
  curve_struct.JacobianP = odefile_handles{4};
  curve_struct.Hessians  = odefile_handles{5};
  curve_struct.HessiansP = odefile_handles{6};
  curve_struct.Der3      = odefile_handles{7};
  curve_struct.Der4      = odefile_handles{8};
  curve_struct.Der5      = odefile_handles{9};


  if length(odefile_handles) > 10      
      curve_struct.user = odefile_handles(10:end);
  else
      curve_struct.user=[];
  end
end