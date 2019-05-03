function loadOdeFile_into_lds(odefile)

  global lds;

  odefile_handles = feval(odefile);

  lds.odefile   = odefile; % todo: carefully remove this duplicate data (the other copy is cds.probfile)
  lds.func      = odefile_handles{2};
  lds.Jacobian  = odefile_handles{3};
  lds.JacobianP = odefile_handles{4};
  lds.Hessians  = odefile_handles{5};
  lds.HessiansP = odefile_handles{6};
  lds.Der3      = odefile_handles{7};
  lds.Der4      = odefile_handles{8};
  lds.Der5      = odefile_handles{9};


  if length(odefile_handles) > 10      
      lds.user = odefile_handles(10:end);
  else
      lds.user=[];
  end
end