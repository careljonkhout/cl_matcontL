function j = cjac_densely_stored(odefile,jacobian,x,p,ap)
  j = full(cjac(odefile,jacobian,x,p,ap));
end