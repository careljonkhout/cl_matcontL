function v = find_tangent_vector(curvefile, x)
  if     isequal(curvefile,@single_shooting)
    v = NewtonPicard.SingleShooting.find_tangent_vector(x);
  elseif isequal(curvefile, @multiple_shooting)
    v = NewtonPicard.MultipleShooting.find_tangent_vector(x);
  else
    error_format_string = ['internal error in' ...
      ' Continuer/+NewtonPicard/find_tangent_vector.m:\n' ...
          'The Newton-Picard method does not support %s'];
    error(error_format_string, func2str(curvefile));
  end
end