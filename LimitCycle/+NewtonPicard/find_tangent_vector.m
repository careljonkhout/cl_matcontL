function v = find_tangent_vector(curvefile, x, v0)

  switch func2str(curvefile)
    case 'single_shooting'
      v = NewtonPicard.SingleShooting.find_tangent_vector(x);
    case 'multiple_shooting'
      v = NewtonPicard.MultipleShooting.find_tangent_vector(x);
    case 'perioddoubling_ss'
      handles  = feval(curvefile);
      jacobian = feval(handles{4}, x);
      kernel   = null(jacobian);
      v        = kernel(:,1);
    otherwise
    error_format_string = ['internal error in' ...
      ' Continuer/+NewtonPicard/find_tangent_vector.m:\n' ...
          'The Newton-Picard method does not support %s'];
    error(error_format_string, func2str(curvefile));
  end
  % v0' * v is the inner product of v and v0, which equals the cosine of the
  % angle between v and v0. If the cosine of the angle is less than zero, then
  % the angle between v and v0 is more than 90 degrees. Therefore we flip the
  % signs of the components of v if the cosine is less than zero. We depend on
  % the consistent orientation of the tangent vector when detecting limit points
  % of cycles (LPC).
  if nargin == 3 && v0'*v < 0
    v = -v;
  end
  v = v / max(abs(v));
end