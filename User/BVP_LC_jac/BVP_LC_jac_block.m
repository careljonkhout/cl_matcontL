function block = BVP_LC_jac_block(i, x, nphase, ncol, mesh, fine_mesh)
  global lds cds;
    
  gl_positions                 = gl_pos(ncol);
  period                       = x(end-1);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = x(end);
  parameters                   = num2cell(parameters);
  block                        = zeros(nphase * ncol, nphase * (ncol + 1));
  for j = 1:(ncol + 1)
    col_indices = (j-1) * nphase + (1:nphase);
    for k = 1:ncol
      row_indices = (k-1) * nphase + (1:nphase);
        block(row_indices, col_indices) =  A(i,j,k);
    end
  end


  function out = A(i,j,k)
    jac = cjac(lds.func, lds.Jacobian, ...
                p(z(i,k),i), parameters, cds.ActiveParams);
    out = eval_lagrange(z(i,k),i,j) * period * jac;
    out = out - eye(size(jac)) * eval_langrange_prime(z(i,k),i,j);
  end
  
  function out = eval_lagrange(t, i, j)
    out = 1;
    for l = 1:(ncol+1)
      if l ~= j
        out = out * (t - s(i,l)) / (s(i,j) - s(i,l));
      end
    end
  end

  % see: https://math.stackexchange.com/questions/1105160/evaluate-derivative-of-lagrange-polynomials-at-construction-points
  function out = eval_langrange_prime(t, i, j)
    sum = 0;
    for l = 1:(ncol+1)
      if l ~= i
        sum = sum + 1/(t - s(i,l));
      end
    end
    out = sum * eval_lagrange(t, i, j);
  end

  function s = s(i,j)
    s = fine_mesh((i-1) * ncol + j);
  end

  function z = z(i,k)
    delta_t = mesh(i+1) - mesh(i);
    z = mesh(i) + delta_t * (gl_positions(k) + 1) / 2;
  end

  function sum = p(t,i)
    sum = zeros(nphase,1);
    for jp = 1:ncol
      x_indices = (1:nphase) + ((i-1)*ncol + jp) * nphase;
      sum = sum + x(x_indices) * eval_lagrange(t,i,jp);
    end
  end
end
% gl_pos