function [t,y] = backward_euler(f,df,y0,t0,tf,n)
  h = (tf-t0)/n;
  t = linspace(t0, tf, n+1);
  y = zeros(n+1, length(y0));
  y(1,:) = y0;
  for i = 1:n
    x0 = y(i,:)';
    done = false;
    while ~ done
      A  = eye(length(y0)) - h/2 * feval(df,t(i),x0);
      x1 = x0 - A \ ( x0   - h/2 * feval( f,t(i),x0)- y(i,:)' );
      x0 = x1;
      done = norm(x1 - x0) < 1e-9;
    end
    y(i+1,:) = x0';
  end
end