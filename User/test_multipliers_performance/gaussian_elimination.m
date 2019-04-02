N = 5;
A = rand(N,N+1);
clc
for i=1:N
  A(i,:) = A(i,:) / A(i,i);
  for j=i+1:N
    A(j,:) = A(j,:) - A(j,i) * A(i,:);
  end
end


  


