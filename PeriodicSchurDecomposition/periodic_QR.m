% seed random number generator for reproducible output
rng(1);
% relatively though rng(2) m = 6 N=100;
format long
clc;
clear all; %#ok<CLALL>
m=6; % number of matrices
N=100; % size of matrix

orthogonal_transformations_tolerance = eps(m*N^2);

G = zeros([N N m]); 
for i=1:m
  G(:,:,i) = randi(2,N) + rand(N);
end

[HH,~]=reduce_to_hessenberg(G);

  prod = eye(N);
  for i=m:-1:1
    prod = HH(:,:,i) * prod;
  end



Q = zeros([N N m]);
for i=1:m
  Q(:,:,i) = eye(N);
end

HH(abs(HH)<eps(100))=0;

M=HH;
iteration_counter = 1;
L = 1; 
U = N-1;
is_deflated = zeros(m,N-1);
tolerance = eps;

iteration_block_size = 15;
total_iteration_counter = 0;
while U-L >= 1 %&& total_iteration_counter < 50
  iteration_counter=0;
  while smallest_element_on_subdiagonal(M,L,U) > tolerance
    total_iteration_counter = total_iteration_counter + 1;
    iteration_counter = iteration_counter + 1; 
    previous_r = smallest_element_on_subdiagonal(M,L,U);
    prod = eye(U-L+2);
    for i=m:-1:1
      prod = M(L:U+1,L:U+1,i) * prod;
    end
    if mod(iteration_counter,iteration_block_size) == iteration_block_size / 4
      dat1 = 0.75;
      dat2 = -0.4375;
      H = prod(end-2:end,end-2:end);
      alpha = abs(H(2,1)) + abs(H(3,2));
      a = dat1*alpha;
      b = dat2*alpha;
      c = alpha;
      d = H(1,1);
      [M,Q]=double_shift(M,Q,L,U+1,a,b,c,d); 
      print_and_pause(M,L,U,previous_r,'DLAQR special',iteration_counter)
    elseif mod(iteration_counter,iteration_block_size) == 2*iteration_block_size / 4
      H = prod(end-2:end,end-2:end);
      a = 1.5 * (abs(H(2,1)) + abs(H(2,2)));
      b = 0;
      c = 0;
      d = 1.5 * (abs(H(3,2)) + abs(H(3,3)));
      [M,Q]=double_shift(M,Q,L,U+1,a,b,c,d); 
      print_and_pause(M,L,U,previous_r,'DHSEQE special',iteration_counter)
    %elseif false && mod(iteration_counter,iteration_block_size) == 30

    else
      if U-L >= 7
        a = M(end-1,end-1,m);
        b = M(end-1,end  ,m);
        c = M(end  ,end-1,m);
        d = M(end  ,end  ,m);
        if (a-d)^2 + 4*b*c < 0
          [M,Q] = double_shift(M,Q,L,U+1,a,b,c,d);
          print_and_pause(M,L,U,previous_r,'ds mod 9',iteration_counter)
        else
          [M,Q] = single_shift(M,Q,d,L,U+1);
          print_and_pause(M,L,U,previous_r,'ss mod 9',iteration_counter)
        end
      else
        
        H = prod(end-3:end,end-3:end);
        eigenvalues = eigs(H);
        [~,i] = min(abs(eigenvalues));
        if isreal(eigenvalues(i))
          [M,Q] = single_shift(M,Q,eigenvalues(i),L,U+1);
          print_and_pause(M,L,U,previous_r,'ss ev',iteration_counter)
        else
          a = real(eigenvalues(i));
          b = imag(eigenvalues(i));
          c = -imag(eigenvalues(i));
          d = real(eigenvalues(i));
          [M,Q] = double_shift(M,Q,L,U+1,a,b,c,d);
          print_and_pause(M,L,U,previous_r,'ds ev',iteration_counter)
        end
      end
    end
    if false && (previous_r < smallest_element_on_subdiagonal(M,L,U))
      dat1 = 0.75;
      dat2 = -0.4375;
      H = M(end-2:end,end-2:end,m);
      alpha = abs(H(2,1)) + abs(H(3,2));
      a = dat1*alpha;
      b = dat2*alpha;
      c = alpha;
      d = H(1,1);
      [M,Q]=double_shift(M,Q,L,U+1,a,b,c,d); 
      print_and_pause(M,L,U,previous_r,'DLAQR after bad shift',iteration_counter)
    end
  end

    
    M(abs(M)<tolerance) = 0;

   
  [M,Q,is_deflated] = deflate_zeros(M,Q,is_deflated);
  

  {'before resetting L U', L, U, iteration_counter}
  M(abs(M)<tolerance) = 0;
  [L,U] = find_largest_unreduced_block(M(:,:,m));
  {'end of main loop', L, U}
end


maximum_error = 0;
for i=2:m
  assert(...
    all(abs(Q(:,:,i)' * HH(:,:,i) * Q(:,:,i-1)-M(:,:,i)) ...
     < orthogonal_transformations_tolerance,'all'),...
    'Assertion failed: Q_%d * H_%d * Q_%d - M_%d exceeds tolerance', ...
    i, i, i-1, i);
  maximum_error_for_H_i = ...
    max(max(abs(Q(:,:,i)' * HH(:,:,i) * Q(:,:,i-1)-M(:,:,i))));
  maximum_error = max(maximum_error_for_H_i, maximum_error); 
end
% note that Q(:,:,0) is defined as Q(:,:,m)
assert(all(abs(Q(:,:,1)' * HH(:,:,1) * Q(:,:,m)-M(:,:,1)) ...
  < orthogonal_transformations_tolerance,'all'), ...
  'Assertion failed: Q_1 * H_1 * Q_m - M_1 exceeds tolerance');
maximum_error_for_H_1 = max(max(abs(Q(:,:,1)' * HH(:,:,1) * Q(:,:,m)-M(:,:,1))));
maximum_error = max(maximum_error_for_H_1, maximum_error); 
fprintf(['The maximum error in the Periodic Schur Decomposition' ...
  ' in this test run is %.5e.\n'], maximum_error);

for i=1:m
  orthogonal = max(max(Q(:,:,i)*Q(:,:,i)'-eye(N))) ...
    < orthogonal_transformations_tolerance;
  assert(orthogonal,['Assertion failed: Q_%d is not orthogonal' ...
    ' .max(max(Q(:,:,i)*Q(:,:,i)-eye(N)))=%.5f'], ...
    i, max(max(Q(:,:,i)*Q(:,:,i)'-eye(N))));
end

% set entries that are very small to zero
% we hope that this way all entries that should be zero become zero
M(abs(M)<eps(N*N))=0;

for i=1:m-1
  % check if M_1 up to M_{m-1} are upper triangular
  assert(istriu(M(:,:,i)),'Assertion failed: M_%d is not upper triangular',i);
end
% check if M_m is Hessenberg
assert(bandwidth(M(:,:,m),'lower') == 1, ...
  'Assertion failed: M_m is not Hessenberg');

function smallest = smallest_element_on_subdiagonal(M,L,U)
  smallest = Inf;
  for i=L:U
    smallest = min([smallest abs(M(i+1,i,end))]);
  end
end

function [L_max,U_max] = find_largest_unreduced_block(M)
  L=1;
  N = size(M,1);
  largest_block_size = 1;
  L_max = 1;
  U_max = 1;
  subdiagonal = diag(M,-1);
  while L < N-2
    while L < N-2 && (~subdiagonal(L) || ~ subdiagonal(L+1))
      L = L + 1;
    end
    if L == N - 2 && ~(subdiagonal(L) && subdiagonal(L+1))
      % periodic QR algorithm is finished
      break
    else 
      % find U
      U = L+1;
      while U < N-1 && subdiagonal(U+1)
        U = U + 1;
      end
    end
    if U - L > largest_block_size
      L_max = L;
      U_max = U;
      largest_block_size = U - L;
    end
    L = U + 1;
  end
end

function print_and_pause(M,L,U,previous_r,shift_type,iteration_counter)
  if iteration_counter*(U-L) > 100*100
    {L,U,smallest_element_on_subdiagonal(M,L,U) ,  smallest_element_on_subdiagonal(M,L,U) - previous_r, shift_type}
    pause
  end

end
