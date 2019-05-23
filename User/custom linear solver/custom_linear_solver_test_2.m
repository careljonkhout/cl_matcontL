
close all
load('arguments_for_multipliers.mat')
global lds

lds.ntst = 5;
lds.nphase = 3;

ntst = lds.ntst;
ncol = lds.ncol;
nphase = lds.nphase;

JJ = zeros((ntst*ncol+1)*nphase +2);
for i=1:lds.ntst
  for j=1:lds.ncol+1
    for k = 1:lds.ncol
      col_indices = (1:nphase) + (ncol*(i-1) + j-1) * nphase;
      row_indices = (1:nphase) + (ncol*(i-1) + k-1) * nphase;
      
      JJ(row_indices, col_indices) = rand(1) * eye(nphase);
    end
  end
end
JJ(1:(ntst*ncol*nphase),end-1:end) = rand(ntst*ncol*nphase,2);
JJ(end-1:end,:) = rand(2,size(JJ,2));
JJ(end-nphase-1:end-2, 1:nphase          ) =   eye(nphase); 
JJ(end-nphase-1:end-2, end-nphase-1:end-2) = - eye(nphase);
b = rand((lds.ntst*lds.ncol+1)*lds.nphase+2,1);

tic; x = linear_solver_collocation_2(JJ,b); toc
tic; x_std = JJ \ b; toc

hold on
abs_diff = abs(x-x_std);
x = abs(x);
x_std = abs(x_std);



fprintf('%.8f\n',abs_diff)
bb = JJ*x;

