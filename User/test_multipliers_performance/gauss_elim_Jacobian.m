close all
load('arguments_for_multipliers.mat')

global lds

ntst   = lds.ntst;
ncol   = lds.ncol;
nphase = lds.nphase;

blocks = zeros(nphase*(ncol+1)+1,nphase*(ncol+1)+3,ntst);
j_row_indices = 1 :  ncol    * nphase;
j_col_indices = 1 : (ncol+1) * nphase;
for i = 1 : ntst
  blocks(1:ncol*nphase, 1:((ncol+1)*nphase),i) = ...
    J(j_row_indices, j_col_indices);
  blocks(1:ncol*nphase, end-2:end-1,i) = ...
           J(j_row_indices, end-1:end);
	blocks(end,1:(ncol +1)*nphase) = J(end,j_col_indices);

  j_row_indices = j_row_indices + ncol*nphase;
  j_col_indices = j_col_indices + ncol*nphase;
end



tic
% Gauss elim on block 1 and 2
for i=1:nphase
  blocks(ncol * nphase+i,                 i, 1   ) =  1;
  blocks(ncol * nphase+i, ncol * nphase + i, ntst) = -1;
end

for block_index = 1:ntst
  for i = 1:ncol*nphase
    blocks(i,:,block_index) = blocks(i,:,block_index) / ...
                                blocks(i,i,block_index);
    for j = i+1:((ncol+1) * nphase + 1)
      if blocks(j,i,block_index)
        blocks(j,:,block_index) = blocks(j,:,block_index) ...
          - blocks(j,i,block_index) * blocks(i,:,block_index);
      end
    end
  end
  if block_index < ntst
    indices = (ncol*nphase)+(1:nphase);
    blocks(     indices, 1:nphase, block_index+1) ...
      =  blocks(indices, indices , block_index);
  end
end

toc
