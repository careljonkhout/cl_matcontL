function gauss_elim_Jacobian_fun(blocks, ntst, ncol, nphase)
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
end

