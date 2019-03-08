
#include<math.h>
#include<mex.h>
#include<matrix.h>

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  
/* q = size(J,1)-1; */
J = J(1:q,1:q);
p = speye(q);
r = lds.col_coords;
number_of_nonzeros_p = lds.ntst*numel(lds.col_coords)^2;
p_row_indices = zeros(number_of_nonzeros_p,1);
p_col_indices = zeros(number_of_nonzeros_p,1);
p_values = zeros(number_of_nonzeros_p,1);
p_element_counter = 0;

for (size_t i=1;i<ntst;i++) {
  /*   sJ = J(r,lds.nphase+r); */
  /*[sl,~] = lu(sJ); */
  /* m = inv(sl); */
  for (size_t j=1;j<ncol_coord;j++) {
  }
  
  for j=r
    range = p_element_counter + (1:numel(r));
    p_row_indices(range) = j * ones(numel(r),1);
    p_col_indices(range) = r;
    p_values(range) = m(j-r(1)+1,:);
    p_element_counter = p_element_counter + numel(r);
  end
  r = r+lds.ncol_coord;
end*
