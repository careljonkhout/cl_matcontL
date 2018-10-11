% 
function [cinf] = partial_cond(border, p)
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

if size(border,1) ~= size(border,2)
   error('the border matrix is NOT square\n'); 
end    

N = size(border,1);
N = N - p;

A = border(1:N,1:N);

last_p_row = border(N+1:end,1:N);
last_p_column = border(1:N,N+1:end);
block = border(N+1:end,N+1:end);

new_last_p_row = A'\last_p_row';
new_last_p_row = new_last_p_row';
new_block = block - last_p_row*last_p_column;

part_matrix = [new_last_p_row, new_block];
part_matrix_inverse = [-block\new_last_p_row, block\eye(p)];

number1 = zeros(1+p,1);
number2 = zeros(1+p,1);

number1(end) = 1;
number2(end) = 1;

for i=1:p
    number1(i) = sum(abs(part_matrix(i,:)));
    number2(i) = sum(abs(part_matrix_inverse(i,:)));
end    

%number1
%number2

cond1 = max(number1);
cond2 = max(number2);

cinf = cond1*cond2;
    
    