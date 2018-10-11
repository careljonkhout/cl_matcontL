% 
function [Index1, Index2] = hopf_locator_construction(border, two_index)
%this is for minumally augmented system of locating Hopf point
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');
N = size(border,2);
p = 2;
N = N - p;

A = border(1:N,1:N);

last_p_row = border(N+1:end,1:N);
last_p_column = border(1:N,N+1:end);
block = border(N+1:end,N+1:end);

new_last_p_row = (A'\last_p_row')';
new_block = block - new_last_p_row*last_p_column;


cinf = inf;
I = 0;

for k=1:6
    
    index1 = two_index(k,1);
    index2 = two_index(k,2);

    new_last_2_row = [new_last_p_row(index1,:);new_last_p_row(index2,:)];
    new_2_block = [new_block(index1,:);new_block(index2,:)];

    part_matrix = [new_last_2_row, new_2_block];
    part_matrix_inverse = -new_2_block\[new_last_2_row, -eye(p)];

    number1 = zeros(1+p,1);
    number2 = zeros(1+p,1);

    number1(end) = 1;
    number2(end) = 1;

    for i=1:p
        number1(i) = sum(abs(part_matrix(i,:)));
        number2(i) = sum(abs(part_matrix_inverse(i,:)));
    end    

    cond1 = max(number1);
    cond2 = max(number2);

    if cinf > cond1*cond2
        cinf = cond1*cond2;
        I = k;
    end

end

Index1 = two_index(I,1);
Index2 = two_index(I,2);
    
    