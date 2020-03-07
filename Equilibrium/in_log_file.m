function in_log_file(eval_r,eval_l)
print_diag(2,'\nEigenvalues: (Subspace)\n')
for j=1:length(eval_r)
    print_diag(2,'%+e',real(eval_r(j)))
    if imag(eval_r(j))
        print_diag(2,' %+e i',imag(eval_r(j)));
    end
    print_diag(2,'\n')
end
print_diag(2,'Eigenvalues: (Reference)\n')
for j=1:length(eval_l)
    print_diag(2,'%+e',real(eval_l(j)))
    if imag(eval_l(j))
        print_diag(2,' %+e i',imag(eval_l(j)));
    end
    print_diag(2,'\n')
end

