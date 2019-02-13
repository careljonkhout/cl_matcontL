
filename = 'fusion_Orb_LC_with_bifurcations_06-Feb-2019_15_17_32';
datafile = fullfile('Data', [filename '.dat']);
matrix_file = fullfile('Data', filename);
load(matrix_file,'s');
singularities25 = s;
[x25, v25, h25, mult25] = loadPoint(datafile);
figure
hold on
plot(x25(end,45:end),x25(end-1,45:end),'b')
for singularity = singularities25
  if ~(strcmp(singularity.label,'00') || strcmp(singularity.label,'99'))
    q_inf = singularity.data.x(end);
    period = singularity.data.x(end-1);
    text(q_inf, period, singularity.label)
  end
end
title(sprintf('fusion N=25 a:%.2f b:%.2f', -1,-0.3));
xlabel('q_{inf}')
ylabel('period')


disp(' ')
disp('singularities detected for N=25:')
fprintf('index \t label \t q_inf \t\t period\n')
for s = singularities25
  fprintf('%d \t %s \t %.6f \t %.6f\n', ...
    s.index, s.label, s.data.x(end), s.data.x(end-1))
end
