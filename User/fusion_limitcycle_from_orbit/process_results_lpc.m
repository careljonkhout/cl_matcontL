run_init_if_needed
if ~ exist('x_lpc', 'var')
  filename = 'fusion_LPC_17-Feb-2019_9_40_03';
  datafile = fullfile('Data', [filename '.dat']);
  matrix_file = fullfile('Data', filename);
  load(matrix_file,'s');
  singularities_lpc = s;
  [x_lpc, ~] = loadPoint(datafile);
end


figure
hold on
plot(x_lpc(end,:),x_lpc(end-1,:),'b')
for singularity = singularities_lpc
  if ~(strcmp(singularity.label,'00') || strcmp(singularity.label,'99'))
    q_inf = singularity.data.x(end);
    a = singularity.data.x(end-1);
    text(q_inf, a, singularity.label)
  end
end
title(sprintf('fusion lpc b=-0.3'));
xlabel('q_{inf}')
ylabel('a')

figure
hold on
plot(x_lpc(end,:),x_lpc(end-2,:),'b')
for singularity = singularities_lpc
  if ~(strcmp(singularity.label,'00') || strcmp(singularity.label,'99'))
    q_inf = singularity.data.x(end);
    period = singularity.data.x(end-2);
    text(q_inf, period, singularity.label)
  end
end
title(sprintf('fusion lpc b=-0.3'));
xlabel('q_{inf}')
ylabel('period')


disp(' ')
disp('singularities detected for lpc:')
fprintf('index \t label \t q_inf \t\t a \t\t period \n')
for s = singularities_lpc
  fprintf('%d \t %s \t %.6f \t %.6f \t %.6f\n', ...
    s.index, s.label, s.data.x(end), s.data.x(end-1), s.data.x(end-2))
end
