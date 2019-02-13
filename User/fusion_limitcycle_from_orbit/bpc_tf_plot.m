
filename = 'fusion_Orb_LC_from_previous_run_13-Feb-2019_14_36_16';
datafile = fullfile('Data', [filename '.dat']);
matrix_file = fullfile('Data', filename);
%load(matrix_file,'s');
%singularities25 = s;
[x, v, h, tvals] = loadPoint(datafile);
figure
plot(tvals(1,30:end-20),'*-')
xlabel('Continuation point')
ylabel('testfunction 1')


