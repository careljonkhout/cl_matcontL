
filename = 'fusion_Orb_LC_N_75_23-Jan-2019_19_04_12';
datafile = fullfile('Data', [filename '.dat']);
matrix_file = fullfile('Data', filename);
load(matrix_file);
[x, v, h, mult] = loadPoint(datafile);