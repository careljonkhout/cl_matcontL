
filename = 'fusion_Orb_LC_05-Dec-2018_15_01_51';
datafile = fullfile('Data', [filename '.dat']);
matrix_file = fullfile('Data', filename);
load(matrix_file);
[x, v, h, mult] = loadPoint(datafile);