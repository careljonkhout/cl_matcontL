
filename = 'fusion_Orb_LC_05-Dec-2018_15_01_51';
datafile = ['Data/' filename '.dat'];
matrix_file = ['Data/' filename];
load(matrix_file);
[x, v, h, mult] = loadPoint(datafile);