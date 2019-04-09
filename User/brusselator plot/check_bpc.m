
filename_lc_lc_1 = 'ortho_coll_bruss_08-Apr-2019_14_18_32';
filename = filename_lc_lc_1;
datafile = fullfile([filename '.dat']);
matrix_file = filename;
load(matrix_file,'s');
singularities_lc_lc_1 = s;
x_lc_lc_1 = loadPoint(datafile);

filename_bpc_lc_1 = 'bruss_bpc_bpc_oc_08-Apr-2019_15_13_38';
filename = filename_bpc_lc_1;
datafile = fullfile([filename '.dat']);
matrix_file = filename;
load(matrix_file,'s');
singularities25 = s;
x_bpc_lc_1 = loadPoint(datafile);

handles = limitcycleL;
jac = feval(handles{4});
feval(jac,singularities_lc_lc_1{7).data);

