load('continue_subspace_data.mat')
%NP_SS_continue_subspace(period,parameters);
cds.n_mesh_intervals = 1;
cds.orbits = cds.cycle_orbit;
contopts.NP_asisTolerance = 1e-3;
NP_MS_continue_subspaces(period, parameters);