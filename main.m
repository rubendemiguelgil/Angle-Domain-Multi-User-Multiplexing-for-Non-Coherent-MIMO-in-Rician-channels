clear, clc, close all;
% Determine where your m-file's folder is.
folder = fileparts(which('NC_simulation.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

params = simulation_parameters()

results_nch = simulate_noncoherent(params)


currDate = strrep(datestr(datetime), ':', '_');
results_dir = ['Results/' currDate ];
mkdir(results_dir)
save([results_dir '/results_nch'], "results_nch")