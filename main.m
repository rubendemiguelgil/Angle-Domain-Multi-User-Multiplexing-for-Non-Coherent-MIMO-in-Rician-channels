clear, clc, close all;
% Determine where your m-file's folder is.
folder = fileparts(which('main.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

params = simulation_parameters()

currDate = strrep(datestr(datetime), ':', '_');
results_dir = ['Results/' currDate];
mkdir(results_dir)
save([results_dir '/params'], "params")

% Legacy NCH_freq joint EEP const (2 users)
results_nch_eep = simulate_noncoherent_joint_const(params)
save([results_dir '/results_nch_eep'], "results_nch_eep")
disp('Nch joint const done')

% NCH_freq
params.diff_decoding_dimension = 'freq';
results_nch_freq = simulate_noncoherent(params)
save([results_dir '/results_nch_freq'], "results_nch_freq")
disp('Nch freq done')


% MMSE
params.perfect_channel_estimation = false;
results_ch_imperfect_mmse = simulate_coherent(params)
save([results_dir '/results_mmse_imperf'], "results_ch_imperfect_mmse")
disp('MMSE imperf done')