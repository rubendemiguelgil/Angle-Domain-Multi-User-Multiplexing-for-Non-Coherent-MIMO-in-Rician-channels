clear, clc, close all;
% Determine where your m-file's folder is.
folder = fileparts(which('NC_simulation.m')); 
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

% MRC
params.combining_method = 'MRC';
params.perfect_channel_estimation = true;
results_ch_perfect_mrc = simulate_coherent(params)
save([results_dir '/results_mrc_perf'], "results_ch_perfect_mrc")
disp('MRC perf done')

params.perfect_channel_estimation = false;
results_ch_imperfect_mrc = simulate_coherent(params)
save([results_dir '/results_mrc_imperf'], "results_ch_imperfect_mrc")
disp('MRC imperf done')

% ZF
params.combining_method = 'ZF';
params.perfect_channel_estimation = true;
results_ch_perfect_zf = simulate_coherent(params)
save([results_dir '/results_zf_perf'], "results_ch_perfect_zf")
disp('ZF perf done')

params.perfect_channel_estimation = false;
results_ch_imperfect_zf = simulate_coherent(params)
save([results_dir '/results_zf_imperf'], "results_ch_imperfect_zf")
disp('ZF imperf done')

% MMSE
params.combining_method = 'MMSE';
params.perfect_channel_estimation = true;
results_ch_perfect_mmse = simulate_coherent(params)
save([results_dir '/results_mmse_perf'], "results_ch_perfect_mmse")
disp('MMSE perf done')

params.perfect_channel_estimation = false;
results_ch_imperfect_mmse = simulate_coherent(params)
save([results_dir '/results_mmse_imperf'], "results_ch_imperfect_mmse")
disp('MMSE imperf done')