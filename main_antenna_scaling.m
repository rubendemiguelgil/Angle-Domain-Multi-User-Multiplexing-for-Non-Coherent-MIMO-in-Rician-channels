clear, clc, close all;
% Determine where your m-file's folder is.
folder = fileparts(which('NC_simulation.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

params = simulation_parameters()
SNR_sweep = params.SNR_sweep;

currDate = strrep(datestr(datetime), ':', '_');
results_dir = ['Results/antenna_scaling'];
mkdir(results_dir)
save([results_dir '/params'], "params")

N_ant_array = [50 100 200 400];

SER_ant_scaling_nch_freq = zeros(length(N_ant_array), length(SNR_sweep));
SINR_ant_scaling_nch_freq = zeros(length(N_ant_array), length(SNR_sweep));

% SER_ant_scaling_perfect_mmse = zeros(length(N_ant_array), length(SNR_sweep));
% SINR_ant_scaling_perfect_mmse = zeros(length(N_ant_array), length(SNR_sweep));

SER_ant_scaling_imperfect_mmse = zeros(length(N_ant_array), length(SNR_sweep));
SINR_ant_scaling_imperfect_mmse = zeros(length(N_ant_array), length(SNR_sweep));

for i = 1:length(N_ant_array)
params.M = N_ant_array(i);


% NCH_freq
params.diff_decoding_dimension = 'freq';
results_nch_freq = simulate_noncoherent(params)
SER_ant_scaling_nch_freq(i,:) = results_nch_freq.SER_total_mtx;
SINR_ant_scaling_nch_freq(i,:) = results_nch_freq.SINR_total_mtx;


% MMSE
% params.combining_method = 'MMSE';
% params.perfect_channel_estimation = true;
% results_ch_perfect_mmse = simulate_coherent(params)
% SER_ant_scaling_perfect_mmse(i,:) = results_ch_perfect_mmse.SER_total_mtx;
% SINR_ant_scaling_perfect_mmse(i,:) = results_ch_perfect_mmse.SINR_total_mtx;

params.perfect_channel_estimation = false;
results_ch_imperfect_mmse = simulate_coherent(params)
SER_ant_scaling_imperfect_mmse(i,:) = results_ch_imperfect_mmse.SER_total_mtx;
SINR_ant_scaling_imperfect_mmse(i,:) = results_ch_imperfect_mmse.SINR_total_mtx;

end

save([results_dir '/workspace'])