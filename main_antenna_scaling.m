clear, clc, close all;
% Determine where your m-file's folder is.
folder = fileparts(which('main_antenna_scaling.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

params = simulation_parameters()
SNR_sweep = params.SNR_sweep;

% Environment
params.N_users = 6; 
params.W = 5;
params.user_angles = [deg2rad(50) deg2rad(-50) deg2rad(30) deg2rad(-30) deg2rad(10) deg2rad(-10)];
params.user_pwr = [1 1 1 1 1 1];
assert(length(params.user_angles) == params.N_users, 'There must be one angle per user.')
assert(length(params.user_pwr) == params.N_users, 'There must be one power per user.')


currDate = strrep(datestr(datetime), ':', '_');
results_dir = ['Results/antenna_scaling'];
mkdir(results_dir)
save([results_dir '/params'], "params")

N_ant_array = [50 100 200 400];

SER_ant_scaling_nch_freq = zeros(length(N_ant_array), length(SNR_sweep));
SINR_ant_scaling_nch_freq = zeros(length(N_ant_array), length(SNR_sweep));

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
params.perfect_channel_estimation = false;
results_ch_imperfect_mmse = simulate_coherent(params)
SER_ant_scaling_imperfect_mmse(i,:) = results_ch_imperfect_mmse.SER_total_mtx;
SINR_ant_scaling_imperfect_mmse(i,:) = results_ch_imperfect_mmse.SINR_total_mtx;

save([results_dir '/workspace'])

end

save([results_dir '/workspace'])
