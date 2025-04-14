clear, clc, close all;
% Determine where your m-file's folder is.
folder = fileparts(which('NC_simulation.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

params = simulation_parameters()
SNR_sweep = params.SNR_sweep;

currDate = strrep(datestr(datetime), ':', '_');
results_dir = ['Results/def_antenna_scaling'];% currDate];
mkdir(results_dir)
save([results_dir '/params'], "params")

N_ant_array = [50 100 200 400 800]

% SER_ant_scaling_leg_single_user = zeros(length(N_ant_array), length(SNR_sweep));
% SINR_ant_scaling_leg_single_user = zeros(length(N_ant_array), length(SNR_sweep));
% 
% SER_ant_scaling_nch_eep = zeros(length(N_ant_array),length(SNR_sweep));
% SINR_ant_scaling_nch_eep = zeros(length(N_ant_array), length(SNR_sweep));

% SER_ant_scaling_nch_freq = zeros(length(N_ant_array), length(SNR_sweep));
% SINR_ant_scaling_nch_freq = zeros(length(N_ant_array), length(SNR_sweep));

SER_ant_scaling_perfect_mmse = zeros(length(N_ant_array), length(SNR_sweep));
SINR_ant_scaling_perfect_mmse = zeros(length(N_ant_array), length(SNR_sweep));

SER_ant_scaling_imperfect_mmse = zeros(length(N_ant_array), length(SNR_sweep));
SINR_ant_scaling_imperfect_mmse = zeros(length(N_ant_array), length(SNR_sweep));

for i = 1:length(N_ant_array)
params.M = N_ant_array(i);

% % Legacy NCH_freq Single user
% params.N_users = 1; 
% params.user_angles = [deg2rad(30)];
% params.user_pwr = [1];
% results_nch_leg_single_user = simulate_noncoherent_joint_const(params)
% SER_ant_scaling_leg_single_user(i,:) = results_nch_leg_single_user.SER_total_mtx;
% SINR_ant_scaling_leg_single_user(i,:) = results_nch_leg_single_user.SINR_total_mtx;


% % Legacy NCH_freq joint EEP const (2 users)
% params.N_users = 2; 
% params.user_angles = [deg2rad(30) deg2rad(-30)];
% params.user_pwr = [1 1];
% results_nch_eep = simulate_noncoherent_joint_const(params)
% SER_ant_scaling_nch_eep(i,:) = results_nch_eep.SER_total_mtx;
% SINR_ant_scaling_nch_eep(i,:) = results_nch_eep.SINR_total_mtx;


% % NCH_freq
% params.diff_decoding_dimension = 'freq';
% results_nch_freq = simulate_noncoherent(params)
% SER_ant_scaling_nch_freq(i,:) = results_nch_freq.SER_total_mtx;
% SINR_ant_scaling_nch_freq(i,:) = results_nch_freq.SINR_total_mtx;


% MMSE
params.combining_method = 'MMSE';
params.perfect_channel_estimation = true;
results_ch_perfect_mmse = simulate_coherent(params)
SER_ant_scaling_perfect_mmse(i,:) = results_ch_perfect_mmse.SER_total_mtx;
SINR_ant_scaling_perfect_mmse(i,:) = results_ch_perfect_mmse.SINR_total_mtx;

params.perfect_channel_estimation = false;
results_ch_imperfect_mmse = simulate_coherent(params)
SER_ant_scaling_imperfect_mmse(i,:) = results_ch_imperfect_mmse.SER_total_mtx;
SINR_ant_scaling_imperfect_mmse(i,:) = results_ch_imperfect_mmse.SINR_total_mtx;

end

save([results_dir '/workspace'])