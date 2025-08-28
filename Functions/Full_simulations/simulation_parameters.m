function [params] = simulation_parameters()
% Creates an object with the simulation parameters
%% Parameters
    % Simulation
    params.L = 102400; % Tx length (in bits)
    params.bps = 2; % 2 bits/symbol in QPSK
    params.L_sym = params.L/params.bps; % Tx length in syms
    params.CP_length = 128; % For now didnt include it due to the narrowband assumption and the use of BER and SINR as KPIs
    params.N_subcarriers = 1024 - 128; % Number of useful subcarriers

    params.L_ofdm_syms = ceil(params.L_sym/params.N_subcarriers); % Length in ofdm symbols
    params.SNR_sweep = -20:20;
    params.n_channel_uses = 20;

    % Environment
    params.N_users = 2; 
    params.M = 50; % Number of Rx antennas (BS)
    params.phase_dist = pi; % Assumed lambda/2 antenna separation
    params.user_angles = [deg2rad(30) deg2rad(-30)];
    params.user_pwr = [1 1];
    assert(length(params.user_angles) == params.N_users, 'There must be one angle per user.')
    assert(length(params.user_pwr) == params.N_users, 'There must be one power per user.')

    % Spatial filter
    params.width = 10;

    % Channel 
    params.K = 10;
    params.N_taps = 32;
    
    % Coherent MIMO
    params.perfect_channel_estimation = false;
    params.combining_method = 'MMSE'; % Choose between 'MRC', 'ZF' or 'MMSE'

    % Non-coherent MIMO
    params.diff_decoding_dimension = 'freq';

end






