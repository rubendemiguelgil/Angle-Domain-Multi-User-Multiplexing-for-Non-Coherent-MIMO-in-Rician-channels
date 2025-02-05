clear, clc, close all;
% Determine where your m-file's folder is.
folder = fileparts(which('Simulation.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
%% Parameters
N_users = 2; 
M = 100; % Number of Rx antennas (BS)
L = 128; % Tx length (in bits)
bps = 2; % 2 bits/symbol in QPSK
L_sym = L/bps; % Tx length in syms
N_subcarriers = 32; % Number of dft points


%% Create bit strings
bits = round(rand(L, N_users));

%% Constellation modulation (QPSK)
syms = QPSK_modulation(bits);

%% Differential OFDM encoding/modulation/carrier alocation
% Pack symbols into N_subcarriers
[ofdm_signal] = OFDM_diff_modulation(syms, N_subcarriers);

%% Carrier modulation/demodulation


%% Differential OFDM decoding/demodulation/carrier dealocation
det_syms_ofdm = OFDM_diff_demodulation(ofdm_signal); 
det_syms = det_syms_ofdm(1:L_sym, :); % Neglect zero padded symbols due to fixed N_subcarriers
det_bits = QPSK_demodulator(det_syms); % Map symbols to bits

%% Metrics (BER, SER, SINR)
    % SER
    error_sym = (det_syms - syms);
    error_sym_flags = (error_sym~=0);
    SER_user = sum(error_sym_flags, 1)/L;
    SER_total = sum(error_sym_flags, 'all')/(L * N_users);
    
    % BER
    error_bits = bits - det_bits;
    BER_user = sum(error_bits)/L;
    BER_total = sum(error_bits, 'all')/(L * N_users);
