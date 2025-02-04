clear, clc, close all;
addpath('Functions/')
%% Parameters
N_users = 2; 
M = 100; % Number of Rx antennas (BS)
L = 100; % Tx length (in bits)

%% Create bit strings
bits = round(rand(L, N_users));

%% Constellation modulation (DQPSK)
syms = QPSK_modulation(bits);
tx_diff_syms = DQPSK_modulation(syms);

%% Constellation demodulation (QPSK) (PROBAR GRAY CODING)
det_diff_syms = QPSK_detector(tx_diff_syms); % Min distance detector
det_syms = DQPSK_demodulation(det_diff_syms); % Remove differential encoding
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
