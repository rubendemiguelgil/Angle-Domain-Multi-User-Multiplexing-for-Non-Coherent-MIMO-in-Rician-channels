function [results] = simulate_coherent(params)
%SIMULATE_COHERENT Summary of this function goes here
%   Detailed explanation goes here
%% Parameters
N_users = params.N_users; 
M = params.M; % Number of Rx antennas (BS)
L = params.L; % Tx length (in bits)
bps = params.bps; % 2 bits/symbol in QPSK
L_sym = params.L_sym; % Tx length in syms
N_subcarriers = params.N_subcarriers; % Number of dft points
% CP_length = params.CP_length; % For now didnt include it due to the narrowband assumption and the use of BER and SINR as KPIs
L_ofdm_syms = params.L_ofdm_syms; % Length in ofdm symbols
pwr = params.user_pwr;

phase_dist = params.phase_dist;
K = params.K;
N_taps = params.N_taps;
angles = params.user_angles;
%% Create bit strings
bits = round(rand(L, N_users));

%% Constellation modulation (QPSK)
syms = QPSK_modulation(bits); 
syms = syms .* repmat(pwr, L_sym, 1);

%% OFDM encoding/modulation/carrier alocation
ofdm_signal = OFDM_modulation(syms, N_subcarriers);

%% Rician channel (for now constant)
H = rician_channel(angles, N_subcarriers, M, N_taps, K, phase_dist);

%% SNR sweep loop
SNR_sweep = params.SNR_sweep;
results.SER_total_mtx = zeros(size(SNR_sweep));
results.BER_total_mtx = zeros(size(SNR_sweep));
results.SINR_total_mtx = zeros(size(SNR_sweep));

for SNR_idx = 1:length(SNR_sweep)

SNR_dB = SNR_sweep(SNR_idx);
N0 = (10.^(-SNR_dB/10)); % Revisar espectrograma

%% Transmission (se tienen que sumar las señales)
y = tx_ofdm_signal(ofdm_signal, H, N0); 

%% OFDM decoding/demodulation/carrier dealocation
y_fft = OFDM_demodulation(y);
% y_fft = 1/sqrt(N_subcarriers) * fft(y, N_subcarriers, 3); % OFDM demodulation (maintanin Parsevals relation)

%% MR combining (Ch estimation missing)
H_hat = H; 
CH_e_N0 = N0;
if params.perfect_channel_estimation == true
    CH_est_errors = zeros(size(H_hat));
else
    CH_est_errors = sqrt(CH_e_N0/2) * (randn(size(H_hat)) + j*randn(size(H_hat)));
end

H_hat_fft = fft(H_hat, N_subcarriers, 2)+ CH_est_errors; 
W = conj(H_hat_fft); 
W_exp = permute(repmat(W, 1, 1, 1, L_ofdm_syms), [4, 1, 2, 3]);

y_filtered = squeeze(sum(W_exp .* repmat(y_fft, 1, 1, 1, N_users), 2));

%% ZF combining

%% MMSE combining

%% QPSK demodulation
rx_syms = reshape(permute(y_filtered, [2, 1, 3]), N_subcarriers * L_ofdm_syms, N_users);
rx_syms = rx_syms(1:L_sym, :);% Neglect zero padded symbols due to fixed N_subcarriers
rx_syms_nm = rx_syms./mean(abs(rx_syms),1); % AGC (set to 1) 
% rx_syms_nm = rx_syms;
det_syms = QPSK_detector(rx_syms_nm); % Min distance QPSK detection 
det_bits = QPSK_demodulator(det_syms); % Map symbols to bits

%% Remove pilots from RX symbols (when calculating throughput)

%% Metrics (BER, SER, SINR, Throughput?)
    % SER
    error_sym = (det_syms - syms);
    error_sym_flags = (error_sym~=0);
    SER_user = sum(error_sym_flags, 1)/L;
    SER_total = sum(error_sym_flags, 'all')/(L * N_users)
    
    % BER
    error_bits = abs(bits - det_bits);
    BER_user = sum(abs(error_bits), 1)/L;
    BER_total = sum(error_bits, 'all')/(L * N_users * bps)
    
    % SINR (from EVM) 
    evm = sqrt(sum(abs(rx_syms_nm - syms).^2, 'all')/(L_sym*N_users));
    SINR_dB = 10*log10(evm);
   

results.SER_total_mtx(SNR_idx) = SER_total;
results.BER_total_mtx(SNR_idx) = BER_total;
results.SINR_total_mtx(SNR_idx) = SINR_dB;
end



end

