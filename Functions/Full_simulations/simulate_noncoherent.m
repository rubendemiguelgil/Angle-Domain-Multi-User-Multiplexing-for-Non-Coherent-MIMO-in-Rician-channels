function [results] = simulate_noncoherent(params)
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

%% Differential OFDM encoding/modulation/carrier alocation
[ofdm_signal] = OFDM_diff_modulation_freq(syms, N_subcarriers);

%% Rician channel (for now constant)
H = rician_channel(angles, N_subcarriers, M, N_taps, K, phase_dist);

rx_phases = repmat([0:M-1]', 1, N_users) * phase_dist .* repmat(sin(angles), M, 1);
[~, user_id]= max(fft(exp(-j*rx_phases), M, 1), [], 1); % Map users for user identification 


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

%% Angular filtering (MRC?) 
[spatial_filter_time, user_mapping] = dft_peaks(y, N_users);
[spatial_filter_time] = user_identification(spatial_filter_time, user_id, user_mapping);

y_filtered_angle =  fft(spatial_filter_time, M, 2) .* fft(y, M, 2);
y_filtered = ifft(y_filtered_angle, M, 2);

%% --comparar con añadir una exponencial antes de sumarlas -- 

%% -- añadir diferencial en el dominio de la frecuencia -- 

%% -- comparar efectos del filtrado, anchura del filtro etc --

%% -- Mirar potencia de señal receptor --
%% Differential OFDM decoding/demodulation/carrier dealocation
rx_syms = OFDM_diff_demodulation_freq(y_filtered); 
rx_syms = rx_syms(1:L_sym, :);% Neglect zero padded symbols due to fixed N_subcarriers
rx_syms_nm = rx_syms./mean(abs(rx_syms),1); % AGC (set to 1) 
det_syms = QPSK_detector(rx_syms_nm); % Min distance QPSK detection 
det_bits = QPSK_demodulator(det_syms); % Map symbols to bits

%% Metrics (BER, SER, SINR)
    % SER
    error_sym = (det_syms - syms);
    error_sym_flags = (error_sym~=0);
    SER_total = sum(error_sym_flags, 'all')/(L * N_users);
    
    % BER
    error_bits = abs(bits - det_bits);
    BER_total = sum(error_bits, 'all')/(L * N_users * bps);
    
    % SINR (from EVM) 
    evm = sqrt(sum(abs(rx_syms_nm - syms).^2, 'all')/(L_sym*N_users));
    SINR_dB = 10*log10(evm);

results.SER_total_mtx(SNR_idx) = SER_total;
results.BER_total_mtx(SNR_idx) = BER_total;
results.SINR_total_mtx(SNR_idx) = SINR_dB;
end

end

