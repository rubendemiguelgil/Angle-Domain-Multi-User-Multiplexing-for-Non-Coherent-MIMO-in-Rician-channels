clear, clc, close all;
% Determine where your m-file's folder is.
folder = fileparts(which('Simulation.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
%% Parameters
N_users = 2; 
M = 100; % Number of Rx antennas (BS)
L = 1280; % Tx length (in bits)
bps = 2; % 2 bits/symbol in QPSK
L_sym = L/bps; % Tx length in syms
N_subcarriers = 1024; % Number of dft points
L_ofdm_syms = ceil(L_sym/N_subcarriers); % Length in ofdm symbols

%% Create bit strings
bits = round(rand(L, N_users));

%% Constellation modulation (QPSK)
syms = QPSK_modulation(bits); 
% pwrs = repmat([1 0.5], length(syms), 1);
% syms = syms .* pwrs;
%% Differential OFDM encoding/modulation/carrier alocation
[ofdm_signal] = OFDM_diff_modulation(syms, N_subcarriers);

%% Rician channel (for now constant)
fc = 2e9;
lambda = 1/fc;
d_antennas = lambda/2;
k = 2*pi/lambda; %  Wave number
    % Rayleigh part
    N_taps = 8;
    H_NLOS= 1/(sqrt(2)) .* (randn(M, N_taps, N_users) + 1j*randn(M, N_taps, N_users));
    H_NLOS_freq = fft(H_NLOS, N_subcarriers, 2);
    H_NLOS_freq = H_NLOS_freq./abs(H_NLOS_freq);
    H_NLOS = ifft(H_NLOS_freq, N_taps, 2);
    % Rician part (determininstic)
        r_cell = 500;
        angles = [pi/4 -pi/3]; % pi * (rand(1, N_users) - 0.5); % ULA (has mirror ambiguity)
        distances = r_cell * randn(1, N_users); % Assumed normal distribution in distance
        init_phases = 2*pi*rand(1,N_users); 
        rx_phases = repmat([0:M-1]', 1, N_users) * k*d_antennas .* repmat(sin(angles), M, 1);
        free_space_loss = 1 * ones(size(distances)); % Same free space loss for all receiving antennas for the same user
        H_LOS = repmat(free_space_loss, M, 1) .* exp(-j * (repmat(init_phases, M, 1) + rx_phases)); 

H_LOS = cat(2, reshape(H_LOS, M, 1, N_users), zeros(M, N_taps-1, N_users)); % Extend H to the length of the OFDM signal (constant channel due to T_ofdm_sym < T_coherence)

K = 0.7;
H = (K) * H_LOS + (1-K) * H_NLOS; % Falta la K

%% Transmission (se tienen que sumar las señales)

y = tx_ofdm_signal(ofdm_signal, H);

%% DFT peak selection
figure(2)
hold on
plot(abs(fft(squeeze(H(:, 1, :))))/max(abs(fft(squeeze(H(:, 1, :))))), 'DisplayName','Channel')
plot(abs(fft(sum(y(2, :, 1), 4)))/max(abs(fft(sum(y(2, :, 1), 4)))), 'DisplayName','Signal')
legend()
%% Angular filtering (MRC) (perfect for now)
spatial_filter = reshape(repmat(exp(-j*rx_phases), N_subcarriers, 1, 1), M, N_subcarriers, N_users); 
spatial_filter = reshape(repelem(spatial_filter, L_ofdm_syms+1, 1, 1), L_ofdm_syms+1, M, N_subcarriers, N_users);

y_mrc_angle = 1/M * fft(spatial_filter, M, 2) .* fft(y, M, 2);
y_mrc = ifft(y_mrc_angle, M, 2);
% y_mrc_conv = zeros(size(y_mrc))
% for user = 1:N_users
%     for ofdm_sym = 
%         y_mrc_conv = conv(spatial_filter, y);
%     end
% end
y_nch = y_mrc;%  squeeze(sum(y_mrc, 2));
plot(abs(fft(spatial_filter(2, :, 1)))/max(abs(fft(spatial_filter(2, :, 1)))), 'DisplayName','Spac_filt 1')
plot(abs(fft(spatial_filter(2, :, 2)))/max(abs(fft(spatial_filter(2, :, 2)))), 'DisplayName','Spac_filt 1')
% figure(2)
% hold on
% plot(abs(squeeze(ofdm_signal(2,:,1))), 'DisplayName', 'Emitted sig')
% % plot(abs(squeeze(y_mrc(2,1,:,1))), 'DisplayName', 'Filtered sig')
% plot(abs(squeeze(y_nch(2, :, 1))), 'DisplayName', 'Nch Filtered sig')
% legend()

%% Differential OFDM decoding/demodulation/carrier dealocation
rx_syms = OFDM_diff_demodulation(y_nch); 
det_syms = QPSK_detector(rx_syms(1:L_sym, :)); % Min distance QPSK detection (Neglect zero padded symbols due to fixed N_subcarriers)
det_bits = QPSK_demodulator(det_syms); % Map symbols to bits

%% Metrics (BER, SER, SINR)
    % SER
    error_sym = (det_syms - syms);
    error_sym_flags = (error_sym~=0);
    SER_user = sum(error_sym_flags, 1)/L;
    SER_total = sum(error_sym_flags, 'all')/(L * N_users)
    
    % BER
    error_bits = bits - det_bits;
    BER_user = sum(error_bits)/L;
    BER_total = sum(error_bits, 'all')/(L * N_users)
    
    % SINR 