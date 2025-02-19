clear, clc, close all;
% Determine where your m-file's folder is.
folder = fileparts(which('CH_simulation.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
%% Parameters
plotting = false;
N_users = 2; 
M = 100; % Number of Rx antennas (BS)
L = 1024000; % Tx length (in bits)
bps = 2; % 2 bits/symbol in QPSK
L_sym = L/bps; % Tx length in syms
N_subcarriers = 1024; % Number of dft points
L_ofdm_syms = ceil(L_sym/N_subcarriers); % Length in ofdm symbols

%% Create bit strings
bits = round(rand(L, N_users));

%% Constellation modulation (QPSK)
syms = QPSK_modulation(bits); 
% pwr = [1 1.1];
% syms = syms .* repmat(pwr, L_sym, 1);
%% OFDM encoding/modulation/carrier alocation
ofdm_signal = OFDM_modulation(syms, N_subcarriers);

%% Rician channel (for now constant)
% Channel Parameters
phase_dist = pi; % Assumed lambda/2 antenna separation
N_taps = 8;
angles = [0.5 -0.5]; % pi * (rand(1, N_users) - 0.5); % ULA (has mirror ambiguity)
% angles = [1.5249   -0.5759   -1.5297   -1.3111]; % Interferencia en coherente
rx_phases = repmat([0:M-1]', 1, N_users) * phase_dist .* repmat(sin(angles), M, 1);
K = 10;

H = rician_channel(angles, N_subcarriers, M, N_taps, K, phase_dist);

H_angle = fft(H, M, 1);
[~, user_id]= max(fft(exp(-j*rx_phases), M, 1), [], 1); % Map users for user identification 


%% SNR sweep loop
SNR_sweep = -20:5;
% SNR_sweep = 15;
SER_total_mtx = zeros(size(SNR_sweep));
BER_total_mtx = zeros(size(SNR_sweep));
SINR_total_mtx = zeros(size(SNR_sweep));

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
CH_est_errors = sqrt(CH_e_N0/2) * (randn(size(H_hat)) + j*randn(size(H_hat)));

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

%% Plotting
if plotting == true
    % Constellation plot
    figure(1)
    clf;
    for user = 1:N_users
    subplot(ceil(sqrt(N_users)), ceil(sqrt(N_users)), user)
        hold on, grid on
        title(['Constellation User ' int2str(user)])
        plot(squeeze(syms(:, user)), 'r+', 'MarkerSize', 4, 'LineWidth', 2)
        plot(squeeze(rx_syms_nm(:, user)), 'b*', 'MarkerSize', 4, 'LineWidth', 2)
        set(gca, 'Children', flipud(get(gca, 'Children')))
        axis('equal')
    end
    % Channel and spatial filter plotting
    figure(2)
    clf;
    hold on
    plot(abs(fft(squeeze(sum(H(:, 1, :), 3))))./max(abs(fft(squeeze(sum(H(:, 1, :), 3)))), [], 'all'), 'DisplayName','Rice Channel antenna domain')
    plot(abs(fft(squeeze(y(1, :, 1))))./max(abs(fft(squeeze(y(1, :, 1)))), [], 'all'), 'DisplayName','Rx signal antenna domain')
    
    legend()
    
    % Signal and noise plotting
    figure(3)
    clf;
    hold on
    plot(abs(squeeze(ofdm_signal(2, :, 1))), 'DisplayName','Tx signal')
    plot(abs(squeeze(y(2, 1, :, 1))), 'DisplayName','Rx signal antenna 1')
    legend()
end

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
    SINR_dB = 10*log10(evm)
   

SER_total_mtx(SNR_idx) = SER_total;
BER_total_mtx(SNR_idx) = BER_total;
SINR_total_mtx(SNR_idx) = SINR_dB;
end

figure(4)
% subplot(3 ,1 ,1)
    % grid on
    % title('BER')
    % plot(SNR_sweep, BER_total_mtx)
    % yscale log

% subplot(3 ,1 ,2)
%     grid on
%     title('SER')
%     plot(SNR_sweep, SER_total_mtx)
%     yscale log
% 
% subplot(3 ,1 ,3)
    grid on
    title('SINR (10*log10(EVM)')
    plot(SNR_sweep, SINR_total_mtx)

%  figure(1)
%  clf
% hold on
% plot(squeeze(abs(ofdm_signal(2, : ,1))), 'DisplayName','Tx sig')
% plot(squeeze(abs(y(2, 1 ,:))), 'DisplayName','Rx sig')
% legend()

% %%
% LS = [1 1e1 1e2 1e3 1e4 1e5 1e6 1e7];
% prod = zeros(length(LS), 1);
% var_1 = zeros(length(LS), 1);
% var_2 = zeros(length(LS), 1);
% 
% for i = 1:length(LS) 
% L = LS(i);
% 
% noise_1 = sqrt(1/(2)).*(randn(L, 1)+j*randn(L, 1));
% noise_2 = sqrt(1/(2)).*(randn(L, 1)+j*randn(L, 1));
% 
% prod(i) = 1/L * abs(noise_1' * noise_2);
% var_1(i) = 1/L * abs(noise_1' * noise_1);
% var_2(i) = 1/L * abs(noise_2' * noise_2);
% end
% 
% plot(LS, prod)
% xscale('log')

