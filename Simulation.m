clear, clc, close all;
% Determine where your m-file's folder is.
folder = fileparts(which('Simulation.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
%% Parameters
N_users = 1; 
M = 100; % Number of Rx antennas (BS)
L = 12800; % Tx length (in bits)
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

    % Rician part (determininstic)
        angles = [pi/5];% -pi/3]; % pi * (rand(1, N_users) - 0.5); % ULA (has mirror ambiguity)
        rx_phases = repmat([0:M-1]', 1, N_users) * k*d_antennas .* repmat(sin(angles), M, 1);

K = 0.5;
H = rician_channel(angles, N_subcarriers, M, N_taps, K);
%% SNR sweep loop
% SNR_sweep = [0 5 10 15 20 25 30 ];
SNR_sweep = 20;
SER_total_mtx = zeros(size(SNR_sweep));
BER_total_mtx = zeros(size(SNR_sweep));
SINR_total_mtx = zeros(size(SNR_sweep));

for SNR_idx = 1:length(SNR_sweep)

SNR_dB = SNR_sweep(SNR_idx);
N0 = (10.^(-SNR_dB/10)); % Revisar espectrograma

%% Transmission (se tienen que sumar las señales)
y = tx_ofdm_signal(ofdm_signal, H, N0);

%% Angular filtering (MRC) (perfect for now)
spatial_filter = reshape(repmat(exp(-j*rx_phases), N_subcarriers, 1, 1), M, N_subcarriers, N_users); 
spatial_filter = reshape(repelem(spatial_filter, L_ofdm_syms+1, 1, 1), L_ofdm_syms+1, M, N_subcarriers, N_users);

y_filtered_angle = 1/(M) .* fft(spatial_filter, M, 2) .* fft(y, M, 2);
y_filtered = ifft(y_filtered_angle, M, 2);


%% DFT peak selection
figure(2)
clf;
subplot(2,1, 1)
hold on
plot(abs(fft(squeeze(H(:, 1, :))))/max(abs(fft(squeeze(H(:, 1, :))))), 'DisplayName','Channel')
plot(abs(fft(sum(y(2, :, 1), 4)))/max(abs(fft(sum(y(2, :, 1), 4)))), 'DisplayName','Signal')
plot(abs(fft(squeeze(spatial_filter(1,:,1)), M, 2)'/M), 'DisplayName','SP filter')
legend()

subplot(2,1, 2)
hold on
plot(abs(fft(sum(y(2, :, 1), 4)))/max(abs(fft(sum(y(2, :, 1), 4)))), 'DisplayName','OG Signal')
plot(abs(fft(sum(y_filtered(2, :, 1), 4)))/max(abs(fft(sum(y_filtered(2, :, 1), 4)))), 'DisplayName','Filt Signal')
legend()

%% Differential OFDM decoding/demodulation/carrier dealocation
rx_syms = OFDM_diff_demodulation(y_filtered); 
rx_syms = rx_syms(1:L_sym, :);% Neglect zero padded symbols due to fixed N_subcarriers
rx_syms = rx_syms./median(abs(rx_syms),1); % AGC (set to 1) 
det_syms = QPSK_detector(rx_syms(1:L_sym, :)); % Min distance QPSK detection 
det_bits = QPSK_demodulator(det_syms); % Map symbols to bits

%% Constellation plot 
figure(1)
clf;
subplot(1, 2, 1)
    hold on, grid on
    title('Constellation User 1')
    plot(squeeze(syms(:, 1)), 'r+', 'MarkerSize', 4, 'LineWidth', 2)
    plot(squeeze(rx_syms(:, 1)), 'b*', 'MarkerSize', 4, 'LineWidth', 2)
    set(gca, 'Children', flipud(get(gca, 'Children')))
    axis('equal')

subplot(1, 2, 2)
    hold on, grid on
    title('Constellation User 2')
    plot(squeeze(syms(:, end)), 'r+', 'MarkerSize', 4, 'LineWidth', 2)
    plot(squeeze(rx_syms(:, end)), 'b*', 'MarkerSize', 4, 'LineWidth', 2)
    set(gca, 'Children', flipud(get(gca, 'Children')) )
    axis('equal')

%% Metrics (BER, SER, SINR)
    % SER
    error_sym = (det_syms - syms);
    error_sym_flags = (error_sym~=0);
    SER_user = sum(error_sym_flags, 1)/L;
    SER_total = sum(error_sym_flags, 'all')/(L * N_users)
    
    % BER
    error_bits = abs(bits - det_bits);
    BER_user = abs(error_bits)/L;
    BER_total = sum(error_bits, 'all')/(L * N_users * bps)
    
    % SINR (from EVM) 
    evm = sqrt(sum(abs(rx_syms - syms).^2, 'all')/(L_sym*N_users));
    SINR_dB = 10*log10(evm)

SER_total_mtx(SNR_idx) = SER_total;
BER_total_mtx(SNR_idx) = BER_total;
SINR_total_mtx(SNR_idx) = SINR_dB;
end
% 
% figure(3)
% % subplot(3 ,1 ,1)
%     grid on
%     title('BER')
%     plot(SNR_sweep, BER_total_mtx)
%     yscale log

% subplot(3 ,1 ,2)
%     grid on
%     title('SER')
%     plot(SNR_sweep, SER_total_mtx)
%     yscale log
% 
% subplot(3 ,1 ,3)
%     grid on
%     title('SINR (10*log10(EVM)')
%     plot(SNR_sweep, SINR_total_mtx)

%  figure(1)
%  clf
% hold on
% plot(squeeze(abs(ofdm_signal(2, : ,1))), 'DisplayName','Tx sig')
% plot(squeeze(abs(y(2, 1 ,:))), 'DisplayName','Rx sig')
% legend()
