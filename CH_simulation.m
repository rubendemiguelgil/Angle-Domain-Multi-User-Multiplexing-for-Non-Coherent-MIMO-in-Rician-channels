clear, clc, close all;
% Determine where your m-file's folder is.
folder = fileparts(which('NC_simulation.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
%% Parameters
N_users = 4; 
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

%% OFDM encoding/modulation/carrier alocation
ofdm_signal = ifft(syms, N_subcarriers, 1);

%% Rician channel (for now constant)
% Channel Parameters
phase_dist = pi; % Assumed lambda/2 antenna separation
N_taps = 8;
angles = pi * (rand(1, N_users) - 0.5); % ULA (has mirror ambiguity)
% angles = [-1.0487   -0.1688    1.3234   -1.1800]; % Caso de filtro espacial demasiado estrecho
% angles =[-0.8084    1.3928   -1.0228    0.9676]; % Caso para enseñar el uso del no coherente
rx_phases = repmat([0:M-1]', 1, N_users) * phase_dist .* repmat(sin(angles), M, 1);
K = 10;

H = rician_channel(angles, N_subcarriers, M, N_taps, K, phase_dist);

H_angle = fft(H, M, 1);
[~, user_id]= max(fft(exp(-j*rx_phases), M, 1), [], 1); % Map users for user identification 


%% SNR sweep loop
% SNR_sweep = [0 5 10 15 20 25 30 ];
SNR_sweep = 20;
SER_total_mtx = zeros(size(SNR_sweep));
BER_total_mtx = zeros(size(SNR_sweep));
SINR_total_mtx = zeros(size(SNR_sweep));

% for SNR_idx = 1:length(SNR_sweep)

SNR_dB = SNR_sweep(1);
N0 = (10.^(-SNR_dB/10)); % Revisar espectrograma

%% Transmission (se tienen que sumar las señales)
y = tx_ofdm_signal(ofdm_signal, H, N0);

%% Ch estimation

H_hat = H; % For now perfect channel knowledge

%% MR combining

%% ZF combining

%% MMSE combining

%% OFDM decoding/demodulation/carrier dealocation
rx_syms = OFDM_diff_demodulation(y_filtered); 
rx_syms = rx_syms(1:L_sym, :);% Neglect zero padded symbols due to fixed N_subcarriers
rx_syms_nm = rx_syms./mean(abs(rx_syms),1); % AGC (set to 1) 
det_syms = QPSK_detector(rx_syms_nm); % Min distance QPSK detection 
det_bits = QPSK_demodulator(det_syms); % Map symbols to bits

%% Constellation plot 
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

%% Metrics (BER, SER, SINR)
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
    evm = sqrt(sum(abs(rx_syms - syms).^2, 'all')/(L_sym*N_users));
    SINR_dB = 10*log10(evm)

SER_total_mtx(SNR_idx) = SER_total;
BER_total_mtx(SNR_idx) = BER_total;
SINR_total_mtx(SNR_idx) = SINR_dB;
% end
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
