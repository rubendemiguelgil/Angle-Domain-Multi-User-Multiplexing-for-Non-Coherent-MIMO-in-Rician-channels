clear, clc, close all;
% Determine where your m-file's folder is.
folder = fileparts(which('NC_simulation.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
%% Parameters
plotting = true;
N_users = 2; 
M = 100; % Number of Rx antennas (BS)
L = 10240; % Tx length (in bits)
bps = 2; % 2 bits/symbol in QPSK
L_sym = L/bps; % Tx length in syms
N_subcarriers = 1024; % Number of dft points
CP_length = 128; % For now didnt include it due to the narrowband assumption and the use of BER and SINR as KPIs
L_ofdm_syms = ceil(L_sym/N_subcarriers); % Length in ofdm symbols
n_ch_uses = 1;

%% Create bit strings
bits = round(rand(L, N_users));

%% Constellation modulation (QPSK)
syms = QPSK_modulation(bits); 
% syms(:,2) = syms(:,2) * exp(j* pi/4); 
pwr = [1 1];
joint_syms = sum(syms, 2);
%% Differential OFDM encoding/modulation/carrier alocation
[ofdm_signal] = OFDM_diff_modulation_freq(syms, N_subcarriers, pwr);

%% Rician channel (for now constant)
% Channel Parameters
phase_dist = pi; % Assumed lambda/2 antenna separation
N_taps = 32;
angles = [deg2rad(32) deg2rad(-30)]; % pi * (rand(1, N_users) - 0.5); % ULA (has mirror ambiguity)

rx_phases = repmat([0:M-1]', 1, N_users) * phase_dist .* repmat(sin(angles), M, 1);
K = 10;

H = rician_channel(angles, N_subcarriers, M, N_taps, K, phase_dist);
% H = transparent_channel(angles, N_subcarriers, M, N_taps, K, phase_dist);

H_angle = fft(H, M, 1);
[~, user_id]= max(fft(exp(-j*rx_phases), M, 1), [], 1); % Map users for user identification 


%% SNR sweep loop
% SNR_sweep = -20:5;
SNR_sweep = 0;
SER_total_mtx = zeros(length(SNR_sweep), n_ch_uses);
BER_total_mtx = zeros(length(SNR_sweep), n_ch_uses);
SINR_total_mtx = zeros(length(SNR_sweep), n_ch_uses);


for SNR_idx = 1:length(SNR_sweep)
    for ch_use = 1:n_ch_uses
    
    SNR_dB = SNR_sweep(SNR_idx);
    N0 = (10.^(-SNR_dB/10)); % Revisar espectrograma
    
    %% Transmission (se tienen que sumar las señales)
    % ofdm_signal = ones(size(ofdm_signal));
    [y, noise] = tx_ofdm_signal(ofdm_signal, H, N0);
    
    %% Angular filtering (MRC?) 
    [spatial_filter_time, user_mapping] = dft_peaks(y, N_users, 50);
    [spatial_filter_time] = user_identification(spatial_filter_time, user_id, user_mapping);
    
    y_filtered_angle =  fft(spatial_filter_time, M, 2) .* fft(y, M, 2);
    y_filtered = ifft(y_filtered_angle, M, 2);
    H_filtered_angle = fft(squeeze(spatial_filter_time(1,:,:,:)), M, 1) .* fft(H, M, 1);
    H_filtered = ifft(H_filtered_angle, M, 1);

    noise_angle = 1/sqrt(M) * fft(noise, M, 2);
    noise_filtered_angle =  fft(spatial_filter_time, M, 2) .* noise_angle;
    noise_filtered = sqrt(M) .* ifft(noise_filtered_angle, M, 2);
    
    figure(8)
    hold on, grid on
    plot(mean(abs(noise_filtered(2, :, :)),3), 'DisplayName','Filtered')
    plot(mean(abs(noise_angle(2, :, :)),3), 'DisplayName','OG')
    legend()
    %% Differential OFDM decoding/demodulation/carrier dealocation
    rx_syms = OFDM_diff_demodulation_freq(y_filtered); 
    rx_syms = rx_syms(1:L_sym, :); % Neglect zero padded symbols due to fixed N_subcarriers
    rx_syms_nm = rx_syms./mean(abs(rx_syms),1); % AGC (set to 1) 
    % rx_syms_nm = rx_syms;
    det_syms = QPSK_detector(rx_syms_nm); % Min distance QPSK detection 
    det_bits = QPSK_demodulator(det_syms); % Map symbols to bits
    % [det_bits, det_joint_syms] = joint_const_detection(rx_syms_nm); % Joint constellation
    %% Constellation plot 
    if plotting == true
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
        
        %% Channel and spatial filter plotting
        figure(2)
        clf;
        hold on
        plot(abs(fft(squeeze(sum(H(:, 1, :), 3))))./max(abs(fft(squeeze(sum(H(:, 1, :), 3)))), [], 'all'), 'DisplayName','Rice Channel antenna domain')
        plot(abs(fft(squeeze(y(2, :, 1))))./max(abs(fft(squeeze(y(2, :, 1)))), [], 'all'), 'DisplayName','Rx signal antenna domain')
        for user = 1:N_users
            plot(abs(fft(squeeze(spatial_filter_time(3,:,1, user)), M, 2)'), 'DisplayName',['Spatial filter user ' int2str(user)], 'LineWidth', 2)
        end
        legend()
        
        %% Signal and noise plotting
        figure(3)
        clf;
        hold on
        % plot(abs(squeeze(ofdm_signal(2, :, 1))), 'DisplayName','Tx signal')
        plot(abs(squeeze(noise(2, 1, :))), 'DisplayName','Noise antenna 1')
        plot(abs(squeeze(noise_filtered(2, 1, :, 1))), 'DisplayName','Filtered Noise')
        % plot(abs(squeeze(y_filtered(2, 1, :))), 'DisplayName','Filtered signal')
        % plot(abs(squeeze(y(2, 1, :, 1))), 'DisplayName','Rx signal antenna 1')
        % plot(mean(abs(squeeze(y(2, :, :))), 1), 'DisplayName', 'Mean antenna Rx signal')
        legend()
    
    end
    
    %% Metrics (BER, SER, SINR)
        % SER
        error_sym = (det_syms - syms);
        error_sym_flags = (error_sym~=0);
        SER_user = sum(error_sym_flags, 1)/L;
        SER_total = sum(error_sym_flags, 'all')/(L/bps * N_users)
        
        % BER
        error_bits = abs(bits - det_bits);
        BER_user = sum(abs(error_bits), 1)/L;
        BER_total = sum(error_bits, 'all')/(L * N_users)
        
        % SINR (from EVM) 
        evm = sqrt(sum(abs(rx_syms_nm - syms).^2, 'all')/(L_sym*N_users));
        SINR_dB = 10*log10(evm)
    
        % % SER Joint
        % error_sym = (transpose(det_joint_syms) - joint_syms);
        % error_sym_flags = (abs(error_sym)>0.1);
        % SER_user = sum(error_sym_flags, 1)/L;
        % SER_total = sum(error_sym_flags, 'all')/(length(joint_syms))
        % 
        % % SINR Joint (from EVM) 
        % evm = sqrt(sum(abs(rx_syms_nm - joint_syms).^2, 'all')/(L_sym*N_users));
        % SINR_dB = 10*log10(evm)

    SER_total_mtx(SNR_idx, ch_use) = SER_total;
    BER_total_mtx(SNR_idx, ch_use) = BER_total;
    SINR_total_mtx(SNR_idx, ch_use) = SINR_dB;
    end
end

% SER_total_mtx = mean(SER_total_mtx, 2);
% BER_total_mtx = mean(BER_total_mtx, 2);
% SINR_total_mtx = mean(SINR_total_mtx, 2);
% 
% 
% figure(4)
% % % subplot(3 ,1 ,1)
%     grid on
%     title('BER')
%     plot(SNR_sweep, BER_total_mtx)
%     yscale log
% % 
% % % subplot(3 ,1 ,2)
% % %     grid on
% % %     title('SER')
% % %     plot(SNR_sweep, SER_total_mtx)
% % %     yscale log
% % % 
% % subplot(3 ,1 ,3)
%     grid on
%     title('SINR (10*log10(EVM)')
%     plot(SNR_sweep, SINR_total_mtx)

%  figure(1)
%  clf
% hold on
% plot(squeeze(abs(ofdm_signal(2, : ,1))), 'DisplayName','Tx sig')
% plot(squeeze(abs(y(2, 1 ,:))), 'DisplayName','Rx sig')
% legend()


% %%
% clear
% 
% F = dftmtx(40);
% 
% prod = F' * F 
% 

% %% Channel shenanigans
% clear
% 
% N_taps = 50;
% decay = 0.1;
% Ch_time = exp(-[1:N_taps]' * decay) .* 1/sqrt(2) .* (randn(N_taps, 1) + j * randn(N_taps, 1));
% 
% plot(abs(Ch_time))
% 
% %% Signal
% N_samples_og = 100;
% N_samples = 1024;
% % 
% % signal_og = 1/sqrt(2) .* (randn(N_samples_og, 1) + j * randn(N_samples_og, 1));
% % signal = interp1([1:N_samples_og]./N_samples_og, signal_og, [1:N_samples]./N_samples, 'spline');
% signal = exp(-j * [1:N_samples] * 2*pi/50);
% %%
% figure(1)
% subplot(2, 2, 1)
%     hold on, grid on
%     title('Signal magnitude')
%     plot(abs(signal))
% 
% subplot(2, 2, 2)
%     hold on, grid on
%     title('Signal phase')
%     plot(unwrap(angle(signal)))
% 
% subplot(2, 2, 3)
%     hold on, grid on
%     title('DFt magnitude')
%     plot(abs(fft(signal, 1024)))
% 
% subplot(2, 2, 4)
%     hold on, grid on
%     title('DFt phase')
%     plot(unwrap(angle(fft(signal, 1024))))
% 
% 
% %%
% figure(2)
% subplot(2, 2, 1)
%     hold on, grid on
%     title('Time impulse response magnitude')
%     plot(abs(Ch_time))
% 
% subplot(2, 2, 2)
%     hold on, grid on
%     title('Time impulse response phase')
%     plot(unwrap(angle(Ch_time)))
% 
% subplot(2, 2, 3)
%     hold on, grid on
%     title('DFt magnitude')
%     plot(abs(fft(Ch_time, 1024)))
% 
% subplot(2, 2, 4)
%     hold on, grid on
%     title('DFt phase')
%     plot(unwrap(angle(fft(Ch_time, 1024))))
% 
% %% Convolution
% 
% y_cconv = cconv(signal, Ch_time, N_samples);
% y_conv = conv(signal, Ch_time, 'full');
% 
% y_prod = ifft(fft(signal) .* transpose(fft(Ch_time, N_samples)));
% 
% figure(3)
% hold on
% plot(real(signal), 'DisplayName','Signal')
% plot(real(y_conv), 'DisplayName','Conv')
% plot(real(y_prod), 'DisplayName','Prod')
% legend()

%%
% figure(5)
% H_fft = fft(H, N_subcarriers, 2);
% for ant = 1:20
%     subplot(2, 1, 1)
%     hold on
%     plot(abs(squeeze(H_fft(ant, :, 1))))
% 
%     subplot(2, 1, 2)
%     hold on
%     % plot((angle(squeeze(H_fft(ant, :, 1)))))
% 
% end
% 
% subplot(2, 1, 2)
% hold on
% % plot(unwrap(angle(squeeze(1/M.*sum(H_fft(:, :, 1),1)))), 'LineWidth',10)
% 
% 
% R_time = 1/M .* squeeze(H(:, 1, :))' * squeeze(H(:, 1, :));
% R_freq = 1/M .* squeeze(H_fft(:, 100, :))' * squeeze(H_fft(:, 101, :));
% 
% 
% H_f_diff = squeeze(H_fft(:, 1:end-1, :)) .* conj(squeeze(H_fft(:, 2:end, :)));
% 
% mean_rot = angle(1/M * sum(H_f_diff));
% 
% for ant = 1:1
%     subplot(2, 1, 2)
%     hold on
%     plot(unwrap(angle(H_f_diff(ant,:))))
% end
% plot(unwrap(mean_rot), 'LineWidth',2)

%%
spatial_filter_time_1 = repmat(squeeze(spatial_filter_time(1,:,:,1)), [1 , 1, N_users]);
H_filtered_angle_1 = fft(spatial_filter_time_1, M, 1) .* fft(H, M, 1);
H_filtered_1 = ifft(H_filtered_angle_1, M, 1);

H_filtered_freq = fft(H_filtered_1, N_subcarriers, 2);
H_filtered_freq_ang = fft(H_filtered_freq, M, 1);
mult_ant_dom =  conj(H_filtered_freq(:,10,1)) .* H_filtered_freq(:,11, 2);
mult_ang_dom = ifft(mult_ant_dom); 

H_freq = fft(H, N_subcarriers, 2);
mult_h_freq =  conj(H_freq(:,10,1)) .* H_freq(:,11, 2);

figure(7)
clf
subplot(2,1,1)
hold on; grid on
plot(1:M, abs(mult_ant_dom), 'DisplayName','Prod');
plot(1:M, H_filtered_freq(:,10,1), 'DisplayName','Fn1')
plot(1:M, H_filtered_freq(:,11,2), 'DisplayName','Fn2')
legend()

subplot(2,1,2)
hold on; grid on
plot(1:M, real(mult_ang_dom),'DisplayName','Conv')
plot(1:M, H_filtered_freq_ang(:,10,1),'DisplayName','Fn1')
plot(1:M, H_filtered_freq_ang(:,11,2), 'DisplayName','Fn2')
legend()

mean_nch_filt = abs(mean(mult_ant_dom))
mean_nch = abs(mean(mult_h_freq))

%%

filt_noise_prod = conj(noise_filtered(1, :,  10, 1)) .* noise_filtered(1, :,  11, 1);
noise_prod = conj(noise(1, :,  10)) .* noise_filtered(1, :,11);


figure(8)
hold on; grid on
plot(1:M, filt_noise_prod, 'DisplayName','Prod filt')
plot(1:M, squeeze(noise_filtered(1, :,  10, 1)), 'DisplayName','Filt')
plot(1:M, squeeze(noise_filtered(1, :,  11, 2)))

plot(1:M, noise_prod, 'DisplayName','Prod normal')
plot(1:M, squeeze(noise(1, :,  10)), 'DisplayName','Normal')
plot(1:M, squeeze(noise(1, :,  11)))

legend()

var_filt = var(filt_noise_prod)
var_normal = var(noise_prod)

mean_filt = abs(mean(filt_noise_prod))
mean_normal = abs(mean(noise_prod))
