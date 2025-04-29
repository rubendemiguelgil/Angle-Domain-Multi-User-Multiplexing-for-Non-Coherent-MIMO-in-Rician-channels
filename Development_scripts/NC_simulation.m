clear, clc, close all;
% Determine where your m-file's folder is.
folder = fileparts(which('NC_simulation.m')); 
folder = [folder '\..'];
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%% Parameters
plotting = true;
N_users = 6; 
M = 50; % Number of Rx antennas (BS)
L = 10220; % Tx length (in bits)
bps = 2; % 2 bits/symbol in QPSK
L_sym = L/bps; % Tx length in syms
N_subcarriers = 1024; % Number of dft points
CP_length = 128; % For now didnt include it due to the narrowband assumption and the use of BER and SINR as KPIs
L_ofdm_syms = ceil(L_sym/N_subcarriers); % Length in ofdm symbols
n_ch_uses = 1;
W = 5;

%% Create bit strings
bits = round(rand(L, N_users));

%% Constellation modulation (QPSK)
syms = QPSK_modulation(bits); 
% syms(:,2) = syms(:,2) * exp(j* pi/4); 
pwr = [1 1 1 1 1 1];
joint_syms = sum(syms, 2);
%% Differential OFDM encoding/modulation/carrier alocation
[ofdm_signal] = OFDM_diff_modulation_freq(syms, N_subcarriers, pwr);

%% SNR sweep loop
% SNR_sweep = -20:20;
SNR_sweep = 2;
SER_total_mtx = zeros(length(SNR_sweep), n_ch_uses);
BER_total_mtx = zeros(length(SNR_sweep), n_ch_uses);
SINR_total_mtx = zeros(length(SNR_sweep), n_ch_uses);


for SNR_idx = 1:length(SNR_sweep)
    for ch_use = 1:n_ch_uses
    
    SNR_dB = SNR_sweep(SNR_idx);
    N0 = (10.^(-SNR_dB/10)); % Revisar espectrograma
    

    %% Rician channel (for now constant)
    % Channel Parameters
    phase_dist = pi; % Assumed lambda/2 antenna separation
    N_taps = 32;
    angles = [deg2rad(50) deg2rad(-50) deg2rad(30) deg2rad(-30) deg2rad(10) deg2rad(-10)]; 
    
    rx_phases = repmat([0:M-1]', 1, N_users) * phase_dist .* repmat(sin(angles), M, 1);
    K = 10;
    
    H_freq = rician_channel(angles, N_subcarriers, M, N_taps, K, phase_dist);
    
    [~, user_id]= max(fft(exp(-j*rx_phases), M, 1), [], 1); % Map users for user identification 



    %% Transmission (se tienen que sumar las señales)
    [y, noise] = tx_ofdm_signal(ofdm_signal, H_freq, N0);
    
    %% Angular filtering
    [spatial_filter_angle, user_mapping] = dft_peaks(y, N_users, W);
    [spatial_filter_angle] = user_identification(spatial_filter_angle, user_id, user_mapping);

    y_filtered_angle =  spatial_filter_angle .* 1/sqrt(M) .* fft(y, M, 2);
    y_filtered = sqrt(M) .* ifft(y_filtered_angle, M, 2);

    %%
    H_freq_angle = 1/sqrt(M) * fft(H_freq, M, 1);
    H_freq_filtered_angle = repmat(squeeze(spatial_filter_angle(1,:,:,1)), 1, 1, N_users) .* H_freq_angle;
    H_freq_filtered = sqrt(M) .* ifft(H_freq_filtered_angle, M, 1);
    % 
    % figure(9)
    % hold on
    % plot(abs(H_freq_angle(:,19,1)),'DisplayName','UnFiltered')
    % plot(abs(H_freq_filtered_angle(:,19,1)),'DisplayName','Filtered')
    % plot(squeeze(spatial_filter_angle(1,:,1,1)), 'Displayname','Spatial_filt')
    % legend()

    % figure(10)
    % hold on, grid on
    % plot(abs(H_freq_filtered(:,19,2)),'DisplayName','Filtered')
    % plot(abs(H_freq(:,19,2)),'DisplayName','UnFiltered')
    % legend()
    % 
    % var_filt = var(H_freq_filtered(:,19,2))
    % var_unfilt = var(H_freq(:,19,2))
    % 
    % mean_ch_filt = mean(conj(H_freq_filtered(:,19,2)) .* H_freq_filtered(:,20,2))
    % mean_ch_unfilt = mean(conj(H_freq(:,19,2)) .* H_freq(:,20,2))
    %%
    noise_angle = 1/sqrt(M) * fft(noise, M, 2);
    noise_filtered_angle = repmat(spatial_filter_angle(1,:,:,:), L_ofdm_syms, 1, 1, 1) .* noise_angle;
    noise_filtered = sqrt(M) .* ifft(noise_filtered_angle, M, 2);
    noise_freq = 1/sqrt(N_subcarriers) .* fft(noise, N_subcarriers, 3);
    noise_filtered_freq = 1/sqrt(N_subcarriers) .* fft(noise_filtered, N_subcarriers, 3);

    % prod_filt = conj(noise_filtered_freq(:, :, 3)) .* noise_filtered_freq(:, :, 4);
    % mean_prod_filt = mean(prod_filt, 2);
    % prod_unfilt = conj(noise_freq(:, :, 3)) .* noise_freq(:, :, 4);
    % mean_prod_unfilt = mean(prod_unfilt, 2);
    % 
    % figure(8)
    % hold on
    % plot(real(mean_prod_unfilt), 'DisplayName','Unfilt')
    % plot(real(mean_prod_filt), 'DisplayName','Filt')
    % 
    % pwr_noise = var(noise_freq, 1, 'all')
    % pwr_unfilt = var(mean_prod_unfilt)
    % pwr_filt = var(mean_prod_filt)
    % 
    % %%
    % cprod_filt = zeros(size(noise_filtered_freq(:, :, 3))); 
    % cprod_unfilt = zeros(size(noise_freq(:, :, 3))); 
    % for i = 1 : L_ofdm_syms
    %     cprod_filt(i, :) = noise_filtered_freq(i, :, 3)' .* H_freq_filtered(: , 4);
    %     cprod_unfilt(i, :) = noise_freq(i, :, 3)' .* H_freq(: , 4);
    % end
    % mean_cprod_filt = mean(cprod_filt, 2);
    % mean_cprod_unfilt = mean(cprod_unfilt, 2);
    % 
    % pwr_unfilt = var(mean_cprod_unfilt)
    % pwr_filt = var(mean_cprod_filt)
    % 
    % figure(8)
    % hold on
    % plot(real(mean_cprod_unfilt), 'DisplayName','Unfilt')
    % plot(real(mean_cprod_filt), 'DisplayName','filt')
    % legend()

    % %%
    % 
    % mean_h_filt = mean(conj(H_freq_filtered(:, 3)) .* H_freq_filtered(: , 3))
    % mean_h_unfilt = mean(conj(H_freq(:, 3)) .* H_freq(: , 3))

  

    %% Differential OFDM decoding/demodulation/carrier dealocation
    rx_syms = OFDM_diff_demodulation_freq(y_filtered); 
    rx_syms = rx_syms(1:L_sym, :); % Neglect zero padded symbols due to fixed N_subcarriers
    % rx_syms_nm = rx_syms./mean(abs(rx_syms),1); % to show scale 1 constellation plots (does not affect the SER due to the use of QPSK) 
    rx_syms_nm = rx_syms;
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
        plot(abs(squeeze(sum(H_freq_filtered(:,1, :), 3)))./max(abs(squeeze(sum(H_freq(:, 1, :), 3))), [], 'all'), 'DisplayName','Rice Channel antenna domain')
        plot(abs(fft(squeeze(y(2, :, 1))))./max(abs(fft(squeeze(y(2, :, 1)))), [], 'all'), 'DisplayName','Rx signal antenna domain')
        for user = 1:N_users
            for t = 1:1
                plot(abs(squeeze(spatial_filter_angle(t,:,10, user))), 'DisplayName',['Spatial filter user ' int2str(user)], 'LineWidth', 2)
            end
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
        evm = (sum(abs(rx_syms_nm - syms).^2, 'all')/(L_sym*N_users));
        SINR_dB = -10*log10(evm)
        sinr_test = var(rx_syms_nm - syms, 1, "all");
        sinr_test_db = -10*log10(sinr_test)
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

SER_total_mtx = mean(SER_total_mtx, 2);
BER_total_mtx = mean(BER_total_mtx, 2);
SINR_total_mtx = mean(SINR_total_mtx, 2);
%%

U = N_users;
sigma = sqrt(10.^(-SNR_sweep/10));

S = (K/(1+K) + (W/M) * 1/(K+1));
I = W/M^2 * (2* U * (U - 1)+ (M/W)* K * (U-1))/(1+K).^2 + ((W/M) * ((U-1)/(1+K)));
C = 2*(W/M^2) .* (U)/(K+1) .* sigma.^2 + 2 * (1/M) .* (K)/(K+1) .* sigma.^2;
N = (W/M^2) .* sigma.^4;

SINR_analytical = S ./ (I + C + N);
SINR_dB_analytical = 10*log10(SINR_analytical);

SINR_analytical_ana = M*U ./ ((U-1)^2 + 2 * U * sigma.^2 + sigma.^4);
SINR_dB_analytical_ana = 10*log10(SINR_analytical_ana);

figure(4)
% subplot(3 ,1 ,1)
%     grid on
%     title('BER')
%     plot(SNR_sweep, BER_total_mtx)
%     yscale log
% 
% subplot(3 ,1 ,2)
    grid on
    title('SER')
    plot(SNR_sweep, SER_total_mtx)
    yscale log
% 
% subplot(3 ,1 ,3)
figure(5)  
hold on, grid on
    title('SINR (10*log10(EVM)')
    plot(SNR_sweep, SINR_total_mtx, 'DisplayName','Simulation')
    plot(SNR_sweep, SINR_dB_analytical, 'DisplayName','Analytical')
        plot(SNR_sweep, SINR_dB_analytical_ana, 'DisplayName','Analytical_ana')
    legend()
%  figure(1)
%  clf
% hold on
% plot(squeeze(abs(ofdm_signal(2, : ,1))), 'DisplayName','Tx sig')
% plot(squeeze(abs(y(2, 1 ,:))), 'DisplayName','Rx sig')
% legend()
% 

