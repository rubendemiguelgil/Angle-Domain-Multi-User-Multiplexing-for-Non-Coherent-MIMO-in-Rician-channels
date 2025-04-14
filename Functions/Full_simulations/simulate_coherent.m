function [results] = simulate_coherent(params)
%SIMULATE_COHERENT this function simulates a coherent scheme with MRC, ZF
%or MMSE precoding
%% Parameters
N_users = params.N_users; 
M = params.M; % Number of Rx antennas (BS)
L = params.L; % Tx length (in bits)
bps = params.bps; % 2 bits/symbol in QPSK
L_sym = params.L_sym; % Tx length in syms
N_subcarriers = params.N_subcarriers; % Number of dft points
L_ofdm_syms = params.L_ofdm_syms; % Length in ofdm symbols
pwr = params.user_pwr;

phase_dist = params.phase_dist;
K = params.K;
N_taps = params.N_taps;
angles = params.user_angles;
n_ch_uses = params.n_channel_uses;

%% Create bit strings
bits = round(rand(L, N_users));

%% Constellation modulation (QPSK)
syms = QPSK_modulation(bits); 
syms = syms .* repmat(pwr, L_sym, 1);

%% OFDM encoding/modulation/carrier alocation
ofdm_signal = OFDM_modulation(syms, N_subcarriers);

%% Rician channel 
H = rician_channel(angles, N_subcarriers, M, N_taps, K, phase_dist);

%% SNR sweep loop
SNR_sweep = params.SNR_sweep;
SER_total_mtx = zeros(length(SNR_sweep), n_ch_uses);
BER_total_mtx = zeros(length(SNR_sweep), n_ch_uses);
SINR_total_mtx = zeros(length(SNR_sweep), n_ch_uses);

for SNR_idx = 1:length(SNR_sweep)
    for ch_use = 1:n_ch_uses
             
        SNR_dB = SNR_sweep(SNR_idx);
        N0 = (10.^(-SNR_dB/10));
        
        %% Transmission 
        y = tx_ofdm_signal(ofdm_signal, H, N0); 
        
        %% OFDM decoding/demodulation/carrier dealocation
        y_fft = OFDM_demodulation(y);
 
        %% Ch estimation
        H_hat = H; 
        CH_e_N0 = N0;
        if params.perfect_channel_estimation == true
            CH_est_errors = zeros(size(H_hat));
        else
            CH_est_errors = sqrt(CH_e_N0/2) * (randn(size(H_hat)) + j*randn(size(H_hat)));
        end
        
        H_hat_fft = fft(H_hat, N_subcarriers, 2)+ CH_est_errors; 
        W = conj(H_hat_fft); 
        
        %% MIMO combining
        switch params.combining_method
            case 'MRC'
                W_exp = permute(repmat(W, 1, 1, 1, L_ofdm_syms), [4, 1, 2, 3]);
                
                y_filtered = squeeze(sum(W_exp .* repmat(y_fft, 1, 1, 1, N_users), 2));
                
            case 'ZF'
                H_fft_zf = permute(H_hat_fft, [1,3,2]); % Subcarrier dim last because it stays unaltered
                H_fft_zf_herm = permute(conj(H_hat_fft), [3,1,2]);
                
                H_sqr_inv = pageinv(pagemtimes(H_fft_zf_herm, H_fft_zf));
                W_zf = pagemtimes(H_sqr_inv, H_fft_zf_herm);
                
                y_filtered = zeros(L_ofdm_syms, N_subcarriers, N_users);
                for sym = 1:L_ofdm_syms
                    y_filtered(sym, :, :) = transpose(squeeze(pagemtimes(W_zf, permute(y_fft(sym, :, :), [2, 1, 3]))));
                end
            
            case 'MMSE'
                H_fft_mmse = permute(H_hat_fft, [1,3,2]); % Subcarrier dim last because it stays unaltered
                H_fft_mmse_herm = permute(conj(H_hat_fft), [3,1,2]);
                
                C_k = pagemtimes(H_fft_mmse, H_fft_mmse_herm) + N0 * eye(M);
                W_mmse = pagemtimes(pageinv(C_k), H_fft_mmse);
                W_mmse_herm = permute(conj(W_mmse), [2,1,3]);
                
                y_filtered = zeros(L_ofdm_syms, N_subcarriers, N_users);
                for sym = 1:L_ofdm_syms
                    y_filtered(sym, :, :) = transpose(squeeze(pagemtimes(W_mmse_herm, permute(y_fft(sym, :, :), [2, 1, 3]))));
                end
        
        end
        %% QPSK demodulation
        rx_syms = reshape(permute(y_filtered, [2, 1, 3]), N_subcarriers * L_ofdm_syms, N_users);
        rx_syms = rx_syms(1:L_sym, :);% Neglect zero padded symbols due to fixed N_subcarriers
        rx_syms_nm = rx_syms./mean(abs(rx_syms),1); 
        det_syms = QPSK_detector(rx_syms_nm); % Min distance QPSK detection 
        det_bits = QPSK_demodulator(det_syms); % Map symbols to bits
        
        %% Metrics (BER, SER, SINR)
            % SER
            error_sym = (det_syms - syms);
            error_sym_flags = (error_sym~=0);
            SER_total = sum(error_sym_flags, 'all')/(L_sym * N_users);
            
            % BER
            error_bits = abs(bits - det_bits);
            BER_total = sum(error_bits, 'all')/(L * N_users);
            
            % SINR (from EVM) 
            evm = sqrt(sum(abs(rx_syms_nm - syms).^2, 'all')/(L_sym*N_users));
            SINR_dB = -10*log10(evm);
           
        
            SER_total_mtx(SNR_idx, ch_use) = SER_total;
            BER_total_mtx(SNR_idx, ch_use) = BER_total;
            SINR_total_mtx(SNR_idx, ch_use) = SINR_dB;
            
            text = strcat("Ch use number: ", int2str(ch_use) ,"; SNR :", int2str(SNR_dB));
            disp(text)
    end
end

results.SER_total_mtx = mean(SER_total_mtx, 2);
results.BER_total_mtx = mean(BER_total_mtx, 2);
results.SINR_total_mtx = mean(SINR_total_mtx, 2);

end

