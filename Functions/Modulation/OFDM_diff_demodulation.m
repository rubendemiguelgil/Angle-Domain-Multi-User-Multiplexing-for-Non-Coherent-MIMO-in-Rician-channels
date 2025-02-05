function [det_syms] = OFDM_diff_demodulation(ofdm_signal)
%OFDM_DIFF_MODULATION removes differential QPSK OFDM modulation from the input
%OFDM signal matrix of the shape [N_ofdm_syms x N_subcarriers x N_users]. The differential encoding
%is done between consecutive symbols in the time-frequency grid, i.e.,
%between symbols separated N_subcarriers.
N_ofdm_syms = size(ofdm_signal, 1);
N_subcarriers = size(ofdm_signal, 2);
N_users = size(ofdm_signal, 3);

ofdm_diff_syms = fft(ofdm_signal, N_subcarriers, 2); % OFDM demodulation

det_diff_ofdm_syms = QPSK_detector(reshape(ofdm_diff_syms, N_ofdm_syms*N_subcarriers, N_users)); % Min distance QPSK detection
det_diff_ofdm_syms = reshape(det_diff_ofdm_syms, N_ofdm_syms, N_subcarriers, N_users);

det_ofdm_syms = det_diff_ofdm_syms(2:end, :, :)./det_diff_ofdm_syms(1:end-1, :, :); % Differential demodulation

det_syms = reshape(det_ofdm_syms, (N_ofdm_syms - 1) * N_subcarriers, N_users); % 1 less symbol due to the differential demodulation


end

