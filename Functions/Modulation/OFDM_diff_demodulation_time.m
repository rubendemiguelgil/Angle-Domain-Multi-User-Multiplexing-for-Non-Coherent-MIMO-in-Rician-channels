function [rx_syms] = OFDM_diff_demodulation_time(rx_ofdm_signal)
%OFDM_DIFF_MODULATION removes differential QPSK OFDM modulation from the input
%OFDM signal matrix of the shape [N_ofdm_syms x N_subcarriers x N_users]. The differential encoding
%is done between consecutive symbols in the time-frequency grid, i.e.,
%between symbols separated N_subcarriers.
N_ofdm_syms = size(rx_ofdm_signal, 1);
N_ant = size(rx_ofdm_signal, 2);
N_subcarriers = size(rx_ofdm_signal, 3);
N_users = size(rx_ofdm_signal, 4);


rx_diff_syms = 1/sqrt(N_subcarriers) * fft(rx_ofdm_signal, N_subcarriers, 3); % OFDM demodulation (maintanin Parsevals relation)
rx_syms = rx_diff_syms(2:end, :, :, :).*conj(rx_diff_syms(1:end-1, :, :, :)); % Differential demodulation
sum_rx_syms = 1/N_ant * sum(rx_syms, 2);
% sum_rx_syms = rx_syms(:, 1, :, :);
rx_syms = reshape(permute(squeeze(sum_rx_syms), [2, 1, 3]), N_subcarriers * (N_ofdm_syms - 1), N_users);

end