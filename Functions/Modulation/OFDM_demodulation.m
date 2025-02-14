function [rx_syms] = OFDM_demodulation(rx_ofdm_signal)
%OFDM_MODULATION removes OFDM modulation from the input
%OFDM signal matrix of the shape [N_ofdm_syms x N_subcarriers x N_users].
N_ofdm_syms = size(rx_ofdm_signal, 1);
N_ant = size(rx_ofdm_signal, 2);
N_subcarriers = size(rx_ofdm_signal, 3);
N_users = size(rx_ofdm_signal, 4);

rx_syms = 1/sqrt(N_subcarriers) * fft(rx_ofdm_signal, N_subcarriers, 3); % OFDM demodulation (maintanin Parsevals relation)
rx_syms = reshape(rx_syms, N_ofdm_syms * N_subcarriers, N_users); % Paralel to Serial (1 less symbol due to the differential demodulation)

end