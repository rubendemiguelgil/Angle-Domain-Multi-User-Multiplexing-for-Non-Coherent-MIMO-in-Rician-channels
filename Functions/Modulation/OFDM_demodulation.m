function [rx_syms] = OFDM_demodulation(rx_ofdm_signal)
%OFDM_MODULATION removes OFDM modulation from each antenna of the input
%OFDM signal matrix of the shape [N_ofdm_syms x N_subcarriers x N_ant x N_users].

N_subcarriers = size(rx_ofdm_signal, 3);


rx_syms = 1/sqrt(N_subcarriers) * fft(rx_ofdm_signal, N_subcarriers, 3); % OFDM demodulation (maintanin Parsevals relation)

end