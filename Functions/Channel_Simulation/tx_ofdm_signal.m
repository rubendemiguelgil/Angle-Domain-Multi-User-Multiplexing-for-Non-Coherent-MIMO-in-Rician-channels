function [rx_signal, noise, rx_signal_freq] = tx_ofdm_signal(ofdm_signal, H_freq, N0)
%TX_OFDM_SIGNAL Summary of this function goes here
%   Detailed explanation goes here

L_ofdm_syms = size(ofdm_signal, 1);
N_subcarriers = size(ofdm_signal, 2);
N_ant = size(H_freq, 1);

ofdm_signal_freq = 1/sqrt(N_subcarriers) .* fft(ofdm_signal, N_subcarriers, 2);

rx_signal_freq = zeros(L_ofdm_syms, N_ant, N_subcarriers);
for sym_idx = 1:L_ofdm_syms 
    rx_signal_freq(sym_idx, :, :) = sum(H_freq .*  ofdm_signal_freq(sym_idx, : , :), 3); % Convolution property allows to just extend ofdm_signal with the antenna trace of H_freq
end

noise = sqrt(N0/(2)).*(randn(L_ofdm_syms, N_ant, N_subcarriers)+j*randn(L_ofdm_syms, N_ant, N_subcarriers)); % total noise variance is N0 per path
noise_freq = sqrt(N0/(2)).*(randn(L_ofdm_syms, N_ant, N_subcarriers)+j*randn(L_ofdm_syms, N_ant, N_subcarriers)); % total noise variance is N0 per path

rx_signal = sqrt(N_subcarriers) * ifft(rx_signal_freq, N_subcarriers, 3) + noise;

rx_signal_freq = rx_signal_freq + noise_freq;

end

