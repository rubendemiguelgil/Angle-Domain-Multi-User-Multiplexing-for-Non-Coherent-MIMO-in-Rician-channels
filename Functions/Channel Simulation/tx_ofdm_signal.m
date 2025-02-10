function [rx_signal] = tx_ofdm_signal(ofdm_signal, H, N0)
%TX_OFDM_SIGNAL Summary of this function goes here
%   Detailed explanation goes here

L_ofdm_syms = size(ofdm_signal, 1);
N_subcarriers = size(ofdm_signal, 2);
N_ant = size(H, 1);

rx_signal_freq = zeros(L_ofdm_syms, N_ant, N_subcarriers);
noise_freq=sqrt(N0*1/(2)).*(randn(L_ofdm_syms, N_ant, N_subcarriers)+j*randn(L_ofdm_syms, N_ant, N_subcarriers)); % total noise variance is N0 per path

for sym_idx = 1:L_ofdm_syms % +1 due to the extra differential symbol
        rx_signal_freq(sym_idx, :, :) = sum(fft(H, N_subcarriers, 2) .* fft(ofdm_signal(sym_idx, :, :), N_subcarriers, 2), 3); % + squeeze(noise_freq(sym_idx, :, :));
end

noise=sqrt(N0*1/(2)).*(randn(L_ofdm_syms, N_ant, N_subcarriers)+j*randn(L_ofdm_syms, N_ant, N_subcarriers)); % total noise variance is N0 per path

rx_signal = ifft(rx_signal_freq, N_subcarriers, 3) + noise;

end

