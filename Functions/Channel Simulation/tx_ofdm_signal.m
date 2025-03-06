function [rx_signal, noise] = tx_ofdm_signal(ofdm_signal, H, N0)
%TX_OFDM_SIGNAL Summary of this function goes here
%   Detailed explanation goes here

L_ofdm_syms = size(ofdm_signal, 1);
N_subcarriers = size(ofdm_signal, 2);
N_ant = size(H, 1);

ofdm_signal_freq = fft(ofdm_signal, N_subcarriers, 2);
H_freq = fft(H, N_subcarriers, 2);
rx_signal_freq = zeros(L_ofdm_syms, N_ant, N_subcarriers);
for sym_idx = 1:L_ofdm_syms % +1 due to the extra differential symbol
    % rx_signal_freq(sym_idx, :, :) = sum(fft(H, N_subcarriers, 2) .* 1/sqrt(N_subcarriers) .* fft(ofdm_signal(sym_idx, :, :), N_subcarriers, 2), 3); 
    rx_signal_freq(sym_idx, :, :) = sum(H_freq .* 1/sqrt(N_subcarriers) .* ofdm_signal_freq(sym_idx, : , :), 3); 
end

noise = sqrt(N0/(2)).*(randn(L_ofdm_syms, N_ant, N_subcarriers)+j*randn(L_ofdm_syms, N_ant, N_subcarriers)); % total noise variance is N0 per path

rx_signal = sqrt(N_subcarriers) * ifft(rx_signal_freq, N_subcarriers, 3);% + noise;

end

