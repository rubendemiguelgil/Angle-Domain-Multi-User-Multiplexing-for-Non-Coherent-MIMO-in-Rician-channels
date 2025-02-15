function [ofdm_signal] = OFDM_diff_modulation(syms, N_subcarriers)
%OFDM_DIFF_MODULATION applies differential QPSK OFDM modulation to the input
%symbol matrix of the shape [N_syms x N_users]. The differential encoding
%is done between consecutive symbols in the time-frequency grid, i.e.,
%between symbols separated N_subcarriers. Note that the differential
%encoding adds one extra symbol (needed to eliminate the effect of the
%channel's phase).
N_syms = size(syms,1);
N_users = size(syms, 2);
N_ofdm_syms = ceil(size(syms, 1)/N_subcarriers);

syms_zero_padded = zeros(N_ofdm_syms*N_subcarriers,N_users);
syms_zero_padded(1:N_syms, :) = syms; 
ofdm_syms = permute(reshape(syms_zero_padded(:), N_subcarriers, N_ofdm_syms, N_users), [2, 1, 3]);

init_diff_symbol = ones(1, N_subcarriers, N_users);
diff_syms = cat(1, init_diff_symbol,  zeros(N_ofdm_syms, N_subcarriers, N_users)); 

for i = 2:N_ofdm_syms + 1
    diff_syms(i, :, :) = diff_syms(i-1, :, :) .* ofdm_syms(i-1, :, :); % Differential modulation
end

tx_diff_syms = diff_syms;
ofdm_signal = sqrt(N_subcarriers) * ifft(tx_diff_syms, N_subcarriers, 2); % OFDM modullation (sqrt(N_subcarriers) maintains Parsevals theorem)

end

