function [ofdm_signal] = OFDM_modulation(syms, N_subcarriers)
%OFDM_MODULATION applies OFDM modulation to the input
%symbol matrix of the shape [N_syms x N_users]. 
N_syms = size(syms,1);
N_users = size(syms, 2);
N_ofdm_syms = ceil(size(syms, 1)/N_subcarriers);

syms_zero_padded = zeros(N_ofdm_syms*N_subcarriers,N_users);
syms_zero_padded(1:N_syms, :) = syms; 
ofdm_syms = reshape(syms_zero_padded, [], N_subcarriers, N_users);

ofdm_signal = sqrt(N_subcarriers) * ifft(ofdm_syms, N_subcarriers, 2); % OFDM modullation (sqrt(N_subcarriers) maintains Parsevals theorem)

end

