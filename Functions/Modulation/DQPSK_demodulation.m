function [det_syms] = DQPSK_demodulation(diff_det_syms)
%DQPSK_MODULATION receives a matrix with the detected DQPSK symbols
% [N_syms x N_users] and removes the differential encoding from the symbols.
N_users = size(diff_det_syms, 2);
init_diff_symbol = ones(1, N_users);

diff_det_syms = [init_diff_symbol; diff_det_syms];
det_syms = diff_det_syms(2:end, :)./diff_det_syms(1:end-1, :);

end

