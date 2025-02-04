function [det_syms] = QPSK_detector(rx_syms)
%QPSK_DETECTOR maps the received symbols to the closest QPSK symbol
%following a minimum distance criterion.
constellation = [1 1j -1 -1j];
bps = 2; % bits per symbol
N_syms = size(rx_syms, 1);
N_users = size(rx_syms, 2);

dists = abs(repelem(rx_syms, 1, 2^bps) - repmat(constellation, N_syms, N_users));
[~, det_sym_idx] = min(reshape(dists, N_syms, 2^bps, N_users), [], 2); 
det_syms = constellation(squeeze(det_sym_idx));

end

