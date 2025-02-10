function [syms] = QPSK_modulation(bits)
%QPSK_MODULATION takes a bidimiensional array with each user's bits as columns [L x N_users]
% and modulates it with a  QPSK constellation. All symbols have
% magnitude one.

constellation = [1 1j -1 -1j];
bps = 2; % bits per symbol
L =  size(bits, 1); 
N_users = size(bits, 2);

assert(rem(L, 2) == 0, "The number of bits must be even and bigger than 2.");

bits_symbol = reshape(bits, [bps, L/bps, N_users]);
sym_idx = squeeze(bit2int(bits_symbol, 2, false)) + 1; % +1 due to Matlab's indexing
syms = reshape(constellation(sym_idx), L/bps, N_users);

end

