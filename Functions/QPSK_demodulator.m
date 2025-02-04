function [bits] = QPSK_demodulator(syms)
%QPSK_DEMODULATOR receives a matrix with the detected QPSK symbols     
% [N_syms X N_users] and maps the symbols to their corresponding bits. It
% uses the fact that the  phases of QPSK symbols are [0, pi/2, pi, 3pi/2]
bps = 2; % bits per symbol
demod_sym_idx = wrapTo2Pi(angle(syms))/(pi/2) + 1;
bits = int2bit((demod_sym_idx - 1), bps, false);

end

