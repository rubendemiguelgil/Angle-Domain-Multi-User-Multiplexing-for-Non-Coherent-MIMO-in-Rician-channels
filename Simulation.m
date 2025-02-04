clear, clc, close all;

%% Parameters
N_users = 2; 
M = 100; % Number of Rx antennas (BS)
L = 100; % Tx length (in bits)

%% Create bit strings
bits = round(rand(L, N_users));

%% Constellation modulation (DQPSK)
bps = 2; % bits per symbol
constellation = [1 1j -1 -1j];
init_diff_symbol = ones(1, N_users);

bits_symbol = reshape(bits, [bps, L/bps, N_users]);
sym_idx = squeeze(bit2int(bits_symbol, 2, false)) + 1; % +1 due to Matlab's indexing
syms = constellation(sym_idx);

diff_syms = [init_diff_symbol;  zeros(L/bps, N_users)]; %% Arreglar mod diferencial
for i = 2:L/bps + 1
    diff_syms(i, :) = diff_syms(i-1, :) .* syms(i-1, :);
end
tx_syms = diff_syms(2:end, :);
%% Constellation demodulation (QPSK) (PROBAR GRAY CODING)
constellation = [1 1j -1 -1j];

dists = abs(repelem(tx_syms, 1, 2^bps) - repmat(constellation, L/bps, N_users));
[~, det_sym_idx] = min(reshape(dists, L/bps, 2^bps, N_users), [], 2); 
det_sym_idx = squeeze(det_sym_idx);

demod_diff_sym = [init_diff_symbol; constellation(det_sym_idx)];
demod_sym = demod_diff_sym(2:end, :)./demod_diff_sym(1:end-1, :);

% Map symbols to int/bits
demod_sym_idx = wrapTo2Pi(angle(demod_sym))/(pi/2) + 1;
demod_bits = int2bit((demod_sym_idx - 1), bps, false);



%% Metrics (BER, SER, SINR)
    % SER
    error_sym = ((demod_sym_idx) - sym_idx);
    error_sym_flags = (error_sym~=0);
    SER_user = sum(error_sym_flags, 1)/L;
    SER_total = sum(error_sym_flags, 'all')/(L * N_users);
    
    % BER
    error_bits = bits - demod_bits;
    BER_user = sum(error_bits)/L;
    BER_total = sum(error_bits, 'all')/(L * N_users);
