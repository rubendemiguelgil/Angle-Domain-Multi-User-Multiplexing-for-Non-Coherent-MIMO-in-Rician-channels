function [tx_diff_syms] = DQPSK_modulation(syms)
%DQPSK_MODULATION takes a bidimiensional array with each user's PSK symbols as columns [N_sym x N_users]
% and applies differential modulation to it.
N_syms = size(syms, 1);
N_users = size(syms, 2);

assert(N_syms > 2, "The number of symbols must bigger than 2 for the differential modulation.");

init_diff_symbol = ones(1, N_users);
diff_syms = [init_diff_symbol;  zeros(N_syms, N_users)]; 

for i = 2:N_syms + 1
    diff_syms(i, :) = diff_syms(i-1, :) .* syms(i-1, :);
end

tx_diff_syms = diff_syms(2:end, :);

end

