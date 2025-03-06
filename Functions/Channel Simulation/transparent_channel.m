function [H] = transparent_channel(user_angles, N_subcarriers, M, N_taps, K, phase_dist_ant)
%RICIAN_CHANNEL Summary of this function goes here
%   Detailed explanation goes here

N_users = length(user_angles);

H = ifft(ones(M, N_subcarriers, N_users), N_subcarriers, 2);
% H = 1/M * ones(M, N_subcarriers, N_users);
end

