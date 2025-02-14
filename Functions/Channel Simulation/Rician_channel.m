function [H] = rician_channel(user_angles, N_subcarriers, M, N_taps, K, phase_dist_ant)
%RICIAN_CHANNEL Summary of this function goes here
%   Detailed explanation goes here

N_users = length(user_angles);

    % Rayleigh part
        % With amplitude changes
        H_NLOS= 1/(sqrt(2)) .* (randn(M, N_taps, N_users) + 1j*randn(M, N_taps, N_users));
        H_NLOS_freq = fft(H_NLOS, N_subcarriers, 2);
        H_NLOS_freq = H_NLOS_freq./abs(H_NLOS_freq);
        H_NLOS = ifft(H_NLOS_freq, N_taps, 2);

        % Without amplitude changes (unitary channel)
        H_NLOS_un= 1/(sqrt(2)) .* (randn(M, N_subcarriers, N_users) + 1j*randn(M, N_subcarriers, N_users));
        H_NLOS_freq_un = fft(H_NLOS, N_subcarriers, 2);
        H_NLOS_freq_un = H_NLOS_freq_un./abs(H_NLOS_freq_un);
        H_NLOS_un = ifft(H_NLOS_freq, N_subcarriers, 2);

    % Rician part (determininstic)
        init_phases = 2*pi*rand(1,N_users); 
        rx_phases = repmat([0:M-1]', 1, N_users) * phase_dist_ant .* repmat(sin(user_angles), M, 1); % Assumed lambda/2 antenna spacing
        free_space_loss = ones(size(user_angles)); % Same free space loss for all receiving antennas for the same user (No loss channel)
        H_LOS = repmat(free_space_loss, M, 1) .* exp(-j * (repmat(init_phases, M, 1) + rx_phases)); 
        H_LOS_un = repmat(free_space_loss, M, 1) .* exp(-j * (repmat(init_phases, M, 1) + rx_phases)); 

H_LOS = cat(2, reshape(H_LOS, M, 1, N_users), zeros(M, N_taps-1, N_users)); % Extend H to the length of the OFDM signal (constant channel due to T_ofdm_sym < T_coherence)
H_LOS_un = cat(2, reshape(H_LOS_un, M, 1, N_users), zeros(M, N_subcarriers-1, N_users)); % Extend H to the length of the OFDM signal (constant channel due to T_ofdm_sym < T_coherence)


% H = sqrt(K/(K+1)) * H_LOS + sqrt(1/(K+1)) * H_NLOS; 
H = sqrt(K/(K+1)) * H_LOS_un + sqrt(1/(K+1)) * H_NLOS_un; 

end

