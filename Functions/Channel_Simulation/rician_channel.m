function [H_freq] = rician_channel(user_angles, N_subcarriers, M, N_taps, K, phase_dist_ant)
%RICIAN_CHANNEL Summary of this function goes here
%   Detailed explanation goes here

N_users = length(user_angles);

    % Rayleigh part (random)
        % Exponentially decaying random PDPs
        decay = 1e-1;
        decrease = exp(-[0:N_taps-1]* decay);
        decrease_mat = repmat(decrease, M, 1, N_users);
        H_NLOS_pdp= 1/(sqrt(2)) .* decrease_mat .* (randn(M, N_taps, N_users) + 1j*randn(M, N_taps, N_users));
        H_NLOS_freq = fft(H_NLOS_pdp, N_subcarriers, 2);
        H_NLOS_freq_nm = H_NLOS_freq./sqrt(repmat(var(H_NLOS_freq, 1, 1), M, 1, 1));
        H_NLOS = ifft(H_NLOS_freq, N_subcarriers, 2);
        

        H_NLOS_rayleigh =  1/(sqrt(2)) .* (randn(M, N_taps, N_users) + 1j*randn(M, N_taps, N_users));
        H_NLOS_rayleigh_freq =  1/(sqrt(N_subcarriers)) .* fft(H_NLOS_rayleigh, N_subcarriers, 2);
        H_NLOS_rayleigh_freq_nm = H_NLOS_rayleigh_freq./sqrt(repmat(var(H_NLOS_rayleigh_freq, 1, 1), M, 1, 1));
        H_NLOS_rayleigh_long =  sqrt(N_subcarriers) .* ifft(H_NLOS_rayleigh_freq, N_subcarriers, 2);
    % Rician part (determininstic)
        init_phases = 2*pi*rand(1,N_users); 
        rx_phases = repmat([0:M-1]', 1, N_users) * phase_dist_ant .* repmat(sin(user_angles), M, 1); % Assumed lambda/2 antenna spacing
        free_space_loss = ones(size(user_angles)); % Same free space loss for all receiving antennas for the same user (No loss channel)
        H_LOS = repmat(free_space_loss, M, 1) .* exp(-j * (repmat(init_phases, M, 1) + rx_phases)); 

H_LOS = cat(2, reshape(H_LOS, M, 1, N_users), zeros(M, N_subcarriers-1, N_users)); % Extend H to the length of the OFDM signal (constant channel due to T_ofdm_sym < T_coherence)
H_LOS_freq = fft(H_LOS, N_subcarriers, 2);

H = sqrt(K/(K+1)) * H_LOS + sqrt(1/(K+1)) * H_NLOS; 
H_freq = sqrt(K/(K+1)) * H_LOS_freq + sqrt(1/(K+1)) * H_NLOS_freq_nm; 

end

