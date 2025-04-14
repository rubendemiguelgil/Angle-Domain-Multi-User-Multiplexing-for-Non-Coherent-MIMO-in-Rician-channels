function [spatial_filter_time, user_mapping] = dft_peaks(rx_signal, N_users, width)
%DFT_PEAKS Creates a spatial filter in the time domain for each user with shape 
% [L_ofdm_sym x N_ant x N_subcarriers x N_users]. To do it, it detects the
% peaks in the DFT in the angular domain (the antenna dimension). This
% function is not capable of identifying users, that problem is assumed
% to be solved. The variable user_mapping acts as previous knowledge about
% the AoA of the users, it is used to correctly associate the spatial filters 
% with their corresponding user.

L_ofdm_syms = size(rx_signal, 1);
N_ant = size(rx_signal, 2);
N_subcarriers = size(rx_signal, 3);
half_wdth = floor(width/2);

angle_domain_signal = fft(rx_signal, N_ant, 2);
mean_angle_domain_signal = mean(angle_domain_signal, 3); % Averages the noise out

spatial_filter_angle = zeros(L_ofdm_syms, N_ant, N_users);
user_mapping = zeros(L_ofdm_syms, N_users);
for sym_idx = 1:L_ofdm_syms
    [~,locs] = findpeaks(abs(mean_angle_domain_signal(sym_idx, :)), 'SortStr','descend');
    for user = 1:N_users
        indexes = locs(user)-half_wdth:locs(user)+half_wdth; 
        indexes_mod = mod(indexes - 1 + N_ant, N_ant) + 1; %  to allow the filter to go to the borders of the array
        spatial_filter_angle(sym_idx, indexes_mod, user) = 1;
        user_mapping(sym_idx, user) = locs(user);
    end
end

spatial_filter_time = ifft(spatial_filter_angle, N_ant, 2);
spatial_filter_time = reshape(repmat(spatial_filter_time, 1, N_subcarriers, 1, 1), L_ofdm_syms, N_ant, N_subcarriers, N_users);

end

