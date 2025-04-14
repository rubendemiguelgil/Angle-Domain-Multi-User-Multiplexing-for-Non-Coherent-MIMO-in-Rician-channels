function [det_bits_final, det_joint_syms] = joint_const_detection(rx_syms)
%JOINT_CONST_DETECTION Summary of this function goes here
%   Detailed explanation goes here
constellation = [1+(0.707 + 0.707j) 1+(-0.707 + 0.707j) 1+(-0.707 - 0.707j) 1+(0.707 - 0.707j) ...
     j+(0.707 + 0.707j) j+(-0.707 + 0.707j) j+(-0.707 - 0.707j) j+(0.707 - 0.707j) ...
    -1+(0.707 + 0.707j) -1+(-0.707 + 0.707j) -1+(-0.707 - 0.707j) -1+(0.707 - 0.707j)...
    -j+(0.707 + 0.707j) -j+(-0.707 + 0.707j) -j+(-0.707 - 0.707j) -j+(0.707 - 0.707j)];

bps = 4; % bits per symbol
N_syms = size(rx_syms, 1);
N_users = 1;

dists = abs(repelem(rx_syms, 1, 2^bps) - repmat(constellation, N_syms, N_users));
[~, det_sym_idx] = min(reshape(dists, N_syms, 2^bps, N_users), [], 2); 
det_joint_syms = constellation(squeeze(det_sym_idx));

det_bits = int2bit(det_sym_idx-1, bps, false);
det_bits_final = zeros(length(det_bits)/2,1);

% Create index masks
idx1 = mod(0:length(det_bits)-1, 4) < 2; % Selects indices 1-2, 5-6, 9-10
idx2 = ~idx1; % Selects indices 3-4, 7-8

% Extract elements
det_bits_final(:,2) = det_bits(idx1);
det_bits_final(:,1) = det_bits(idx2);

end

