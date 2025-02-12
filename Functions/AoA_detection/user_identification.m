function [spatial_filter_time] = user_identification(spatial_filter_time, user_id, user_mapping)
%USER_IDENTIFICATION rearranges the spatial filter so it associates each
%signal to the corresponding user. The variable user_mapping represents
%previous knowledge about the AoA of the users.

L_ofdm_syms = size(user_mapping, 1);
N_users = size(user_mapping, 2);

id_dists = abs(reshape(repelem(user_mapping, 1, N_users), L_ofdm_syms, N_users, N_users) - repmat(user_id, L_ofdm_syms, 1, N_users));
[~, user_idx] = min(id_dists, [], 2); 
spatial_filter_time_aux = spatial_filter_time;
for sym = 1:L_ofdm_syms
    for user = 1:N_users
        spatial_filter_time(sym, :, :, user_idx(sym, user)) = spatial_filter_time_aux(sym, :, :, user);
    end
end

end

