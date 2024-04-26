function user_rate = get_user_rate_cvx(gamma_0, num_antenna, p_max, z_user, num_user, z_target_l, num_target, E, sensing_th, A, user_rate_ISAC_sum)

    user_rate_comm = (rel_entr(z_user / gamma_0* num_antenna* p_max, z_user / gamma_0* num_antenna* p_max + 1) + rel_entr(z_user / gamma_0* num_antenna* p_max + 1, z_user / gamma_0* num_antenna* p_max)) / log(2);
    
    for k = 1 : num_user
        
        z_user_expanded = repmat(z_user(k,:), num_target, 1);
        user_rate_ISAC_tmp = log(1 + gamma_0 * (num_antenna * p_max - z_target_l * sensing_th) .* inv_pos(z_user_expanded)) / log(2);

        E_tmp = reshape(E(k,:,:), [size(E(k,:,:),2) size(E(k,:,:),3)]);
        user_rate_ISAC_sum(k,:) = sum(E_tmp .* (user_rate_ISAC_tmp - user_rate_comm(k,:)));
    end
    
    user_rate = A .* user_rate_comm + user_rate_ISAC_sum;
end