function user_rate = get_user_rate(gamma_0, num_antenna, p_max, distance_user, num_user, distance_target, num_target, E, sensing_th, A, user_rate_ISAC_sum)

    user_rate_comm = log2(1 + (gamma_0 * num_antenna * p_max ./ (distance_user .* distance_user)));
    
    for k = 1 : num_user
        
        user_rate_ISAC_tmp = log2(1 + gamma_0 * (num_antenna * p_max - distance_target.^2 * sensing_th) ./ (distance_user(k,:).^2));

        E_tmp = reshape(E(k,:,:), [size(E(k,:,:),2) size(E(k,:,:),3)]);
        user_rate_ISAC_sum(k,:) = sum(E_tmp .* (user_rate_ISAC_tmp - user_rate_comm(k,:)));
    end
    
    user_rate = A .* user_rate_comm + user_rate_ISAC_sum;
end