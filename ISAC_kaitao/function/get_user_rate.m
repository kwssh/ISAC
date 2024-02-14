function user_rate = get_user_rate(channel_gain, noise_power, num_antenna, p_max, distance_user, num_user, distance_target, num_target, E, sensing_th, A, user_rate_ISAC_sum)

    gamma_0 = channel_gain / noise_power;

    user_rate_comm = log2(1 + (gamma_0 * num_antenna * p_max ./ (distance_user .* distance_user)));
    
    for k = 1 : num_user
        
        user_rate_ISAC_tmp = log2(1 + ((gamma_0 * num_antenna * p_max - (distance_target .* distance_target * sensing_th)) ./ (distance_user(k,:) .* distance_user(k,:))));

        user_rate_ISAC_sum(k,:) = sum(E(num_target * k - num_target + 1 : num_target * k,:) .* (user_rate_ISAC_tmp - user_rate_comm(k,:)));
    end
    
    user_rate = A .* user_rate_comm + user_rate_ISAC_sum;
end