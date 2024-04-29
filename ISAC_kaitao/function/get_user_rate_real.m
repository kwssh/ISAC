function user_rate = get_user_rate_real(distance_user, distance_target, num_user, num_target, N, gamma_0, p_max, num_antenna, sensing_th, A_opt, E_opt, uav, PARAM, channel_gain)

    user_rate = zeros(num_user, N);
    user_rate_tmp = zeros(num_user, N);
    user_rate_comm = zeros(num_user, N);
    user_rate_ISAC = zeros(num_user, num_target, N);
    
    psi = zeros(num_user, num_target, N);
    gamma_tilda = zeros(num_user, num_target, N);
    gamma_tmp = zeros(num_user, num_target, N);

    user_rate_comm = log2(1 + gamma_0 * num_antenna * p_max ./ (distance_user.^2));

    steering_user = get_steering(distance_user, num_antenna, uav, PARAM.USER);

    channel_gain_tmp = sqrt(channel_gain ./ (distance_user.^2));
    channel_gain_repeat = repmat(channel_gain_tmp, [1, 1, num_antenna]);
    channel_gain_permute = permute(channel_gain_repeat, [1, 3, 2]);

    channel_user = channel_gain_permute ./ steering_user;

    steering_target = get_steering(distance_target, num_antenna, uav, PARAM.TARGET);

    distance_target_repeat = repmat(distance_target, [1, 1, num_antenna]);
    distance_target_permute = permute(distance_target_repeat, [1, 3, 2]);

    channel_target = steering_target ./ distance_target_permute;
    
    for n = 1 : N
        for k = 1 : num_user
            for j = 1 : num_target
                psi(k,j,n) = acos((norm(channel_user(k,:,n) * channel_target(j,:,n)')) / (norm(channel_user(k,:,n)) * norm(channel_target(j,:,n)')));
                gamma_tilda(k,j,n) = num_antenna * p_max * cos(psi(k,j,n))^2 / (distance_target(j,n)^2);

                if gamma_tilda(k,j,n) >= sensing_th
                    user_rate_ISAC(k,j,n) = log2(1+gamma_0 * num_antenna * p_max / (distance_user(k,n)^2));
                else
                    gamma_tmp(k,j,n) = gamma_0 * distance_target(j,n)^2 / (distance_user(k,n)^2) * (sqrt(sensing_th) * cos(psi(k,j,n)) + sqrt(num_antenna * p_max / distance_target(j,n)^2-sensing_th) * sin(psi(k,j,n)))^2;
                    user_rate_ISAC(k,j,n) = log2(1+gamma_tmp(k,j,n));
                end

                user_rate_tmp(k,n) = user_rate_tmp(k,n) + E_opt(k,j,n) * (user_rate_ISAC(k,j,n) - user_rate_comm(k,n));
    
            end

            user_rate(k,n) = A_opt(k,n) * user_rate_comm(k,n) + user_rate_tmp(k,n);
        end
    end
end