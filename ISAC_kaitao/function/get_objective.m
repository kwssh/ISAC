function user_rate = get_objective(distance_user, distance_target, num_user, num_target, N, gamma_0, p_max, num_antenna, sensing_th, A_opt, E_opt)

    user_rate = zeros(num_user, N);
    user_rate_comm = zeros(num_user, N);
    user_rate_ISAC = zeros(num_user, N);
    user_rate_ISAC_tmp = zeros(num_user, num_target, N);

    for n = 1 : N
        for k = 1 : num_user

            user_rate_comm(k,n) = log2(1 + gamma_0 * num_antenna * p_max / (distance_user(k,n)^2));

            for j = 1 : num_target

                user_rate_ISAC_tmp(k,j,n) = log2(1 + gamma_0 * (num_antenna * p_max - distance_target(j,n)^2 * sensing_th) / (distance_user(k,n)^2));

                user_rate_ISAC(k,n) = user_rate_ISAC(k,n) + E_opt(k,j,n) * (user_rate_ISAC_tmp(k,j,n) - user_rate_comm(k,n));
            end

            user_rate(k,n) = A_opt(k,n) * user_rate_comm(k,n) + user_rate_ISAC(k,n);

        end
    end
    

end