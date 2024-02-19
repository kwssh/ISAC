function uav = get_uav_trajectory(A_opt, E_opt, A_bar_opt, E_bar_opt, distance_user, distance_target, num_user, num_target, N, gamma_0, p_max, num_antenna, sensing_th)

    A = gamma_0 * p_max * num_antenna;
    
    cvx_begin

        cvx_solver Mosek

        variable z_user(num_user, N)
        variable z_target(num_target, N)

        expressions user_rate_ISAC(num_user, N)

        user_rate_comm_first = log2(1 + A ./ distance_user) - A * (z_user - distance_user) ./ ((distance_user.^2 + A * distance_user) * log(2));
        user_rate_comm = (rel_entr(A + z_user, z_user) + rel_entr(z_user, A + z_user)) * log2(exp(1)) / A;
        
        for k = 1 : num_user

            z_user_tmp = repmat(z_user(k, :), num_target, 1);
        
            user_rate_ISAC_tmp = log2(1 + ((gamma_0 * num_antenna * p_max - (z_target * sensing_th)) * inv_pos(z_user_tmp)));
    
            user_rate_ISAC_sum(k,:) = sum(E(num_target * k - num_target + 1 : num_target * k,:) .* (user_rate_ISAC_tmp - user_rate_comm(k,:)));
        end



end