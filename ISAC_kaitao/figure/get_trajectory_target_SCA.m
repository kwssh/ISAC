function z_target = get_trajectory_target_SCA(uav, z_user, gamma_0, num_antenna, p_max, A_opt, E_opt, N, num_user, num_target, sensing_th, PARAM, isac_duration, rate_th)

    user_rate_comm = log2(1+gamma_0 * num_antenna * p_max ./ z_user);

    cvx_begin

        cvx_solver Mosek

        variable z_target(num_target, N)

        expressions user_rate(num_user, N)
        expressions user_rate_ISAC(num_user, N)
        expressions user_rate_ISAC_target(num_user, num_target, N)

        for n = 1 : N
            for k = 1 : num_user
                
                for j = 1 : num_target
                    user_rate_ISAC_target(k,j,n) = log(1 + gamma_0 * (num_antenna * p_max - z_target(j,n) * sensing_th) / z_user(k,n));

                    user_rate_ISAC(k,n) = user_rate_ISAC(k,n) + E_opt(k,j,n) * (user_rate_ISAC_target(k,j,n) - user_rate_comm(k,n));

                    z_target(j,n) >= pow_pos(norm([PARAM.TARGET(j, 1) - uav(n, 1), PARAM.TARGET(j, 2) - uav(n, 2), PARAM.UAV_Z]), 2);
                end

                user_rate(k,n) = A_opt(k,n) * user_rate_comm(k,n) + user_rate_ISAC(k,n);
            end
        end

        maximize(sum(sum(user_rate)))

        for i = 1:isac_duration:N
            sum(user_rate(:,i:i+isac_duration-1),2)/ isac_duration >= rate_th;
        end
    cvx_end
end