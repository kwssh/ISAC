function [uav, z_user, user_rate] = get_trajectory_user_SCA(distance_user, distance_target, num_user, num_target, N, gamma_0, p_max, num_antenna, sensing_th, PARAM, uav_t, V_max, delta_t, A_opt, E_opt, rate_th, isac_duration)

    num_episode_SCA = 10^6;
    user_rate_episode_SCA = zeros(num_user, N, num_episode_SCA);

    z_user_l = distance_user.^2;
    z_target_l = distance_target.^2;

    for episode_SCA = 1 : num_episode_SCA

        cvx_begin
    
            cvx_solver Mosek
    
            variable z_user(num_user, N)
            variable uav(N, 2)
    
            expressions user_rate(num_user, N)
            expressions user_rate_comm(num_user, N)
            expressions user_rate_ISAC(num_user, N)
            expressions user_rate_ISAC_target(num_user, num_target, N)
            expressions sensing_constraint(num_target, N)
            expressions distance_target_new(num_target, N)
            expressions distance_user_new(num_user, N)
            expressions user_rate_th(num_user, N)
            expressions user_rate_ISAC_sum(num_user, N)
    
            distance_target_new = get_distance_cvx(PARAM.TARGET, uav, PARAM.UAV_Z, distance_target_new);

            for n = 1 : N
                for k = 1 : num_user
    
                    user_rate_comm_tmp = -gamma_0 * num_antenna * p_max * (z_user(k,n) - z_user_l(k,n)) / (z_user_l(k,n)^2 + gamma_0 * num_antenna * p_max * z_user_l(k,n)) / log(2);
                    user_rate_comm(k,n) = user_rate_comm_tmp + log2(1 + gamma_0 * num_antenna * p_max / z_user_l(k,n));
    
                    for j = 1 : num_target
                        
                        user_rate_ISAC_tmp = -gamma_0 * (num_antenna * p_max - z_target_l(j,n) * sensing_th) * (z_user(k,n) - z_user_l(k,n)) / (z_user_l(k,n)^2 + gamma_0 * (num_antenna * p_max - z_target_l(j,n) * sensing_th) * z_user_l(k,n)) / log(2);
                        user_rate_ISAC_target(k,j,n) = user_rate_ISAC_tmp + log2(1 + gamma_0 * (num_antenna * p_max - z_target_l(j,n) * sensing_th) / z_user_l(k,n));
    
                        user_rate_comm_tmp_tmp = (rel_entr(z_user(k,n) / (gamma_0* num_antenna* p_max), z_user(k,n) / (gamma_0* num_antenna* p_max) + 1) + rel_entr(z_user(k,n) / (gamma_0* num_antenna* p_max) + 1, z_user(k,n) / (gamma_0* num_antenna* p_max))) / log(2);
                        
                        user_rate_ISAC(k,n) = user_rate_ISAC(k,n) + E_opt(k,j,n) * (user_rate_ISAC_target(k,j,n) - user_rate_comm_tmp_tmp);
    
                    end
    
                    user_rate(k,n) = A_opt(k,n) * user_rate_comm(k,n) + user_rate_ISAC(k,n);
    
                    z_user(k,n) >= pow_pos(norm([PARAM.USER(k, 1) - uav(n, 1), PARAM.USER(k, 2) - uav(n, 2), PARAM.UAV_Z]), 2);
    
                end

                for j = 1 : num_target
                    for i = 1 : num_user
                        sensing_constraint(j,n) = sensing_constraint(j,n) + E_opt(i,j,n) * (num_antenna * p_max - pow_pos(distance_target_new(j,n), 2) * sensing_th);
                    end
                end

                sensing_constraint >= 0;
    
                if n == 1
                    uav(1,1) == uav_t(1,1);
                    uav(1,2) == uav_t(1,2);
    
                    uav(N,1) == uav_t(N,1);
                    uav(N,2) == uav_t(N,2);
                else
                    velocity_UAV = norm([uav(n,1) - uav(n-1,1), uav(n,2) - uav(n-1,2)]);
                    velocity_UAV <= V_max * delta_t;
                end
            end

            maximize(sum(sum(user_rate)))

            for i = 1:isac_duration:N
                sum(user_rate(:,i:i+isac_duration-1),2)/ isac_duration >= rate_th;
            end
    
        cvx_end

        user_rate_episode_SCA(:,:,episode_SCA) = user_rate;

        % break

        if episode_SCA > 1
            if abs(sum(sum(user_rate_episode_SCA(:,:,episode_SCA))) - sum(sum(user_rate_episode_SCA(:,:,episode_SCA-1)))) <= 0.01
                break
            end
        end
        
        z_user_l = z_user;
    end
end