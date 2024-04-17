function uav = get_uav_trajectory_BCD(distance_user, distance_target, num_user, num_target, N, gamma_0, p_max, num_antenna, sensing_th)


    cvx_begin

        cvx_solver Mosek

        variable z_user(num_user, N)
        expressions user_rate_ISAC_user(num_user, num_target, N)
        
        for n = 1 : N
            for k = 1 : num_user
                for j = 1 : num_target
                    user_rate_ISAC_user(k,j,n) = log(1 + gamma_0 * (num_antenna * p_max - distance_target(j,n) * sensing_th) * inv_pos(z_user(k,n))) / log(2);
                end

                1000 * z_user(k,n) >= 1000 * distance_user(k,n);
            end
        end

        maximize(sum(sum(sum(user_rate_ISAC_user))))

    cvx_end

    cvx_begin

        cvx_solver Mosek

        variable z_target(num_target, N)
        expressions user_rate_ISAC_target(num_user, num_target, N)
        
        for n = 1 : N
            for k = 1 : num_user
                for j = 1 : num_target
                    user_rate_ISAC_target(k,j,n) = log(1 + gamma_0 * (num_antenna * p_max - z_target(j,n) * sensing_th) * inv_pos(distance_user(k,n))) / log(2);

                    z_target(j,n) >= distance_target(j,n);
                end

               
            end
        end

        maximize(sum(sum(sum(user_rate_ISAC_target))))

    cvx_end

            



end