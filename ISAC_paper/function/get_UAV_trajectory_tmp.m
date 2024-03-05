function uav = get_UAV_trajectory_tmp(uav_t, W_opt, R_opt, channel, channel_her, steering, steering_her, N, num_user, num_target, distance_user, noise_power, user, uav_z)

    E = zeros(num_user, N);
    F = zeros(num_user, N);
    
    cvx_begin
    
        cvx_solver Mosek

        variable uav(N, 2)
        variable eta(num_user, N)

        expressions distance_uav(num_user, N)
        expressions gamma_tilda(num_user, N)
        expressions gamma_low(num_user, N)
        
        for n = 1 : N
    
            for k = 1 : num_user

                distance_uav(k,n) = get_distance(uav(n, :), user, uav_z);
    
                E(k,n) = real(trace(channel(:, k, n) * channel_her(k, :, n) * W_opt(:,:,k,n)));
    
                for i = 1 : num_user
    
                    if i == k
                        continue;
                    end
    
                    F(k,n) = F(k,n) + real(trace(channel(:, k, n) * channel_her(k, :, n) * W_opt(:,:,i,n)));
    
                end
    
                F(k,n) = F(k,n) + real(trace(channel(:, k, n) * channel_her(k, :, n) * R_opt(:,:,n)));
    
                first_val = log2((E(k,n) + F(k,n)) / distance_user(k,n) + noise_power);
                first_dev_tmp1 = -log2(exp(1)) * (E(k,n) + F(k,n)) / distance_user(k,n)^2;
                first_dev_tmp2 = (E(k,n) + F(k,n)) / distance_user(k,n) + noise_power;
                first_dev = first_dev_tmp1 / first_dev_tmp2;

                gamma_low(k,n) = first_val + first_dev * (distance_uav(k,n) - distance_user(k,n));

                gamma_tilda(k,n) = gamma_low(k,n) - log_sum_exp([log(F(k,n)) + eta(k,n) log(noise_power)]);
            end
    
    
        end

end