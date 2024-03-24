function uav = get_UAV_trajectory_tmp_no_interference(uav_t, W_opt, R_opt, N, num_user, num_target, noise_power, user, uav_z, target, sensing_th, V_max, delta_t, PARAM)

    distance_user = zeros(PARAM.NUM_USER, PARAM.N);
    distance_user_tmp = zeros(PARAM.NUM_USER, PARAM.N);
    channel = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N);
    channel_her = zeros(PARAM.NUM_USER, PARAM.NUM_ANTENNA, PARAM.N);
    distance_target_t_tmp = zeros(PARAM.NUM_TARGET, PARAM.N);
    steering_target = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_TARGET, PARAM.N);
    steering_target_her = zeros(PARAM.NUM_TARGET, PARAM.NUM_ANTENNA, PARAM.N);

    E = zeros(num_user, N);
    F = zeros(num_user, N);

    for n = 1 : PARAM.N
        for k = 1 : PARAM.NUM_USER
            distance_user_tmp(k, n) = get_distance(uav_t(n, :), PARAM.USER(k,:), PARAM.UAV_Z);
            distance_user(k, n) = distance_user_tmp(k, n).^2;
            channel(:, k, n) = get_channel_tmp(uav_t(n, :), PARAM.USER(k,:), PARAM.SCALING * 10, PARAM.UAV_Z, PARAM.NUM_ANTENNA);
            channel_her(k, :, n) = transpose(conj(channel(:, k, n)));
        end
    
        for j = 1 : PARAM.NUM_TARGET
            distance_target_t_tmp(j, n) = get_distance(uav_t(n, :), PARAM.TARGET(j,:), PARAM.UAV_Z);
            steering_target(:, j, n) = get_steering(distance_target_t_tmp(j, n), PARAM.SCALING * 10, PARAM.UAV_Z, PARAM.NUM_ANTENNA);
            steering_target_her(j, :, n) = transpose(conj(steering_target(:, j, n)));
        end
    end

    cvx_begin
    
        cvx_solver Mosek

        variable uav(N, 2)
        variable eta(num_user, N)

        expressions distance_user_uav(num_user, N)
        expressions distance_user_uav(num_target, N)
        expressions sensing_power(num_target, N)
        expressions gamma_tilda(num_user, N)
        expressions gamma_low(num_user, N)
        
        for n = 1 : N

            W_sum = sum(W_opt(:,:,1:num_user,n), 3);
    
            for k = 1 : num_user

                distance_user_uav(k,n) = get_distance(uav(n, :), user, uav_z);
    
                E(k,n) = real(trace(channel(:, k, n) * channel_her(k, :, n) * W_opt(:,:,k,n)));
    
                for i = 1 : num_user
    
                    if i == k
                        continue;
                    end
    
                    % F(k,n) = F(k,n) + real(trace(channel(:, k, n) * channel_her(k, :, n) * W_opt(:,:,i,n)));
    
                end
    
                % F(k,n) = F(k,n) + real(trace(channel(:, k, n) * channel_her(k, :, n) * R_opt(:,:,n)));
    
                first_val = log2((E(k,n) + F(k,n)) / distance_user(k,n) + PARAM.NOISE_POWER_SCALING * 100);
                first_dev_tmp1 = -log2(exp(1)) * (E(k,n) + F(k,n)) / distance_user(k,n)^2;
                first_dev_tmp2 = (E(k,n) + F(k,n)) / distance_user(k,n) + PARAM.NOISE_POWER_SCALING * 100;
                first_dev = first_dev_tmp1 / first_dev_tmp2;

                gamma_low(k,n) = first_val + first_dev * (distance_user_uav(k,n) - distance_user(k,n));

                % gamma_tilda(k,n) = gamma_low(k,n) - (log_sum_exp([log(F(k,n)) + eta(k,n) log(PARAM.NOISE_POWER_SCALING)])) / log(2);

                gamma_tilda(k,n) = gamma_low(k,n) - log2(PARAM.NOISE_POWER_SCALING * 100);
                
                subject to
                    PARAM.SCALING * (1 / exp(eta(k,n))) <= PARAM.SCALING * (norm([uav_t(n,1) - user(k,1), uav_t(n,2) - user(k,2)])^2 + sum(2 * (uav_t(n,:) - user(k,:)) .* (uav(n,:) - uav_t(n,:))));
            end

            for j = 1 : num_target
                distance_target_uav(j,n) = get_distance(uav(n, :), target, uav_z);
                sensing_power(j, n) = real(steering_target_her(j,:,n) * (W_sum + R_opt(:,:,n)) * steering_target(:,j,n));
                sensing_power(j, n) >= distance_target_uav(j,n) * sensing_th;
            end

            -1000 <= uav(n,1) <= 1000;
            -1000 <= uav(n,2) <= 1000;

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

        maximize(sum(sum(gamma_tilda)));

    cvx_end
end