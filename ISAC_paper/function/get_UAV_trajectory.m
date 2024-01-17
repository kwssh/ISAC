function uav = get_UAV_trajectory(uav_t, W_opt, R_opt, user, num_user, channel_gain, noise_power, sensing_th, num_target, target, trust_region, V_max, N, uav_z, delta_t)
    
    num_antenna = size(W_opt, 1);
    c = zeros(num_user, N);
    c_1_tmp = zeros(num_user, N);
    c_2_tmp = zeros(num_user, N);

    d =zeros(1, 2, num_user, N);
    d_1_tmp = zeros(1, 2, num_user, N);
    d_2_tmp = zeros(1, 2, num_user, N);

    h_target = zeros(num_target, N);
    i_target = zeros(1, 2, num_target, N);
    distance_target_UAV_t = zeros(num_target, N);

    G_opt = zeros(num_antenna, num_antenna, N);

    for n = 1:N
        for k = 1:num_user
            G_opt(:,:,n) = G_opt(:,:,n) + W_opt(:,:,k,n);
        end
    end

    G_opt = G_opt + R_opt;
    
    cvx_begin
    
        cvx_solver Mosek

        variable uav(N, 2)
        expressions user_rate(num_user, N)
        expressions sensing_constraint(num_target, N)
        expressions distance_target_UAV(num_target, N)

        for n = 1:N

            for k = 1:num_user
        
                distance_k = get_distance(user(k,:), uav_t(n,:), uav_z);
        
                for i = 1:num_user
                    c_1_tmp(k,n) = c_1_tmp(k,n) + get_beam_gain_UAV(W_opt(:,:,i,n), distance_k);
                    d_1_tmp(:,:,k,n) = d_1_tmp(:,:,k,n) + get_beam_gain_diff_UAV(W_opt(:,:,i,n), distance_k, uav_t(n,:), user(k,:));
        
                    if i == k
                        continue;
                    end
        
                    c_2_tmp(k,n) = c_2_tmp(k,n) + get_beam_gain_UAV(W_opt(:,:,i,n), distance_k);
                    d_2_tmp(:,:,k,n) = d_2_tmp(:,:,k,n) + get_beam_gain_diff_UAV(W_opt(:,:,i,n), distance_k, uav_t(n,:), user(k,:));
                end
        
                c_1_tmp_tmp = log2(c_1_tmp(k,n) + get_beam_gain_UAV(R_opt(:,:,n), distance_k) + noise_power * distance_k^2 / channel_gain);
                c_2_tmp_tmp = log2(c_2_tmp(k,n) + get_beam_gain_UAV(R_opt(:,:,n), distance_k) + noise_power * distance_k^2 / channel_gain);
        
                c(k,n) = c_1_tmp_tmp - c_2_tmp_tmp;
        
                e_tmp = c_1_tmp(k,n) + get_beam_gain_UAV(R_opt(:,:,n), distance_k) + noise_power * distance_k^2 / channel_gain;
                f_tmp = c_2_tmp(k,n) + get_beam_gain_UAV(R_opt(:,:,n), distance_k) + noise_power * distance_k^2 / channel_gain;
        
                d_1_tmp_tmp = d_1_tmp(:,:,k,n) + get_beam_gain_diff_UAV(R_opt(:,:,n), distance_k, uav_t(n,:), user(k,:)) + noise_power * (uav_t(n,:) - user(k,:)) / channel_gain;
                d_1_tmp_tmp = d_1_tmp_tmp / (e_tmp * log(2));
        
                d_2_tmp_tmp = d_2_tmp(:,:,k,n) + get_beam_gain_diff_UAV(R_opt(:,:,n), distance_k, uav_t(n,:), user(k,:)) + noise_power * (uav_t(n,:) - user(k,:)) / channel_gain;
                d_2_tmp_tmp = d_2_tmp_tmp / (f_tmp * log(2));
        
                d(:,:,k,n) = d_1_tmp_tmp - d_2_tmp_tmp;
    
                user_rate(k,n) = c(k,n) + sum(d(:,:,k,n) .* (uav(n,:) - uav_t(n,:)));
            end
    
            subject to
                    
                for j = 1:num_target
    
                    distance_target_UAV(j) = get_distance(target(j,:), uav(n,:), uav_z);
                    distance_target_UAV_t(j,n) = get_distance(target(j,:), uav_t(n,:), uav_z);
    
                    h_target(j,n) = get_beam_gain_UAV(G_opt(:,:,n), distance_target_UAV_t(j,n));
                    i_target(:,:,j,n) = get_beam_gain_diff_UAV(G_opt(:,:,n), distance_target_UAV_t(j,n), uav_t(n,:), target(j,:));
    
                    sensing_constraint = h_target(j,n) + sum(i_target(:,:,j,n) .* (uav(n,:) - uav_t(n,:))) - sensing_th * pow_pos(distance_target_UAV(j),2);
                    sensing_constraint >= 0;
                end
    
                distance_UAV = norm([uav(n,1) - uav_t(n,1), uav(n,2) - uav_t(n,2)]);
                distance_UAV <= trust_region;
    
                0 <= uav(n,1) <= 1000;
                0 <= uav(n,2) <= 1000;

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
        
        maximize(sum(sum(user_rate)));

    cvx_end
end