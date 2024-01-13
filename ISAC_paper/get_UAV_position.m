function uav = get_UAV_position(uav_t, W_opt, R_opt, user, num_user, channel_gain, noise_power, sensing_th, num_target, target, trust_region)
    
    c = zeros(num_user, 1);
    c_1_tmp = zeros(num_user, 1);
    c_2_tmp = zeros(num_user, 1);

    d =zeros(1, 2, num_user);
    d_1_tmp = zeros(1, 2, num_user);
    d_2_tmp = zeros(1, 2, num_user);

    h_target = zeros(num_target, 1);
    i_target = zeros(1, 2, num_target);
    distance_target_UAV_t = zeros(num_target, 1);

    G_opt = 0;

    for k = 1:num_user
        G_opt = G_opt + W_opt(:,:,k);
    end

    G_opt = G_opt + R_opt;
    
    cvx_begin
    
        cvx_solver Mosek

        variable uav(1,2)
        expressions user_rate(num_user, 1)
        expressions sensing_constraint(num_target, 1)
        expressions distance_target_UAV(num_target, 1)

        for k = 1:num_user
    
            distance_k = get_distance(user(k,:), uav_t);
    
            for i = 1:num_user
                c_1_tmp(k) = c_1_tmp(k) + get_beam_gain_UAV(W_opt(:,:,i), distance_k);
                d_1_tmp(:,:,k) = d_1_tmp(:,:,k) + get_beam_gain_diff_UAV(W_opt(:,:,i), distance_k, uav_t, user(k,:));
    
                if i == k
                    continue;
                end
    
                c_2_tmp(k) = c_2_tmp(k) + get_beam_gain_UAV(W_opt(:,:,i), distance_k);
                d_2_tmp(:,:,k) = d_2_tmp(:,:,k) + get_beam_gain_diff_UAV(W_opt(:,:,i), distance_k, uav_t, user(k,:));
            end
    
            c_1_tmp_tmp = log2(c_1_tmp(k) + get_beam_gain_UAV(R_opt, distance_k) + noise_power * distance_k^2 / channel_gain);
            c_2_tmp_tmp = log2(c_2_tmp(k) + get_beam_gain_UAV(R_opt, distance_k) + noise_power * distance_k^2 / channel_gain);
    
            c(k) = c_1_tmp_tmp - c_2_tmp_tmp;
    
            e_tmp = c_1_tmp(k) + get_beam_gain_UAV(R_opt, distance_k) + noise_power * distance_k^2 / channel_gain;
            f_tmp = c_2_tmp(k) + get_beam_gain_UAV(R_opt, distance_k) + noise_power * distance_k^2 / channel_gain;
    
            d_1_tmp_tmp = d_1_tmp(:,:,k) + get_beam_gain_diff_UAV(R_opt, distance_k, uav_t, user(k,:)) + noise_power * (uav_t - user(k,:)) / channel_gain;
            d_1_tmp_tmp = d_1_tmp_tmp / (e_tmp * log(2));
    
            d_2_tmp_tmp = d_2_tmp(:,:,k) + get_beam_gain_diff_UAV(R_opt, distance_k, uav_t, user(k,:)) + noise_power * (uav_t - user(k,:)) / channel_gain;
            d_2_tmp_tmp = d_2_tmp_tmp / (f_tmp * log(2));
    
            d(:,:,k) = d_1_tmp_tmp - d_2_tmp_tmp;

            user_rate(k) = c(k) + sum(d(:,:,k) .* (uav - uav_t));
        end

        maximize(sum(user_rate));

        subject to
                
            for j = 1:num_target

                distance_target_UAV(j) = get_distance(target(j,:), uav);
                distance_target_UAV_t(j) = get_distance(target(j,:), uav_t);

                h_target(j) = get_beam_gain_UAV(G_opt, distance_target_UAV_t(j));
                i_target(:,:,j) = get_beam_gain_diff_UAV(G_opt, distance_target_UAV_t(j), uav_t, target(j,:));

                sensing_constraint(j) = h_target(j) + sum(i_target(:,:,j) .* (uav - uav_t)) - sensing_th * pow_pos(distance_target_UAV(j),2);
                sensing_constraint(j) >= 0;
            end

            distance_UAV = norm([uav(1,1) - uav_t(1,1), uav(1,2) - uav_t(1,2)]);
            distance_UAV <= trust_region;

            0 <= uav(1,1) <= 1000;
            0 <= uav(1,2) <= 1000;
    cvx_end
end