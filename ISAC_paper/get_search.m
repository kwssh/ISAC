function get_search(W_opt, W_t, R_opt, R_t, p_max, sensing_th, num_target, num_antenna, noise_power, user, target, num_user, scaling)

    sum_rate_xy = zeros(1000, 1000);

    for x = 1:1000

        for y = 1:1000

            uav = [x y];

            for k = 1:num_user
                distance_user(k) = get_distance(uav, user(k,:));
                channel(:,k) = get_channel(uav, user(k,:), scaling);
                channel_her(k,:) = transpose(conj(channel(:,k)));
            end
        
            for j = 1:num_target
                distance_target(j) = get_distance(uav, target(j,:));
                steering_target(:,j) = get_steering(distance_target(j), scaling);
                steering_target_her(j,:) = transpose(conj(steering_target(:,j)));
            end
        
            [sum_rate_current_UAV, ~, ~] = get_test(W_opt, W_t, R_opt, R_t, p_max, sensing_th, num_target, num_antenna, channel, channel_her, noise_power, steering_target, steering_target_her, distance_target);
            
            sum_rate_xy(x,y) = sum(sum_rate_current_UAV);
        end
    end

    disp("qweqwe")
















end