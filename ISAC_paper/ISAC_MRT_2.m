function sum_rate_final = ISAC_paper_Mobile_UAV_comm()
    clear
    format long
    rng(123);

    scaling = 1000;
    
    % num_target = 18;
    num_target = 0;

    % num_user = 8;
    num_user = 1;

    num_antenna = 12;
    noise_power = 10^-14 * scaling^2;
    p_max = 0.5;
    channel_gain = 10^(-6);

    % sensing_th = 10^(-4.3) * scaling^2;
    sensing_th = 0 * scaling^2;

    V_max = 30000000000000;
    N = 1;

    num_episdoe = 10^6;
    
    % user = [370 400;
    %         380 340; 
    %         420 300; 
    %         470 270; 
    %         530 270; 
    %         580 300; 
    %         620 340; 
    %         630 400];

    user = [370 400];
    
    start_position = [450, 525];
    end_position = [550 525];
    
    target_x_1 = randi([450, 550], num_target / 2, 1);
    target_x_2 = 1000 - target_x_1;
   
    target_y = randi([590, 610], num_target / 2, 1);
    
    target_1 = [target_x_1 target_y];
    target_2 = [target_x_2 target_y];
   
    target = [target_1; target_2];
    
    channel_t = zeros(num_antenna, num_user, N);
    channel_her_t = zeros(num_user, num_antenna, N);

    channel = zeros(num_antenna, num_user, N);
    channel_her = zeros(num_user, num_antenna, N);

    steering_target_t = zeros(num_antenna, num_target, N);
    steering_target_her_t = zeros(num_target, num_antenna, N);
    steering_target = zeros(num_antenna, num_target, N);
    steering_target_her = zeros(num_target, num_antenna, N);

    distance_target_t = zeros(num_target, N);
    distance_user_t = zeros(num_user, N);
    distance_target = zeros(num_target, N);
    distance_user = zeros(num_user, N);

    sum_rate_episode = zeros(num_user, N, num_episdoe);

    W_opt = zeros(num_antenna, num_antenna, num_user, N);
    R_opt = zeros(num_antenna, num_antenna, N);

    %-------------------------------------------------------------------------------------------------------------------------------------------------------------

    for episode = 1:num_episdoe

        if episode == 1

            % [uav_t, W_t, R_t] = get_init_trajectory(target, num_antenna, num_user, num_target, sensing_th, p_max, scaling, V_max, N, start_position, end_position);
            % [uav_t, W_t, R_t] = get_init_trajectory_SF(target, num_antenna, num_user, num_target, sensing_th, p_max, scaling, V_max, N, start_position(1), end_position(1), start_position(2));

            uav_t = [450, 525];
            sum_rate_current = 0;
        end

        sum_rate_prev = sum_rate_current;

        for n = 1 : N

            for k = 1:num_user
                distance_user_t(k, n) = get_distance(uav_t(n, :), user(k,:));
                channel_t(:,k,n) = get_channel(uav_t(n, :), user(k,:), scaling);
                channel_her_t(k,:,n) = transpose(conj(channel_t(:,k,n)));

                W_opt_tmp = (channel_t(:,k,n) / (sqrt(channel_her_t(k,:,n) * channel_t(:,k,n)))) / sqrt(2);
                W_opt(:,:,k,n) = W_opt_tmp * transpose(conj(W_opt_tmp));

            end
        
        end

%-------------------------------------------------------------------------------------------------------------------------------------------------------%

        % get_search(W_opt, W_t, R_opt, R_t, p_max, sensing_th, num_target, num_antenna, noise_power, user, target, num_user, scaling);

        [sum_rate_prev_UAV, error_prev_UAV] = get_test_trajectory(W_opt, R_opt, p_max, sensing_th, num_target, channel_t, channel_her_t, noise_power, steering_target_t, steering_target_her_t, distance_target_t, N);
     
        trust_region = 0.5;

        while(1)

            % uav = get_UAV_trajectory_2(uav_t *scaling, W_opt, R_opt, user *scaling, num_user, channel_gain, noise_power, sensing_th, num_target, target, trust_region, V_max, N, scaling);
            uav = get_UAV_trajectory(uav_t, W_opt, R_opt, user, num_user, channel_gain, noise_power / scaling^2, sensing_th / scaling^2, num_target, target, trust_region, V_max, N);
            
            % disp(uav)
            % disp(q)

            disp(episode);
            disp(trust_region);
            disp(sum(sum(sum_rate_prev_UAV)));

            for n = 1:N
            
                for k = 1:num_user
                    distance_user(k, n) = get_distance(uav(n, :), user(k, :));
                    channel(:, k, n) = get_channel(uav(n, :), user(k, :), scaling);
                    channel_her(k, :, n) = transpose(conj(channel(:, k, n)));
                end
    
                for j = 1:num_target
                    distance_target(j, n) = get_distance(uav(n, :), target(j, :));
                    steering_target(:, j, n) = get_steering(distance_target(j, n), scaling);
                    steering_target_her(j, :, n) = transpose(conj(steering_target(:, j, n)));
                end
            end

            [sum_rate_current_UAV, error_current_UAV] = get_test_trajectory(W_opt, R_opt, p_max, sensing_th, num_target, channel, channel_her, noise_power, steering_target, steering_target_her, distance_target, N);
                                                           
            disp(sum(sum(sum_rate_current_UAV)));
            disp(sum(sum(sum_rate_current_UAV)) - sum(sum(sum_rate_prev_UAV)));
            disp(uav_t);

            if sum(sum(sum_rate_current_UAV)) > sum(sum(sum_rate_prev_UAV))
                uav_t = uav;
                trust_region = 0.5;

                sum_rate_prev_UAV = sum_rate_current_UAV;
            else
                trust_region = trust_region / 2;
            end

            if trust_region < 10^(-6)
                break
            end
        end

%--------------------------------------------------------------------------------------------%

        for n = 1:N
            for k = 1:num_user
                distance_user_t(k,n) = get_distance(uav_t(n,:), user(k,:));
                channel_t(:,k,n) = get_channel(uav_t(n,:), user(k,:), scaling);
                channel_her_t(k,:,n) = transpose(conj(channel_t(:,k,n)));
            end
        
            for j = 1:num_target
                distance_target_t(j,n) = get_distance(uav_t(n,:), target(j,:));
                steering_target_t(:,j,n) = get_steering(distance_target_t(j,n), scaling);
                steering_target_her_t(j,:,n) = transpose(conj(steering_target_t(:,j,n)));
            end
        end

        [sum_rate_current, error] = get_test_trajectory(W_opt, R_opt, p_max, sensing_th, num_target, channel_t, channel_her_t, noise_power, steering_target_t, steering_target_her_t, distance_target_t, N);
     
        disp(error);
            
        sum_rate_episode(:,:,episode) = sum_rate_current;

        if sum(sum(sum_rate_current)) - sum(sum(sum_rate_prev)) < 1e-6
            break;
        end

        R_t = R_opt;
        W_t = W_opt;
    end

    % get_received_BEAM_GAIN(W_opt, R_opt, num_user, num_antenna, sensing_th, scaling, distance_user_t, distance_target_t);

end