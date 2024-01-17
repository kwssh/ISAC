function sum_rate_final = ISAC_MRT()
    clear
    format long
    rng(123);

    %-----------------------------setting parameter-----------------------------------------------------------------------------------------------------------------------------%
    PARAM.SCALING = 1000;

    PARAM.NUM_USER = 1;
    PARAM.NUM_TARGET = 0;
    PARAM.NUM_ANTENNA = 12;
    PARAM.NUM_EPISODE = 100;

    PARAM.USER = [370 400];
    PARAM.UAV_START = [450 525];
    PARAM.UAV_END = [550 525];
    PARAM.UAV_Z = 1;
    PARAM.TARGET = [randi([450, 550], PARAM.NUM_TARGET, 1) randi([590, 610], PARAM.NUM_TARGET, 1)];

    PARAM.NOISE_POWER = 10^-14;
    PARAM.NOISE_POWER_SCALING = PARAM.NOISE_POWER  * PARAM.SCALING^2;

    PARAM.SENSING_TH = 5 * 10^-5 * 0;
    PARAM.SENSING_TH_SCALING = PARAM.SENSING_TH * PARAM.SCALING^2;

    PARAM.P_MAX = 0.5;
    PARAM.CHANNEL_GAIN = 10^(-6);

    PARAM.T = 15;
    PARAM.N = 3;
    PARAM.DELTA_T = PARAM.T / PARAM.N;
    PARAM.V_MAX = 1000;
    PARAM.TRUST_REGION = 100;
    %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
    %-----------------------------initializing variable-----------------------------------------------------------------------------------------------------------------------------%
    channel_t = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N);
    channel_her_t = zeros(PARAM.NUM_USER, PARAM.NUM_ANTENNA, PARAM.N);

    channel = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N);
    channel_her = zeros(PARAM.NUM_USER, PARAM.NUM_ANTENNA, PARAM.N);

    steering_target_t = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_TARGET, PARAM.N);
    steering_target_her_t = zeros(PARAM.NUM_TARGET, PARAM.NUM_ANTENNA, PARAM.N);

    steering_target = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_TARGET, PARAM.N);
    steering_target_her = zeros(PARAM.NUM_TARGET, PARAM.NUM_ANTENNA, PARAM.N);

    distance_target_t = zeros(PARAM.NUM_TARGET, PARAM.N);
    distance_user_t = zeros(PARAM.NUM_USER, PARAM.N);

    distance_target = zeros(PARAM.NUM_TARGET, PARAM.N);
    distance_user = zeros(PARAM.NUM_USER, PARAM.N);

    user_rate_episode = zeros(PARAM.NUM_USER, PARAM.N, PARAM.NUM_EPISODE);

    W_opt = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N);
    R_opt = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.N);
    %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
    %-----------------------------Epsiode Start-----------------------------------------------------------------------------------------------------------------------------%   
    for episode = 1 : PARAM.NUM_EPISODE

        if episode == 1
            %-----------------------------initializing precoder-----------------------------------------------------------------------------------------------------------------------------%
            [uav_t, ~, ~] = get_init_trajectory_SF(PARAM.TARGET, PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.NUM_TARGET, PARAM.SENSING_TH_SCALING, PARAM.P_MAX, PARAM.SCALING, PARAM.V_MAX, PARAM.N, PARAM.UAV_START(1), PARAM.UAV_END(1), PARAM.UAV_START(2), PARAM.UAV_Z);
            %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
      
            %-----------------------------get channel and steering------------------------------------------------------------------------------------------------------------------------------%
            for n = 1 : PARAM.N
                for k = 1 : PARAM.NUM_USER
                    distance_user_t(k, n) = get_distance(uav_t(n, :), PARAM.USER(k,:), PARAM.UAV_Z);
                    channel_t(:, k, n) = get_channel(uav_t(n, :), PARAM.USER(k,:), PARAM.SCALING, PARAM.UAV_Z, PARAM.NUM_ANTENNA);
                    channel_her_t(k, :, n) = transpose(conj(channel_t(:, k, n)));

                    W_opt_tmp = (channel_t(:,k,n) / (sqrt(channel_her_t(k,:,n) * channel_t(:,k,n)))) / sqrt(2);
                    W_opt(:,:,k,n) = W_opt_tmp * transpose(conj(W_opt_tmp));
                end
            
                for j = 1 : PARAM.NUM_TARGET
                    distance_target_t(j, n) = get_distance(uav_t(n, :), PARAM.TARGET(j,:), PARAM.UAV_Z);
                    steering_target_t(:, j, n) = get_steering(distance_target_t(j, n), PARAM.SCALING, PARAM.UAV_Z, PARAM.NUM_ANTENNA);
                    steering_target_her_t(j, :, n) = transpose(conj(steering_target_t(:, j, n)));
                end
            end
            %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
            [user_rate_current, ~] = get_test_trajectory(W_opt, R_opt, PARAM.P_MAX, PARAM.SENSING_TH_SCALING, PARAM.NUM_TARGET, channel_t, channel_her_t, PARAM.NOISE_POWER_SCALING, steering_target_t, steering_target_her_t, distance_target_t, PARAM.N);
            user_rate_episode(:,:,1) = user_rate_current;
        end

        %-----------------------------get MRT precoder------------------------------------------------------------------------------------------------------------------------------%
        for n = 1 : PARAM.N
            for k = 1 : PARAM.NUM_USER
                W_opt_tmp = (channel_t(:,k,n) / (sqrt(channel_her_t(k,:,n) * channel_t(:,k,n)))) / sqrt(2);
                W_opt(:,:,k,n) = W_opt_tmp * transpose(conj(W_opt_tmp));
            end
        end
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
        user_rate_prev = user_rate_current;

        [user_rate_prev_UAV, error_prev_UAV] = get_test_trajectory(W_opt, R_opt, PARAM.P_MAX, PARAM.SENSING_TH_SCALING, PARAM.NUM_TARGET, channel_t, channel_her_t, PARAM.NOISE_POWER_SCALING, steering_target_t, steering_target_her_t, distance_target_t, PARAM.N);
        
        %-----------------------------optimize UAV-----------------------------------------------------------------------------------------------------------------------------%
       
        trust_region = PARAM.TRUST_REGION;

        while(1)

            uav = get_UAV_trajectory(uav_t, W_opt, R_opt, PARAM.USER, PARAM.NUM_USER, PARAM.CHANNEL_GAIN, PARAM.NOISE_POWER, PARAM.SENSING_TH, PARAM.NUM_TARGET, PARAM.TARGET, trust_region, PARAM.V_MAX, PARAM.N, PARAM.UAV_Z, PARAM.DELTA_T);
   
            %-----------------------------get channel and steering with new UAV-----------------------------------------------------------------------------------------------------------------------------%
            for n = 1 : PARAM.N
                for k = 1 : PARAM.NUM_USER
                    distance_user(k, n) = get_distance(uav(n, :), PARAM.USER(k,:), PARAM.UAV_Z);
                    channel(:, k, n) = get_channel(uav(n, :), PARAM.USER(k,:), PARAM.SCALING, PARAM.UAV_Z, PARAM.NUM_ANTENNA);
                    channel_her(k, :, n) = transpose(conj(channel(:, k, n)));
                end
            
                for j = 1 : PARAM.NUM_TARGET
                    distance_target(j, n) = get_distance(uav(n, :), PARAM.TARGET(j,:), PARAM.UAV_Z);
                    steering_target(:, j, n) = get_steering(distance_target(j, n), PARAM.SCALING, PARAM.UAV_Z, PARAM.NUM_ANTENNA);
                    steering_target_her(j, :, n) = transpose(conj(steering_target(:, j, n)));
                end
            end
            %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

            [user_rate_current_UAV, error_current_UAV] = get_test_trajectory(W_opt, R_opt, PARAM.P_MAX, PARAM.SENSING_TH_SCALING, PARAM.NUM_TARGET, channel, channel_her, PARAM.NOISE_POWER_SCALING, steering_target, steering_target_her, distance_target, PARAM.N);
                                                          
            %-----------------------------display part-----------------------------------------------------------------------------------------------------------------------------%
            disp(['Episode           : ', num2str(episode)]);
            disp(['Trust Region      : ', num2str(trust_region)]);
            disp(['Previous Sum rate : ', num2str(sum(user_rate_prev_UAV))]);
            disp(['Current Sum rate  : ', num2str(sum(user_rate_current_UAV))]);
            disp(['Diff Sum rate     : ', num2str(sum(user_rate_current_UAV) - sum(user_rate_prev_UAV))]);
            get_display(distance_user, 'Distance user     : ');
            get_display(uav_t, 'UAV position      : ');
            %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

            if sum(sum(user_rate_current_UAV)) > sum(sum(user_rate_prev_UAV))
                uav_t = uav;
                trust_region = PARAM.TRUST_REGION;

                user_rate_prev_UAV = user_rate_current_UAV;
            else
                trust_region = trust_region / 2;
            end

            if trust_region < 10^(-6)
                break
            end
        end

        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

        for n = 1 : PARAM.N
            for k = 1 : PARAM.NUM_USER
                distance_user_t(k, n) = get_distance(uav_t(n, :), PARAM.USER(k,:), PARAM.UAV_Z);
                channel_t(:, k, n) = get_channel(uav_t(n, :), PARAM.USER(k,:), PARAM.SCALING, PARAM.UAV_Z, PARAM.NUM_ANTENNA);
                channel_her_t(k, :, n) = transpose(conj(channel_t(:, k, n)));
            end
        
            for j = 1 : PARAM.NUM_TARGET
                distance_target_t(j, n) = get_distance(uav_t(n, :), PARAM.TARGET(j,:), PARAM.UAV_Z);
                steering_target_t(:, j, n) = get_steering(distance_target_t(j, n), PARAM.SCALING, PARAM.UAV_Z, PARAM.NUM_ANTENNA);
                steering_target_her_t(j, :, n) = transpose(conj(steering_target_t(:, j, n)));
            end
        end

        [user_rate_current, ~] = get_test_trajectory(W_opt, R_opt, PARAM.P_MAX, PARAM.SENSING_TH_SCALING, PARAM.NUM_TARGET, channel_t, channel_her_t, PARAM.NOISE_POWER_SCALING, steering_target_t, steering_target_her_t, distance_target_t, PARAM.N);
       
        user_rate_episode(:,:,episode) = user_rate_current;

        if sum(sum(user_rate_current)) <= sum(sum(user_rate_prev))
            break;
        end

    end
    %-----------------------------Epsiode End-----------------------------------------------------------------------------------------------------------------------------%
    
    % get_received_BEAM_GAIN(W_opt, R_opt, num_user, num_antenna, sensing_th, scaling, distance_user_t, distance_target_t);

end