function sum_rate_final = ISAC_paper_Mobile_UAV()
    clear
    format long
    %-----------------------------setting parameter-----------------------------------------------------------------------------------------------------------------------------%
    PARAM.SCALING = 1000;

    PARAM.NUM_USER = 8;
    PARAM.NUM_TARGET = 0;
    PARAM.NUM_ANTENNA = 12;
    PARAM.NUM_EPISODE = 10^6;

    PARAM.USER = [370 400;
                  380 340; 
                  420 300; 
                  470 270; 
                  530 270; 
                  580 300; 
                  620 340; 
                  630 400];
    % PARAM.USER = [370 400];
    PARAM.UAV_START = [450 525];
    PARAM.UAV_END = [550 525];
    PARAM.UAV_Z = 100;
    PARAM.TARGET = get_target(PARAM.NUM_TARGET);
    
    PARAM.NOISE_POWER = 10^-14;
    PARAM.NOISE_POWER_SCALING = PARAM.NOISE_POWER  * PARAM.SCALING^2;

    PARAM.SENSING_TH = 5 * 10^-5;
    PARAM.SENSING_TH_SCALING = PARAM.SENSING_TH * PARAM.SCALING^2;

    PARAM.P_MAX = 0.5;
    PARAM.CHANNEL_GAIN = 10^(-6);

    PARAM.T = 15;
    PARAM.N = 3;
    PARAM.DELTA_T = PARAM.T / PARAM.N;
    PARAM.V_MAX = 1000;
    PARAM.TRUST_REGION = PARAM.DELTA_T * PARAM.V_MAX;
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
    %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
    %-----------------------------Epsiode Start-----------------------------------------------------------------------------------------------------------------------------%   
    for episode = 1 : PARAM.NUM_EPISODE

        if episode == 1
            %-----------------------------initializing precoder-----------------------------------------------------------------------------------------------------------------------------%
            [uav_t, W_t, R_t] = get_init_trajectory_SF(PARAM.TARGET, PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.NUM_TARGET, PARAM.SENSING_TH_SCALING, PARAM.P_MAX, PARAM.SCALING, PARAM.V_MAX, PARAM.N, PARAM.UAV_START(1), PARAM.UAV_END(1), PARAM.UAV_START(2), PARAM.UAV_Z);
            %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
      
            %-----------------------------get channel and steering-----------------------------------------------------------------------------------------------------------------------------%
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
            %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
            [user_rate_current, ~] = get_test_trajectory(W_t, R_t, PARAM.P_MAX, PARAM.SENSING_TH_SCALING, PARAM.NUM_TARGET, channel_t, channel_her_t, PARAM.NOISE_POWER_SCALING, steering_target_t, steering_target_her_t, distance_target_t, PARAM.N);
            user_rate_episode(:,:,1) = user_rate_current;
        end
        
        user_rate_prev = user_rate_current;

        %-----------------------------initializing precoder optimize variable-----------------------------------------------------------------------------------------------------------------------------%
        alpha = zeros(PARAM.NUM_USER, PARAM.N);
        alpha_tmp = zeros(PARAM.NUM_USER, PARAM.N);
        B = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N);

        for n = 1 : PARAM.N
        
            for k = 1 : PARAM.NUM_USER
    
                for i = 1 : PARAM.NUM_USER
                    if i == k
                        continue;
                    end
    
                    alpha_tmp(k, n) = alpha_tmp(k, n) + real(trace(channel_t(:,k,n) * channel_her_t(k,:,n) * W_t(:,:,i,n)));
                end
    
                alpha(k,n) = alpha_tmp(k,n) + real(trace(channel_t(:,k,n) * channel_her_t(k,:,n) * R_t(:,:,n))) + PARAM.NOISE_POWER_SCALING;
                alpha(k,n) = log(alpha(k,n)) / log(2);
        
                B(:,:,k,n) = channel_t(:,k,n) * channel_her_t(k,:,n);
                B(:,:,k,n) = B(:,:,k,n) / (alpha_tmp(k,n) + real(trace(channel_t(:,k,n) * channel_her_t(k,:,n) * R_t(:,:,n))) + PARAM.NOISE_POWER_SCALING);
                B(:,:,k,n) = B(:,:,k,n) / log(2);
    
            end
        end
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
        %-----------------------------precoder CVX start-----------------------------------------------------------------------------------------------------------------------------%
        cvx_begin

            cvx_solver Mosek

            expressions sum_rate(PARAM.NUM_USER, PARAM.N)
            expressions objective_1(PARAM.NUM_USER, PARAM.N)
            expressions objective_1_tmp(PARAM.NUM_USER, PARAM.N)
            expressions tmp(PARAM.NUM_USER, PARAM.N)
            expressions tmp_tmp(PARAM.NUM_USER, PARAM.N)

            expressions objective_2(PARAM.NUM_USER, PARAM.N)
            expressions objective_2_tmp(PARAM.NUM_USER, PARAM.N)

            expressions sensing_constraint(PARAM.NUM_TARGET, PARAM.N)

            variable W(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N) complex
            variable R(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.N) complex

            for n = 1 :PARAM.N

                sensing_constraint_tmp = 0;
                power_constraint_tmp = 0;

                for k = 1 : PARAM.NUM_USER
    
                    for i = 1:PARAM.NUM_USER
                        objective_1_tmp(k,n) = objective_1_tmp(k,n) + real(trace(channel_t(:,k,n) * channel_her_t(k,:,n) * W(:,:,i,n)));
    
                        if i == k
                            continue;
                        end
    
                        objective_2_tmp(k,n) = objective_2_tmp(k,n) + real(trace(B(:,:,k,n) * (W(:,:,i,n) - W_t(:,:,i,n))));
                    end
    
                    tmp(k,n) = objective_1_tmp(k,n) + real(trace(channel_t(:,k,n) * channel_her_t(k,:,n) * R(:,:,n))) + PARAM.NOISE_POWER_SCALING;
                    tmp_tmp(k,n) = -rel_entr(1, tmp(k,n));
                    objective_1(k,n) = tmp_tmp(k,n) / log(2);
        
                    objective_2(k,n) = alpha(k,n) + objective_2_tmp(k,n) + real(trace(B(:,:,k,n) * (R(:,:,n) - R_t(:,:,n))));
        
                    sum_rate(k,n) = objective_1(k,n) - objective_2(k,n);
            
                    sensing_constraint_tmp = sensing_constraint_tmp + W(:,:,k,n);
                    power_constraint_tmp = power_constraint_tmp + real(trace(W(:,:,k,n)));
                end
    
                power_constraint = power_constraint_tmp + real(trace(R(:,:,n)));
    
                for j = 1 : PARAM.NUM_TARGET
                    sensing_constraint(j,n) = real(steering_target_her_t(j,:,n) * (sensing_constraint_tmp + R(:,:,n)) * steering_target_t(:,j,n));
                end
    
                subject to
    
                    for k = 1 : PARAM.NUM_USER
                        W(:,:,k,n) == hermitian_semidefinite(PARAM.NUM_ANTENNA);
                    end
    
                    R(:,:,n) == hermitian_semidefinite(PARAM.NUM_ANTENNA);
    
                    power_constraint <= PARAM.P_MAX;
    
                    for j = 1 : PARAM.NUM_TARGET
                        sensing_constraint(j,n) >= PARAM.SENSING_TH_SCALING * distance_target_t(j,n)^2;
                    end
            end

            maximize(sum(sum(sum_rate)));

        cvx_end
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
        
        [W_opt, R_opt] = get_precoder_opt_trajectory(channel_t, channel_her_t, W, R, PARAM.NUM_USER, PARAM.NUM_ANTENNA, PARAM.N);
      
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
            get_display(distance_user_t, 'Distance user     : ');
            get_display(uav_t, 'UAV position      : ');
            %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

            if sum(user_rate_current_UAV(:, 2:PARAM.N-1)) > sum(user_rate_prev_UAV(:, 2:PARAM.N-1))
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

        if sum(sum(user_rate_current)) - sum(sum(user_rate_prev)) < 1e-6
            break;
        end

        R_t = R_opt;
        W_t = W_opt;
    end
    %-----------------------------Epsiode End-----------------------------------------------------------------------------------------------------------------------------%
    
    % get_received_BEAM_GAIN(W_opt, R_opt, num_user, num_antenna, sensing_th, scaling, distance_user_t, distance_target_t);

end