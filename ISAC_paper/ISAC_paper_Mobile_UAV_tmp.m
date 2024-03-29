function sum_rate_final = ISAC_paper_Mobile_UAV()
    clear
    format long
    
    [save_idx, folder_name] = init_ISAC();
    
    %-----------------------------setting parameter-----------------------------------------------------------------------------------------------------------------------------%
    PARAM.SCALING = 1000;
    PARAM.SCALING_TMP = 1;

    PARAM.NUM_USER = 4;
    PARAM.NUM_TARGET = 2;
    PARAM.NUM_ANTENNA = 12;
    PARAM.NUM_EPISODE = 10^6;

    PARAM.USER = [-300 -300; -70 -400; 70 -400; 300 -300];
    PARAM.UAV_START = [-10 0];
    PARAM.UAV_END = [10 0];
    PARAM.UAV_Z = 30;
    % PARAM.TARGET = get_target(PARAM.NUM_TARGET);
    PARAM.TARGET = [-100 100; 100 100];
    
    PARAM.NOISE_POWER = 10^-14;
    PARAM.NOISE_POWER_SCALING = PARAM.NOISE_POWER  * PARAM.SCALING^2;

    PARAM.SENSING_TH = 10^(-5);
    PARAM.SENSING_TH_SCALING = PARAM.SENSING_TH * PARAM.SCALING^2;

    PARAM.P_MAX = 0.5;
    PARAM.CHANNEL_GAIN = 10^(-6);

    PARAM.T = 5;
    PARAM.N = 5;
    PARAM.DELTA_T = PARAM.T / PARAM.N;
    PARAM.V_MAX = 10000;
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
    sensing_error_episode = zeros(PARAM.NUM_EPISODE, 1);
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
       
        
            %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
        
            uav_t = get_UAV_trajectory_tmp(uav_t, W_opt, R_opt, PARAM.N, PARAM.NUM_USER, PARAM.NUM_TARGET, PARAM.NOISE_POWER, PARAM.USER, PARAM.UAV_Z, PARAM.TARGET, PARAM.SENSING_TH, PARAM.V_MAX, PARAM.DELTA_T, PARAM);
            % uav_t = [450 525; 500 400; 550 525];
              %-----------------------------display part-----------------------------------------------------------------------------------------------------------------------------%
     
            get_display(uav_t, 'UAV position      : ');
            %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

            
            %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

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

        [user_rate_current, sensing_error] = get_test_trajectory(W_opt, R_opt, PARAM.P_MAX, PARAM.SENSING_TH_SCALING, PARAM.NUM_TARGET, channel_t, channel_her_t, PARAM.NOISE_POWER_SCALING, steering_target_t, steering_target_her_t, distance_target_t, PARAM.N);
       
        user_rate_episode(:,:,episode) = user_rate_current;
        sensing_error_episode(episode) = sensing_error;

        if abs(sum(sum(user_rate_current)) - sum(sum(user_rate_prev))) < 1e-2
            break;
        end

        R_t = R_opt;
        W_t = W_opt;
    end
    %-----------------------------Epsiode End-----------------------------------------------------------------------------------------------------------------------------%
    
    if save_idx == 1

        save_path = strcat(pwd, '\figure\', folder_name);

        if ~exist(save_path, 'dir')
            mkdir(save_path);
        end

        file_txt = fopen(strcat(save_path, '\result.txt'), 'a');

        fprintf(file_txt, 'user rate : ');
        fclose(file_txt);

        writematrix(user_rate_prev, fullfile(save_path, 'result.txt'), 'Delimiter', '\t', 'WriteMode', 'append');

        file_txt = fopen(strcat(save_path, '\result.txt'), 'a');
        
        for i = 1 : episode
        
            fprintf(file_txt, 'sum rate : %f, ', sum(sum(user_rate_episode(:,:,i))));
            fprintf(file_txt, 'number of sensing error : %d\n', sensing_error_episode(episode));
        end

        fclose(file_txt);

        fig3 = plot_UAV_trajectory(uav_t, PARAM);
        saveas(fig3, strcat(save_path, '\UAV_trajectory'));
        close(fig3);

        for n = 1 : PARAM.N
            fig1 = get_received_BEAM_GAIN(W_opt(:,:,:,n), R_opt(:,:,n), PARAM.NUM_USER, PARAM.NUM_ANTENNA, PARAM.SENSING_TH_SCALING, PARAM.SCALING, distance_user_t(:, n), distance_target_t(:, n), PARAM.UAV_Z);
            saveas(fig1, strcat(save_path, '\Beampattern_gain_degree_', num2str(n)));
            close(fig1);

            fig4 = get_received_BEAM_GAIN_no_distance(W_opt(:,:,:,n), R_opt(:,:,n), PARAM.NUM_USER, PARAM.NUM_ANTENNA, PARAM.SENSING_TH_SCALING, PARAM.SCALING, distance_user_t(:, n), distance_target_t(:, n), PARAM.UAV_Z);
            saveas(fig4, strcat(save_path, '\Beampattern_gain_degree_no_', num2str(n)));
            close(fig4);
            
            fig2 = get_received_BEAM_GAIN_eleavtion(W_opt(:,:,:,n), R_opt(:,:,n), PARAM.USER, uav_t(n,:), PARAM.TARGET, PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.UAV_Z, PARAM.NUM_TARGET);
            saveas(fig2, strcat(save_path, '\Beampattern_gain_distance_', num2str(n)));
            close(fig2);
        end
    end
end