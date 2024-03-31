function sum_rate_final = ISAC_paper_Mobile_UAV()
    clear
    format long
    
    [save_idx, folder_name] = init_ISAC();
    
    %-----------------------------setting parameter-----------------------------------------------------------------------------------------------------------------------------%
    PARAM.SCALING = 1000;
    PARAM.SCALING_TMP = 1;

    PARAM.NUM_USER = 4;
    PARAM.NUM_TARGET = 0;
    PARAM.NUM_ANTENNA = 12;
    PARAM.NUM_EPISODE = 10^6;

    PARAM.USER = [-300 -300; -70 -400; 70 -400; 300 -300];
    PARAM.UAV_START = [-10 0];
    PARAM.UAV_END = [10 0];
    PARAM.UAV_Z = 30;
    PARAM.TARGET = [-100 100; 100 100];
    
    PARAM.NOISE_POWER = 10^-14;
    PARAM.NOISE_POWER_SCALING = PARAM.NOISE_POWER  * PARAM.SCALING^2;

    PARAM.SENSING_TH = 10^(-4.6);
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
            [channel_t, channel_her_t, steering_target_t, steering_target_her_t, distance_target_t, distance_user_t] = get_channel_steering(PARAM, uav_t);
            %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
            [user_rate_current, ~] = get_test_trajectory(W_t, R_t, PARAM.P_MAX, PARAM.SENSING_TH_SCALING, PARAM.NUM_TARGET, channel_t, channel_her_t, PARAM.NOISE_POWER_SCALING, steering_target_t, steering_target_her_t, distance_target_t, PARAM.N);
            user_rate_episode(:,:,1) = user_rate_current;
        end
        
        user_rate_prev = user_rate_current;

        %-----------------------------optimize precoder-----------------------------------------------------------------------------------------------------------------------------%
        [W_opt, R_opt] = get_precoder(PARAM, channel_t, channel_her_t, W_t, R_t);
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
        %-----------------------------optimize UAV-----------------------------------------------------------------------------------------------------------------------------%
        uav_t = get_UAV_trajectory_tmp(uav_t, W_opt, R_opt, PARAM.N, PARAM.NUM_USER, PARAM.NUM_TARGET, PARAM.NOISE_POWER, PARAM.USER, PARAM.UAV_Z, PARAM.TARGET, PARAM.SENSING_TH, PARAM.V_MAX, PARAM.DELTA_T, PARAM);
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

        %-----------------------------get channel and steering-----------------------------------------------------------------------------------------------------------------------------%
        [channel_t, channel_her_t, steering_target_t, steering_target_her_t, distance_target_t, distance_user_t] = get_channel_steering(PARAM, uav_t);
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

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