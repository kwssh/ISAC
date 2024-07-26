function sum_rate_final = qwer()
    clear
    format long
    rng('shuffle');
    %-----------------------------setting parameter-----------------------------------------------------------------------------------------------------------------------------%
    PARAM.SCALING = 1000;
    PARAM.NUM_ANTENNA = 16;
    PARAM.NUM_EPISODE = 10^(6);

    PARAM.USER = [250 400; 350 450; 450 450; 550 400];
    PARAM.UAV_START = [270 200];
    PARAM.UAV_END = [530 200];
    PARAM.TARGET = [320 160; 360 120; 440 120; 480 160];

    % PARAM.USER = [-100 50; 0 300; 100 50];
    % PARAM.UAV_START = [-100 0];
    % PARAM.UAV_END = [100 0];
    % PARAM.TARGET = [0 0];

    PARAM.NUM_USER = size(PARAM.USER,1);
    PARAM.NUM_TARGET = size(PARAM.TARGET,1);

    PARAM.UAV_Z = 40;

    PARAM.NOISE_POWER = 10^(-10);
    PARAM.NOISE_POWER_SCALING = PARAM.NOISE_POWER  * PARAM.SCALING^2;

    PARAM.SENSING_TH_db = -7;
    % PARAM.SENSING_TH = 10^(0.1 * PARAM.SENSING_TH_db) * 10^(-3);
    PARAM.SENSING_TH = 12 * 10^(-5);
    PARAM.SENSING_TH_SCALING = PARAM.SENSING_TH * PARAM.SCALING^2;

    PARAM.RATE_TH = 1.5;
    PARAM.RATE_TH_SCALING = PARAM.RATE_TH * PARAM.SCALING^2;

    PARAM.P_MAX = 0.1;
    PARAM.CHANNEL_GAIN = 10^(-3);
    PARAM.GAMMA = PARAM.CHANNEL_GAIN / PARAM.NOISE_POWER;

    PARAM.TOTAL_TIME = 14;                                                    % T
    PARAM.TOTAL_DURATION = 1;                                              % delta_t
    PARAM.TOTAL_TIME_SLOT = PARAM.TOTAL_TIME / PARAM.TOTAL_DURATION;          % N

    PARAM.ISAC_TIME = 7;                                                     % T_L
    PARAM.ISAC_TIME_SLOT_NUM = PARAM.TOTAL_TIME / PARAM.ISAC_TIME;            % L
    PARAM.ISAC_DURATION = PARAM.TOTAL_TIME_SLOT / PARAM.ISAC_TIME_SLOT_NUM;   % N_L

    PARAM.V_MAX = 171;
    PARAM.ETA = 10^(-1);  % 8번
    PARAM.ETA_MIN = 10^(-10);
    % PARAM.ETA = 9.536743164062501e-08;
    PARAM.Z = 0.99;
    PARAM.EPISILON_SCA = 1;
    PARAM.EPISILON_BCD = 0.01;
    %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
    objective_val_episode = zeros(PARAM.NUM_USER, PARAM.TOTAL_TIME_SLOT, PARAM.NUM_EPISODE);
    user_rate_ISAC_sum = zeros(PARAM.NUM_USER, PARAM.TOTAL_TIME_SLOT);

    for episode_association = 1 : PARAM.NUM_EPISODE
        for episode = 1 : PARAM.NUM_EPISODE
           
            if episode_association == 1 && episode == 1
               % [old_A_opt, old_E_opt, old_A_bar_opt, old_E_bar_opt, old_uav] = get_init(PARAM.UAV_START(1), PARAM.UAV_END(1), PARAM.UAV_START(2), PARAM.TOTAL_TIME_SLOT, PARAM.NUM_USER, PARAM.NUM_TARGET, PARAM.ISAC_DURATION, PARAM);
               [old_A_opt, old_E_opt, old_A_bar_opt, old_E_bar_opt, old_uav] = get_init_fix(PARAM.UAV_START(1), PARAM.UAV_END(1), PARAM.UAV_START(2), PARAM.TOTAL_TIME_SLOT, PARAM.NUM_USER, PARAM.NUM_TARGET, PARAM.ISAC_DURATION, PARAM);
                
               % old_A_opt = old_A_bar_opt;
               % old_E_opt = old_E_bar_opt;
               % 
               % [old_A_bar_opt, old_E_bar_opt] = get_slack_variable(old_A_opt, old_E_opt);
    
               % old_A_bar_opt = old_A_opt;
               % old_E_bar_opt = old_E_opt;
    
               % old_A_opt = randi([0, 1], [PARAM.NUM_USER, PARAM.TOTAL_TIME_SLOT]);
               % old_C_opt = randi([0, 1], [PARAM.NUM_TARGET, PARAM.TOTAL_TIME_SLOT]);
               % % 
               % old_A_bar_opt = rand([PARAM.NUM_USER, PARAM.TOTAL_TIME_SLOT]);
               % old_C_bar_opt = rand([PARAM.NUM_TARGET, PARAM.TOTAL_TIME_SLOT]);
               % 
               % for n = 1 : PARAM.TOTAL_TIME_SLOT
               %     for i = 1 : PARAM.NUM_USER
               %         for j = 1: PARAM.NUM_TARGET
               %              old_E_opt(i,j,n) = old_A_opt(i,n) * old_C_opt(j,n);
               %              old_E_bar_opt(i,j,n) = old_A_bar_opt(i,n) * old_C_bar_opt(j,n);
               %         end
               %     end
               % end
    
               % [old_A_bar_opt, old_E_bar_opt] = get_slack_variable(old_A_opt, old_E_opt);
            end
    
            distance_user = get_distance(PARAM.USER, old_uav, PARAM.UAV_Z);
            distance_target = get_distance(PARAM.TARGET, old_uav, PARAM.UAV_Z);
        
            if episode_association == 1369
                disp("qwewq")
            end

            [new_A_bar_opt, new_E_bar_opt] = get_slack_variable(old_A_opt, old_E_opt);
            [new_A_opt, new_E_opt] = get_association(old_A_bar_opt, old_E_bar_opt, PARAM.NUM_ANTENNA, PARAM.P_MAX, distance_user, PARAM.NUM_USER, distance_target, PARAM.NUM_TARGET, PARAM.SENSING_TH, PARAM.TOTAL_TIME_SLOT, PARAM.ETA, PARAM.GAMMA, PARAM.ISAC_DURATION, PARAM.RATE_TH);
    
            % new_A_opt = old_A_opt;
            % new_E_opt = old_E_opt;

            % new_A_opt(new_A_opt > 0.99) = 1;
            % new_A_opt(new_A_opt < 0.01) = 0;
            % 
            % new_E_opt(new_E_opt > 0.99) = 1;
            % new_E_opt(new_E_opt < 0.01) = 0;

            % new_uav = old_uav;
            [new_uav, user_rate] = get_uav_trajectory_BCD_SCA(distance_user, distance_target, PARAM.NUM_USER, PARAM.NUM_TARGET, PARAM.TOTAL_TIME_SLOT, PARAM.GAMMA, PARAM.P_MAX, PARAM.NUM_ANTENNA, PARAM.SENSING_TH, PARAM, old_uav, PARAM.V_MAX, PARAM.TOTAL_DURATION, new_A_opt, new_E_opt, PARAM.RATE_TH, PARAM.ISAC_DURATION, PARAM.EPISILON_SCA);
    
            new_distance_user = get_distance(PARAM.USER, new_uav, PARAM.UAV_Z);
            new_distance_target = get_distance(PARAM.TARGET, new_uav, PARAM.UAV_Z);
    
            % objective_val_episode(episode) = get_objective_val(new_distance_user, new_distance_target, PARAM.NUM_USER, PARAM.TOTAL_TIME_SLOT, PARAM.GAMMA, PARAM.P_MAX, PARAM.NUM_ANTENNA, PARAM.SENSING_TH, new_A_opt, new_E_opt, new_A_bar_opt, new_E_bar_opt, PARAM.ETA, PARAM.NUM_TARGET);
            objective_val_episode(:,:,episode) = get_user_rate(PARAM.GAMMA, PARAM.NUM_ANTENNA, PARAM.P_MAX, new_distance_user, PARAM.NUM_USER, new_distance_target, PARAM.NUM_TARGET, new_E_opt, PARAM.SENSING_TH, new_A_opt, user_rate_ISAC_sum, PARAM.ISAC_TIME_SLOT_NUM);

            old_A_opt = new_A_opt;
            old_E_opt = new_E_opt;
    
            old_A_bar_opt = new_A_bar_opt;
            old_E_bar_opt = new_E_bar_opt;
    
            old_uav = new_uav;

            if episode > 1
                if sum(sum(objective_val_episode(episode))) - sum(sum(objective_val_episode(episode-1))) < PARAM.EPISILON_BCD
                    break
                end
            end
        end

        if PARAM.ETA <= PARAM.ETA_MIN
            PARAM.ETA = PARAM.ETA_MIN;
        else
            PARAM.ETA = PARAM.ETA * PARAM.Z;
        end

        A_diff_1 = abs(old_A_opt .* (1 - old_A_bar_opt));
        A_diff_2 = abs(old_A_opt - old_A_bar_opt);

        E_diff_1 = abs(old_E_opt .* (1 - old_E_bar_opt));
        E_diff_2 = abs(old_E_opt - old_E_bar_opt);

        max_A_diff_1 = max(max(A_diff_1));
        max_A_diff_2 = max(max(A_diff_2));

        max_E_diff_1 = max(max(max(E_diff_1)));
        max_E_diff_2 = max(max(max(E_diff_2)));

        max_diff = max([max_A_diff_1, max_A_diff_2, max_E_diff_1, max_E_diff_2]);

        if episode_association == 8
            disp("faefea")
        end
        
        if max_diff <= 0.0001
            break
        end
    end
end