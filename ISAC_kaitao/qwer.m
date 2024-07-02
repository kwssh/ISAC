function sum_rate_final = qwer()
    clear
    format long
    rng(123);
    
    %-----------------------------setting parameter-----------------------------------------------------------------------------------------------------------------------------%
    PARAM.SCALING = 1000;

    PARAM.NUM_USER = 4;
    PARAM.NUM_TARGET = 4;
    PARAM.NUM_ANTENNA = 16;
    PARAM.NUM_EPISODE = 10^(6);

    PARAM.USER = [250 400; 350 450; 450 450; 550 400];
    PARAM.UAV_START = [270 200];
    PARAM.UAV_END = [530 200];
    PARAM.UAV_Z = 40;
    PARAM.TARGET = [320 160; 360 120; 440 120; 480 160];

    PARAM.NOISE_POWER = 10^(-10);
    PARAM.NOISE_POWER_SCALING = PARAM.NOISE_POWER  * PARAM.SCALING^2;

    PARAM.SENSING_TH_db = -7;
    % PARAM.SENSING_TH = 10^(0.1 * PARAM.SENSING_TH_db) * 10^(-3);
    PARAM.SENSING_TH = 12 * 10^(-5) * 0;
    PARAM.SENSING_TH_SCALING = PARAM.SENSING_TH * PARAM.SCALING^2;

    PARAM.RATE_TH = 0.25;
    PARAM.RATE_TH_SCALING = PARAM.RATE_TH * PARAM.SCALING^2;

    PARAM.P_MAX = 0.1;
    PARAM.CHANNEL_GAIN = 10^(-3);
    PARAM.GAMMA = PARAM.CHANNEL_GAIN / PARAM.NOISE_POWER;

    PARAM.TOTAL_TIME = 40;                                                    % T
    PARAM.TOTAL_DURATION = 0.25;                                              % delta_t
    PARAM.TOTAL_TIME_SLOT = PARAM.TOTAL_TIME / PARAM.TOTAL_DURATION;          % N

    PARAM.ISAC_TIME = 40;                                                     % T_L
    PARAM.ISAC_TIME_SLOT_NUM = PARAM.TOTAL_TIME / PARAM.ISAC_TIME;            % L
    PARAM.ISAC_DURATION = PARAM.TOTAL_TIME_SLOT / PARAM.ISAC_TIME_SLOT_NUM;   % N_L

    PARAM.V_MAX = 30;
    PARAM.ETA = 10^(-4);
    PARAM.RATE_TH = 0.25;
    %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
    objective_val_episode = zeros(1, PARAM.NUM_EPISODE);

    for episode = 1 : PARAM.NUM_EPISODE
        
        if episode == 1
           [old_A_opt, old_E_opt, old_A_bar_opt, old_E_bar_opt, old_uav] = get_init(PARAM.UAV_START(1), PARAM.UAV_END(1), PARAM.UAV_START(2), PARAM.TOTAL_TIME_SLOT, PARAM.NUM_USER, PARAM.NUM_TARGET, PARAM.ISAC_DURATION, PARAM);
           % [old_A_opt, old_E_opt, old_A_bar_opt, old_E_bar_opt, old_uav] = get_init_fix(PARAM.UAV_START(1), PARAM.UAV_END(1), PARAM.UAV_START(2), PARAM.TOTAL_TIME_SLOT, PARAM.NUM_USER, PARAM.NUM_TARGET, PARAM.ISAC_DURATION, PARAM);
        
           old_A_opt = old_A_bar_opt;
           old_E_opt = old_E_bar_opt;

           [old_A_bar_opt, old_E_bar_opt] = get_slack_variable(old_A_opt, old_E_opt);
        end

        distance_user = get_distance(PARAM.USER, old_uav, PARAM.UAV_Z);
        distance_target = get_distance(PARAM.TARGET, old_uav, PARAM.UAV_Z);
    
        [new_A_bar_opt, new_E_bar_opt] = get_slack_variable(old_A_opt, old_E_opt);
        [new_A_opt, new_E_opt] = get_association(old_A_bar_opt, old_E_bar_opt, PARAM.NUM_ANTENNA, PARAM.P_MAX, distance_user, PARAM.NUM_USER, distance_target, PARAM.NUM_TARGET, PARAM.SENSING_TH, PARAM.TOTAL_TIME_SLOT, PARAM.ETA, PARAM.GAMMA, PARAM.ISAC_DURATION, PARAM.RATE_TH);

        [new_uav, user_rate] = get_uav_trajectory_BCD_SCA(distance_user, distance_target, PARAM.NUM_USER, PARAM.NUM_TARGET, PARAM.TOTAL_TIME_SLOT, PARAM.GAMMA, PARAM.P_MAX, PARAM.NUM_ANTENNA, PARAM.SENSING_TH, PARAM, old_uav, PARAM.V_MAX, PARAM.TOTAL_DURATION, old_A_opt, old_E_opt, PARAM.RATE_TH, PARAM.ISAC_DURATION);

        new_distance_user = get_distance(PARAM.USER, new_uav, PARAM.UAV_Z);
        new_distance_target = get_distance(PARAM.TARGET, new_uav, PARAM.UAV_Z);

        objective_val_episode(episode) = get_objective_val(new_distance_user, new_distance_target, PARAM.NUM_USER, PARAM.TOTAL_TIME_SLOT, PARAM.GAMMA, PARAM.P_MAX, PARAM.NUM_ANTENNA, PARAM.SENSING_TH, new_A_opt, new_E_opt, new_A_bar_opt, new_E_bar_opt, PARAM.ETA, PARAM.NUM_TARGET);
        if episode > 1
            if abs(sum(sum(objective_val_episode(episode))) - sum(sum(objective_val_episode(episode-1)))) < 0.01
                break
            end
        end

        old_A_opt = new_A_opt;
        old_E_opt = new_E_opt;

        old_A_bar_opt = new_A_bar_opt;
        old_E_bar_opt = new_E_bar_opt;

        old_uav = new_uav;
    end
end