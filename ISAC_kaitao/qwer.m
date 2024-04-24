function sum_rate_final = qwer()
    clear
    clc
    format long
    rng(123);
    
    %-----------------------------setting parameter-----------------------------------------------------------------------------------------------------------------------------%
    PARAM.SCALING = 1000;

    PARAM.NUM_USER = 1;
    PARAM.NUM_TARGET = 1;
    PARAM.NUM_ANTENNA = 16;
    PARAM.NUM_EPISODE = 100;

    PARAM.USER = [10 0];
    PARAM.UAV_START = [-10 0];
    PARAM.UAV_END = [10 0];
    PARAM.UAV_Z = 100;
    % PARAM.TARGET = [320 160; 360 120; 440 120; 480 160];
    PARAM.TARGET = [20 0];

    PARAM.NOISE_POWER = 10^(-10);
    PARAM.NOISE_POWER_SCALING = PARAM.NOISE_POWER  * PARAM.SCALING^2;

    PARAM.SENSING_TH = 10^(-4.3);
    PARAM.SENSING_TH_SCALING = PARAM.SENSING_TH * PARAM.SCALING^2;

    PARAM.RATE_TH = 0.25;
    PARAM.RATE_TH_SCALING = PARAM.RATE_TH * PARAM.SCALING^2;

    PARAM.P_MAX = 0.000001;
    PARAM.CHANNEL_GAIN = 10^(-3);
    %PARAM.GAMMA = PARAM.CHANNEL_GAIN / PARAM.NOISE_POWER;
    PARAM.GAMMA = 1;

    PARAM.T = 80;
    PARAM.T_L = 20;
    PARAM.DELTA_T = 10;
    PARAM.L = PARAM.T / PARAM.T_L;
    PARAM.N = 1;
    PARAM.N_L = PARAM.N / PARAM.L;
    PARAM.V_MAX = 30;
    PARAM.ETA = 10^-1;
    %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
    for episode = 1 : 10^6
        
        %if episode == 1
         %   [old_A_opt, old_E_opt, old_A_bar_opt, old_E_bar_opt, old_uav] = get_init(PARAM.UAV_START(1), PARAM.UAV_END(1), PARAM.UAV_START(2), PARAM.N, PARAM.N_L, PARAM.NUM_USER, PARAM.NUM_TARGET);
        %end
        
        old_uav = [0 0];
            
        distance_user = get_distance(PARAM.USER, old_uav, PARAM.UAV_Z);
        distance_target = get_distance(PARAM.TARGET, old_uav, PARAM.UAV_Z);
    
        %[new_A_bar_opt, new_E_bar_opt] = get_slack_variable(old_A_opt, old_E_opt);
        %[new_A_opt, new_E_opt] = get_period_tmp(old_A_bar_opt, old_E_bar_opt, PARAM.NUM_ANTENNA, PARAM.P_MAX, distance_user, PARAM.NUM_USER, distance_target, PARAM.NUM_TARGET, PARAM.SENSING_TH, PARAM.N, PARAM.ETA, PARAM.N_L, PARAM.L, PARAM.RATE_TH, PARAM.GAMMA);
        % [new_A_opt, new_E_opt] = get_period(old_A_bar_opt, old_E_bar_opt, PARAM.NUM_ANTENNA, PARAM.P_MAX, distance_user, PARAM.NUM_USER, distance_target, PARAM.NUM_TARGET, PARAM.SENSING_TH, PARAM.N, PARAM.ETA, PARAM.N_L, PARAM.L, PARAM.RATE_TH, PARAM.GAMMA);
    
        %new_uav = get_uav_trajectory(old_A_opt, old_E_opt, old_A_bar_opt, old_E_bar_opt, distance_user, distance_target, PARAM.NUM_USER, PARAM.NUM_TARGET, PARAM.N, PARAM.GAMMA, PARAM.P_MAX, PARAM.NUM_ANTENNA, PARAM.SENSING_TH);
        new_uav = get_uav_trajectory_BCD(distance_user, distance_target, PARAM.NUM_USER, PARAM.NUM_TARGET, PARAM.N, PARAM.GAMMA, PARAM.P_MAX, PARAM.NUM_ANTENNA, PARAM.SENSING_TH, PARAM, old_uav);
    
    
    
    
    
    end
    
end