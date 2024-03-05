function sum_rate_final = qwer()
    clear
    clc
    format long
    rng(123);
    
    %-----------------------------setting parameter-----------------------------------------------------------------------------------------------------------------------------%
    PARAM.SCALING = 1000;

    PARAM.NUM_USER = 4;
    PARAM.NUM_TARGET = 2;
    PARAM.NUM_ANTENNA = 16;
    PARAM.NUM_EPISODE = 100;

    PARAM.USER = [250 400; 350 450; 450 450; 550 400];
    PARAM.UAV_START = [280 200];
    PARAM.UAV_END = [520 200];
    PARAM.UAV_Z = 40;
    % PARAM.TARGET = [320 160; 360 120; 440 120; 480 160];
    PARAM.TARGET = [370 160; 430 160];

    PARAM.NOISE_POWER = 10^(-10);
    PARAM.NOISE_POWER_SCALING = PARAM.NOISE_POWER  * PARAM.SCALING^2;

    PARAM.SENSING_TH = 6 * 10^(-5);
    PARAM.SENSING_TH_SCALING = PARAM.SENSING_TH * PARAM.SCALING^2;

    PARAM.RATE_TH = 0.25;
    PARAM.RATE_TH_SCALING = PARAM.RATE_TH * PARAM.SCALING^2;

    PARAM.P_MAX = 0.1;
    PARAM.CHANNEL_GAIN = 10^(-3);
    PARAM.GAMMA = PARAM.CHANNEL_GAIN / PARAM.NOISE_POWER;

    PARAM.T = 80;
    PARAM.T_L = 20;
    PARAM.DELTA_T = 10;
    PARAM.L = PARAM.T / PARAM.T_L;
    PARAM.N = PARAM.T / PARAM.DELTA_T;
    PARAM.N_L = PARAM.N / PARAM.L;
    PARAM.V_MAX = 30;
    PARAM.ETA = 10^-1;

    PARAM.LoS_C = 10;
    PARAM.LoS_D = 0.6;
    PARAM.LoS_K = 0.2;

    PARAM.RCS = 25* sqrt(2) * (1 + 1i);
    PARAM.PEAK = 10^(-0.5);
    PARAM.PSI = zeros(PARAM.N);
    PARAM.ZETA = zeros(PARAM.N);
    %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
    for episode = 1 : 10^6
        
        if episode == 1

            uav_old = [100 100 150; 110 110 160; 120 120 170];
            W_old = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N);
            R_old = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_TARGET, PARAM.N);
            V_old = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N);

        end

        distance_user_old = get_distance(PARAM.USER, uav_old);
        distance_target_old = get_distance(PARAM.TARGET, uav_old);
    
        [channel_user, channel_target] = get_channel(PARAM.LoS_C, PARAM.LoS_D, PARAM.LoS_K, PARAM.CHANNEL_GAIN, PARAM.NUM_ANTENNA, distance_user_old, uav_old, distance_target_old, PARAM.RCS);
        
    
    
    
    end
    
end