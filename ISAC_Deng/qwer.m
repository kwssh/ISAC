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

    PARAM.SENSING_TH = 10^(-0.3);
    PARAM.SENSING_TH_SCALING = PARAM.SENSING_TH * PARAM.SCALING^2;

    PARAM.RATE_TH = 0.25;
    PARAM.RATE_TH_SCALING = PARAM.RATE_TH * PARAM.SCALING^2;

    PARAM.P_MAX = 10^(0.7);
    PARAM.CHANNEL_GAIN = 10^(-3);
    PARAM.GAMMA = PARAM.CHANNEL_GAIN / PARAM.NOISE_POWER;

    PARAM.T = 25;
    PARAM.N = 25;
    PARAM.TAU = PARAM.T / PARAM.N;

    PARAM.LoS_C = 10;
    PARAM.LoS_D = 0.6;
    PARAM.LoS_K = 0.2;

    PARAM.RCS = 25* sqrt(2) * (1 + 1i);
    PARAM.PEAK = 10^(-0.5);
    PARAM.PSI = zeros(PARAM.N);
    PARAM.ZETA = zeros(PARAM.N);

    PARAM.RATE_TH_DL = 10^(0.03);
    PARAM.RATE_TH_UL = 10^(0.03);

    PARAM.P_UAV = 20000;
    PARAM.P_0 = 80;
    PARAM.U_TIP = 120;
    PARAM.P_1 = 31.43;
    PARAM.C_0 = 0.0046;
    PARAM.V_0 = 4;
    PARAM.G_0 = 10;
    %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
    for episode = 1 : 10^6
        
        if episode == 1

            uav_old = [100 100 150; 110 110 160; 120 120 170];
            W_old = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N);
            R_old = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_TARGET, PARAM.N);
            V_old = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N);

            X_DL_old = ones(PARAM.NUM_USER, PARAM.N);
            X_UL_old = ones(PARAM.NUM_USER, PARAM.N);
        end

        distance_user_old = get_distance(PARAM.USER, uav_old);
        distance_target_old = get_distance(PARAM.TARGET, uav_old);
    
        [channel_user_DL, channel_user_UL, channel_target, channel_target_diff] = get_channel(PARAM.LoS_C, PARAM.LoS_D, PARAM.LoS_K, PARAM.CHANNEL_GAIN, PARAM.NUM_ANTENNA, distance_user_old, uav_old, distance_target_old, PARAM.RCS);
        W_new = get_transmit_precoder_com(channel_user_DL, channel_user_UL, channel_target, W_old, R_old, V_old, PARAM.PSI, PARAM.NOISE_POWER, PARAM.PEAK, PARAM.TAU / 2, PARAM.RATE_TH_DL, PARAM.RATE_TH_UL, PARAM.P_MAX, uav_old, PARAM.P_UAV, PARAM.P_0, PARAM.U_TIP, PARAM.P_1, PARAM.C_0, PARAM.V_0, PARAM.G_0);
        V_new = get_receive_precoder_com(channel_user_DL, channel_user_UL, channel_target, W_new, R_old, V_old, PARAM.PSI, PARAM.NOISE_POWER, PARAM.PEAK, PARAM.TAU / 2, PARAM.RATE_TH_UL);
        [R_new, X_DL_new, X_UL_new] = get_transmit_precoder_sensing(channel_user_DL, channel_user_UL, channel_target, W_new, R_old, V_new, PARAM.PSI, PARAM.NOISE_POWER, PARAM.PEAK, PARAM.TAU / 2, PARAM.RATE_TH_DL, PARAM.RATE_TH_UL, PARAM.P_MAX, uav_old, PARAM.P_UAV, PARAM.P_0, PARAM.U_TIP, PARAM.P_1, PARAM.C_0, PARAM.V_0, PARAM.G_0, channel_target_diff, PARAM.SENSING_TH, PARAM.RCS, distance_target_old, X_DL_old, X_UL_old);
    
    
    end
    
end