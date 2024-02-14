function sum_rate_final = qwer()
    clear
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
    PARAM.TARGET = [320 160; 480 160];

    PARAM.NOISE_POWER = 10^(-10);
    PARAM.NOISE_POWER_SCALING = PARAM.NOISE_POWER  * PARAM.SCALING^2;

    PARAM.SENSING_TH = 6 * 10^(-5);
    PARAM.SENSING_TH_SCALING = PARAM.SENSING_TH * PARAM.SCALING^2;

    PARAM.RATE_TH = 0.25;
    PARAM.RATE_TH_SCALING = PARAM.RATE_TH * PARAM.SCALING^2;

    PARAM.P_MAX = 0.1;
    PARAM.CHANNEL_GAIN = 10^(-3);

    PARAM.T = 80;
    PARAM.T_L = 20;
    PARAM.DELTA_T = 20;
    PARAM.L = PARAM.T / PARAM.T_L;
    PARAM.N = PARAM.T / PARAM.DELTA_T;
    PARAM.N_L = PARAM.N / PARAM.L;
    PARAM.V_MAX = 30;
    PARAM.ETA = 100;
    %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
    [A_bar_opt, E_bar_opt, uav_init] = get_init(PARAM.UAV_START(1), PARAM.UAV_END(1), PARAM.UAV_START(2), PARAM.N, PARAM.N_L, PARAM.NUM_USER, PARAM.NUM_TARGET);

    distance_user = get_distance(PARAM.USER, uav_init, PARAM.UAV_Z);
    distance_target = get_distance(PARAM.TARGET, uav_init, PARAM.UAV_Z);

    [A_opt, E_opt] = get_period(A_bar_opt, E_bar_opt, PARAM.CHANNEL_GAIN, PARAM.NOISE_POWER, PARAM.NUM_ANTENNA, PARAM.P_MAX, distance_user, PARAM.NUM_USER, distance_target, PARAM.NUM_TARGET, PARAM.SENSING_TH, PARAM.N, PARAM.ETA, PARAM.N_L, PARAM.L, PARAM.RATE_TH);


end