function sum_rate_final = qwer()
    clear
    clc
    format long
    rng(123);
    
    %-----------------------------setting parameter-----------------------------------------------------------------------------------------------------------------------------%
    PARAM.SCALING = 0.1;
    PARAM.SCALING_2 = 10^2;

    PARAM.NUM_USER = 4;
    PARAM.NUM_TARGET = 2;
    PARAM.NUM_ANTENNA = 6;
    PARAM.NUM_EPISODE = 100;

    PARAM.USER = [-100 -100; -30 -100; 30 -100; 100 -100];
    PARAM.UAV_START = [-100 0];
    PARAM.UAV_END = [100 0];
    PARAM.UAV_Z_MIN = 40;
    PARAM.UAV_Z_MAX = 100;
    % PARAM.TARGET = [320 160; 360 120; 440 120; 480 160];
    PARAM.TARGET = [-50 50; 50 50];

    PARAM.NOISE_POWER = 10^(-14);
    PARAM.NOISE_POWER_SCALING = PARAM.NOISE_POWER  * PARAM.SCALING;

    PARAM.SENSING_TH = 10^(-0.3);

    PARAM.P_MAX = 10^(0.7);
    PARAM.CHANNEL_GAIN = 10^(-8);

    PARAM.T = 25;
    PARAM.N = 3;
    PARAM.TAU = PARAM.T / PARAM.N;

    PARAM.LoS_C = 10;
    PARAM.LoS_D = 0.6;
    PARAM.LoS_K = 0.2;

    PARAM.RCS = 2 * 10^(-5) * (1 + 1i);
    % PARAM.RCS = 25 * sqrt(2) * (1 + 1i);

    PARAM.PEAK = 10^(-0.5);
    PARAM.PSI = zeros(PARAM.N, 1);
    PARAM.PSI(1:2:end) = 1;

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

            channel_n_LoS_user_DL = (randn(PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N) + 1i * randn(PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N)) / sqrt(2);
            channel_n_LoS_user_UL = (randn(PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N) + 1i * randn(PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N)) / sqrt(2);
            channel_n_LoS_target = (randn(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_TARGET, PARAM.N) + 1i * randn(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_TARGET, PARAM.N)) / sqrt(2);

            uav_old_tmp = linspace(PARAM.UAV_START(1), PARAM.UAV_END(1), PARAM.N);
            uav_old = [uav_old_tmp' ones(PARAM.N, 1) * PARAM.UAV_START(2) repmat(PARAM.UAV_Z_MIN, length(uav_old_tmp), 1)];
            v_xy_old = sqrt(sum(diff(uav_old(:, 1:2), 1, 1).^2, 2)) ./ PARAM.TAU;
            theta_old = sqrt(sqrt(1 + v_xy_old.^4 / (4 * PARAM.V_0^4)) - v_xy_old.^2 / (2 * PARAM.V_0^2));

            distance_user_old = get_distance(PARAM.USER, uav_old);
            distance_target_old = get_distance(PARAM.TARGET, uav_old);

            [channel_user_DL, channel_user_UL, channel_target, channel_target_diff, ~, ~, ~] = get_channel(PARAM.LoS_C, PARAM.LoS_D, PARAM.LoS_K, PARAM.CHANNEL_GAIN, PARAM.NUM_ANTENNA, distance_user_old, uav_old, distance_target_old, PARAM.RCS, channel_n_LoS_user_DL, channel_n_LoS_user_UL, channel_n_LoS_target, PARAM.SCALING);

            [W_old, R_old, V_old] = get_init(channel_user_DL, channel_user_UL, channel_target, channel_target_diff, PARAM.PSI, distance_target_old, PARAM.NOISE_POWER_SCALING, PARAM.SENSING_TH, PARAM.RCS, PARAM.RATE_TH_DL, PARAM.RATE_TH_UL, PARAM.PEAK, PARAM.P_MAX, PARAM.SCALING);
        
            % W_old = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N);
            % V_old = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N);
            % R_old = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_TARGET, PARAM.N);
        end

        distance_user_old = get_distance(PARAM.USER, uav_old);
        distance_target_old = get_distance(PARAM.TARGET, uav_old);
    
        [channel_user_DL, channel_user_UL, channel_target, channel_target_diff, channel_user_hat_DL, channel_user_hat_UL, channel_target_hat] = get_channel(PARAM.LoS_C, PARAM.LoS_D, PARAM.LoS_K, PARAM.CHANNEL_GAIN, PARAM.NUM_ANTENNA, distance_user_old, uav_old, distance_target_old, PARAM.RCS, channel_n_LoS_user_DL, channel_n_LoS_user_UL, channel_n_LoS_target, PARAM.SCALING);
        W_new = get_transmit_precoder_com(channel_user_DL, channel_user_UL, channel_target, W_old, R_old, V_old, PARAM.PSI, PARAM.NOISE_POWER_SCALING, PARAM.PEAK, PARAM.TAU / 2, PARAM.RATE_TH_DL, PARAM.RATE_TH_UL, PARAM.P_MAX, uav_old, PARAM.P_UAV, PARAM.P_0, PARAM.U_TIP, PARAM.P_1, PARAM.C_0, PARAM.V_0, PARAM.G_0);
        [V_new, X_DL_old, X_UL_old] = get_receive_precoder_com(channel_user_DL, channel_user_UL, channel_target, W_new, R_old, V_old, PARAM.PSI, PARAM.NOISE_POWER_SCALING, PARAM.PEAK, PARAM.TAU / 2, PARAM.RATE_TH_UL);
        R_new = get_transmit_precoder_sensing(channel_user_DL / PARAM.SCALING * PARAM.SCALING_2, channel_user_UL / PARAM.SCALING * PARAM.SCALING_2, channel_target / PARAM.SCALING * PARAM.SCALING_2, W_new, R_old, V_new, PARAM.PSI, PARAM.NOISE_POWER_SCALING / PARAM.SCALING * PARAM.SCALING_2, PARAM.PEAK, PARAM.TAU / 2, PARAM.RATE_TH_DL, PARAM.RATE_TH_UL, PARAM.P_MAX, uav_old, PARAM.P_UAV, PARAM.P_0, PARAM.U_TIP, PARAM.P_1, PARAM.C_0, PARAM.V_0, PARAM.G_0, channel_target_diff / PARAM.SCALING * PARAM.SCALING_2, PARAM.SENSING_TH, PARAM.RCS, distance_target_old, X_DL_old, X_UL_old);
        uav_new = get_trajectory(channel_user_hat_DL, channel_user_hat_UL, channel_target_hat, W_old, R_old, V_old, PARAM.PSI, PARAM.NOISE_POWER_SCALING, PARAM.PEAK, PARAM.TAU / 2, PARAM.RATE_TH_DL, PARAM.RATE_TH_UL, PARAM.P_MAX, uav_old, PARAM.P_UAV, PARAM.P_0, PARAM.U_TIP, PARAM.P_1, PARAM.C_0, PARAM.V_0, PARAM.G_0, distance_user_old, distance_target_old, PARAM.USER, PARAM.TARGET, theta_old, channel_target_diff, PARAM.RCS, PARAM.SENSING_TH);

    
    end
    
end