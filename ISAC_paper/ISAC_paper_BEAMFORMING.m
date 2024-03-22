function sum_rate_final = ISAC_paper_BEAMFORMING()
    clear
    format long
    rng(123);

    %-----------------------------setting parameter-----------------------------------------------------------------------------------------------------------------------------%
    PARAM.SCALING = 1000;

    PARAM.NUM_USER = 2;
    PARAM.NUM_TARGET = 0;
    PARAM.NUM_ANTENNA = 12;
    PARAM.NUM_EPISODE = 100;

    PARAM.USER = [370 400; 630 400];
    PARAM.UAV_T = [370 400];
    PARAM.UAV_Z = 10;
    PARAM.TARGET = get_target(PARAM.NUM_TARGET);
    PARAM.TARGET = [500 525];

    PARAM.NOISE_POWER = 10^-14;
    PARAM.NOISE_POWER_SCALING = PARAM.NOISE_POWER  * PARAM.SCALING^2;

    PARAM.SENSING_TH = 10^(-4);
    PARAM.SENSING_TH_SCALING = PARAM.SENSING_TH * PARAM.SCALING^2;

    PARAM.P_MAX = 0.5;
    %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
    %-----------------------------initializing variable-----------------------------------------------------------------------------------------------------------------------------%
    channel_t = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_USER);
    channel_her_t = zeros(PARAM.NUM_USER, PARAM.NUM_ANTENNA);

    steering_target_t = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_TARGET);
    steering_target_her_t = zeros(PARAM.NUM_TARGET, PARAM.NUM_ANTENNA);
    distance_target_t = zeros(PARAM.NUM_TARGET, 1);
    distance_user_t = zeros(PARAM.NUM_USER, 1);

    user_rate_episode = zeros(PARAM.NUM_USER, PARAM.NUM_EPISODE);
    %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
    for episode = 1 : PARAM.NUM_EPISODE

        %-----------------------------get channel and steering-----------------------------------------------------------------------------------------------------------------------------%
        for k = 1 : PARAM.NUM_USER
            distance_user_t(k) = get_distance(PARAM.UAV_T, PARAM.USER(k,:), PARAM.UAV_Z);
            channel_t(:,k) = get_channel(PARAM.UAV_T, PARAM.USER(k,:), PARAM.SCALING, PARAM.UAV_Z, PARAM.NUM_ANTENNA);
            channel_her_t(k,:) = transpose(conj(channel_t(:,k)));
        end
    
        for j = 1 : PARAM.NUM_TARGET
            distance_target_t(j) = get_distance(PARAM.UAV_T, PARAM.TARGET(j,:), PARAM.UAV_Z);
            steering_target_t(:,j) = get_steering(distance_target_t(j), PARAM.SCALING, PARAM.UAV_Z, PARAM.NUM_ANTENNA);
            steering_target_her_t(j,:) = transpose(conj(steering_target_t(:,j)));
        end
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
        %-----------------------------initializing precoder-----------------------------------------------------------------------------------------------------------------------------%
        if episode == 1

            [W_t, R_t] = get_init(PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.NUM_TARGET, PARAM.SENSING_TH_SCALING, PARAM.P_MAX, steering_target_t, steering_target_her_t, distance_target_t);

            [user_rate_current, ~, ~] = get_test(W_t, W_t, R_t, R_t, PARAM.P_MAX, PARAM.SENSING_TH_SCALING, PARAM.NUM_TARGET, PARAM.NUM_ANTENNA, channel_t, channel_her_t, PARAM.NOISE_POWER_SCALING, steering_target_t, steering_target_her_t, distance_target_t);
            user_rate_episode(:,1) = user_rate_current;

        end
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
        user_rate_prev = user_rate_current;

        %-----------------------------initializing precoder optimize variable-----------------------------------------------------------------------------------------------------------------------------%
        alpha = zeros(PARAM.NUM_USER, 1);
        alpha_tmp = zeros(PARAM.NUM_USER, 1);
        B = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_USER);

        sensing_constraint_tmp = 0;
        power_constraint_tmp = 0;
        
        for k = 1 : PARAM.NUM_USER

            for i = 1 : PARAM.NUM_USER
                if i == k
                    continue;
                end

                alpha_tmp(k) = alpha_tmp(k) + real(trace(channel_t(:,k) * channel_her_t(k,:) * W_t(:,:,i)));
            end

            alpha(k) = alpha_tmp(k) + real(trace(channel_t(:,k) * channel_her_t(k,:) * R_t)) + PARAM.NOISE_POWER_SCALING;
            alpha(k) = log(alpha(k)) / log(2);

            B(:,:,k) = channel_t(:,k) * channel_her_t(k,:);
            B(:,:,k) = B(:,:,k) / (alpha_tmp(k) + real(trace(channel_t(:,k) * channel_her_t(k,:) * R_t)) + PARAM.NOISE_POWER_SCALING);
            B(:,:,k) = B(:,:,k) / log(2);

        end
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
        %-----------------------------precoder CVX start-----------------------------------------------------------------------------------------------------------------------------%
        cvx_begin

            cvx_solver Mosek

            expressions sum_rate(PARAM.NUM_USER, 1)
            expressions objective_1(PARAM.NUM_USER, 1)
            expressions objective_1_tmp(PARAM.NUM_USER, 1)
            expressions tmp(PARAM.NUM_USER, 1)
            expressions tmp_tmp(PARAM.NUM_USER, 1)

            expressions objective_2(PARAM.NUM_USER, 1)
            expressions objective_2_tmp(PARAM.NUM_USER, 1)

            expressions sensing_constraint(PARAM.NUM_TARGET, 1)

            variable W(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_USER) complex
            variable R(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA) complex

            for k = 1 : PARAM.NUM_USER

                for i = 1:PARAM.NUM_USER
                    objective_1_tmp(k) = objective_1_tmp(k) + real(trace(channel_t(:,k) * channel_her_t(k,:) * W(:,:,i)));

                    if i == k
                        continue;
                    end

                    objective_2_tmp(k) = objective_2_tmp(k) + real(trace(B(:,:,k) * (W(:,:,i) - W_t(:,:,i))));
                end

                tmp(k) = objective_1_tmp(k) + real(trace(channel_t(:,k) * channel_her_t(k,:) * R)) + PARAM.NOISE_POWER_SCALING;
                tmp_tmp(k) = -rel_entr(1, tmp(k));
                objective_1(k) = tmp_tmp(k) / log(2);

                objective_2(k) = alpha(k) + objective_2_tmp(k) + real(trace(B(:,:,k) * (R - R_t)));

                sum_rate(k) = objective_1(k) - objective_2(k);

                sensing_constraint_tmp = sensing_constraint_tmp + W(:,:,k);
                power_constraint_tmp = power_constraint_tmp + real(trace(W(:,:,k)));
            end

            power_constraint = power_constraint_tmp + real(trace(R));

            for j = 1 : PARAM.NUM_TARGET
                sensing_constraint(j) = real(steering_target_her_t(j,:) * (sensing_constraint_tmp + R) * steering_target_t(:,j));
            end

            maximize(sum(sum_rate));

            subject to

                for k = 1 : PARAM.NUM_USER
                    W(:,:,k) == hermitian_semidefinite(PARAM.NUM_ANTENNA);
                end

                R == hermitian_semidefinite(PARAM.NUM_ANTENNA);

                power_constraint <= PARAM.P_MAX;

                for j = 1 : PARAM.NUM_TARGET
                    sensing_constraint(j) >= PARAM.SENSING_TH_SCALING * distance_target_t(j)^2;
                end

        cvx_end
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

        [W_opt, R_opt] = get_precoder_opt(channel_t, channel_her_t, W, R, PARAM.NUM_USER, PARAM.NUM_ANTENNA);

        [user_rate_current, ~, error] = get_test(W_opt, W_t, R_opt, R_t, PARAM.P_MAX, PARAM.SENSING_TH_SCALING, PARAM.NUM_TARGET, PARAM.NUM_ANTENNA, channel_t, channel_her_t, PARAM.NOISE_POWER_SCALING, steering_target_t, steering_target_her_t, distance_target_t);

        user_rate_episode(:,episode) = user_rate_current;

        if sum(user_rate_current) - sum(user_rate_prev) < 1e-6
            break;
        end

        R_t = R_opt;
        W_t = W_opt;
    end

    % get_received_BEAM_GAIN(W_opt, R_opt, PARAM.NUM_USER, PARAM.NUM_ANTENNA, PARAM.SENSING_TH_SCALING, PARAM.SCALING, distance_user_t, distance_target_t, PARAM.UAV_Z);
    % get_received_BEAM_GAIN_eleavtion(W_opt, R_opt, PARAM.USER, PARAM.UAV_T, PARAM.TARGET, PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.UAV_Z);
end