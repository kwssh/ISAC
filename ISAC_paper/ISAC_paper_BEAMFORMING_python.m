function [W_opt, R_opt] = ISAC_paper_BEAMFORMING_python(BEAM_PARAM)

    %-----------------------------setting parameter-----------------------------------------------------------------------------------------------------------------------------%
    num_user = BEAM_PARAM.NUM_USER;
    num_target = BEAM_PARAM.NUM_TARGET;
    num_antenna = BEAM_PARAM.NUM_ANTENNA;
    user = BEAM_PARAM.USER;
    target = BEAM_PARAM.TARGET;
    uav = BEAM_PARAM.UAV;
    uav_z = BEAM_PARAM.UAV_Z;
    v_max = BEAM_PARAM.V_MAX;
    p_max = BEAM_PARAM.P_MAX;
    channel_gain = BEAM_PARAM.CHANNEL_GAIN;
    scaling = BEAM_PARAM.SCALING;
    noise_power = BEAM_PARAM.NOISE_POWER;
    noise_power_scaling = BEAM_PARAM.NOISE_POWER_SCALING;
    sensing_threshold = BEAM_PARAM.SENSING_THRESHOLD;
    sensing_threshold_scaling = BEAM_PARAM.SENSING_THRESHOLD_SCALING;
    %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
    channel = zeros(num_antenna, num_user);
    channel_her = zeros(num_user, num_antenna);

    steering_target = zeros(num_antenna, num_target);
    steering_target_her = zeros(num_target, num_antenna);
    distance_target = zeros(num_target, 1);
    distance_user = zeros(num_user, 1);

    for k = 1:num_user
        distance_user(k) = get_distance(uav, user(k,:), uav_z);
        channel(:,k) = get_channel(uav, user(k,:), scaling, uav_z, num_antenna);
        channel_her(k,:) = transpose(conj(channel(:,k)));
    end

    for j = 1:num_target
        distance_target(j) = get_distance(uav, target(j,:), uav_z);
        steering_target(:,j) = get_steering(distance_target(j), scaling, uav_z, num_antenna);
        steering_target_her(j,:) = transpose(conj(steering_target(:,j)));
    end

    [W_t, R_t] = get_init(num_antenna, num_user, num_target, sensing_threshold_scaling, p_max, steering_target, steering_target_her, distance_target);

    [sum_rate_prev, ~] = get_test(W_t, W_t, R_t, R_t, p_max, sensing_threshold_scaling, num_target, num_antenna, channel, channel_her, noise_power_scaling, steering_target, steering_target_her, distance_target);
    
    while(1)

        alpha = zeros(num_user, 1);
        alpha_tmp = zeros(num_user, 1);
        B = zeros(num_antenna, num_antenna, num_user);
    
        for k = 1:num_user
    
            for i = 1:num_user
                if i == k
                    continue;
                end
    
                alpha_tmp(k) = alpha_tmp(k) + real(trace(channel(:,k) * channel_her(k,:) * W_t(:,:,i)));
            end
    
            alpha(k) = alpha_tmp(k) + real(trace(channel(:,k) * channel_her(k,:) * R_t)) + noise_power_scaling;
            alpha(k) = log(alpha(k)) / log(2);
    
            B(:,:,k) = channel(:,k) * channel_her(k,:);
            B(:,:,k) = B(:,:,k) / (alpha_tmp(k) + real(trace(channel(:,k) * channel_her(k,:) * R_t)) + noise_power_scaling);
            B(:,:,k) = B(:,:,k) / log(2);
    
        end
    
        %-----------------------------precoder CVX start-----------------------------------------------------------------------------------------------------------------------------%
        cvx_begin
    
            cvx_solver Mosek
    
            expressions sum_rate(num_user, 1)
            expressions objective_1(num_user, 1)
            expressions objective_1_tmp(num_user, 1)
            expressions tmp(num_user, 1)
            expressions tmp_tmp(num_user, 1)
    
            expressions objective_2(num_user, 1)
            expressions objective_2_tmp(num_user, 1)
    
            expressions sensing_constraint(num_target, 1)
    
            variable W(num_antenna, num_antenna, num_user) complex
            variable R(num_antenna, num_antenna) complex
    
            sensing_constraint_tmp = 0;
            power_constraint_tmp = 0;
    
            for k = 1:num_user
    
                for i = 1:num_user
                    objective_1_tmp(k) = objective_1_tmp(k) + real(trace(channel(:,k) * channel_her(k,:) * W(:,:,i)));
    
                    if i == k
                        continue;
                    end
    
                    objective_2_tmp(k) = objective_2_tmp(k) + real(trace(B(:,:,k) * (W(:,:,i) - W_t(:,:,i))));
                end
    
                tmp(k) = objective_1_tmp(k) + real(trace(channel(:,k) * channel_her(k,:) * R)) + noise_power_scaling;
                tmp_tmp(k) = -rel_entr(1, tmp(k));
                objective_1(k) = tmp_tmp(k) / log(2);
    
                objective_2(k) = alpha(k) + objective_2_tmp(k) + real(trace(B(:,:,k) * (R - R_t)));
    
                sum_rate(k) = objective_1(k) - objective_2(k);
    
                sensing_constraint_tmp = sensing_constraint_tmp + W(:,:,k);
                power_constraint_tmp = power_constraint_tmp + real(trace(W(:,:,k)));
            end
    
            power_constraint = power_constraint_tmp + real(trace(R));
    
            for j = 1:num_target
                sensing_constraint(j) = real(steering_target_her(j,:) * (sensing_constraint_tmp + R) * steering_target(:,j));
            end
    
            subject to
    
                for k = 1:num_user
                    W(:,:,k) == hermitian_semidefinite(num_antenna);
                end
    
                R == hermitian_semidefinite(num_antenna);
    
                power_constraint <= p_max;
    
                for j = 1 : num_target
                    sensing_constraint(j) >= sensing_threshold_scaling * distance_target(j)^2;
                end
    
            maximize(sum(sum(sum_rate)));
    
        cvx_end
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
        
        w_opt = zeros(num_antenna, num_user);
        W_opt = zeros(num_antenna, num_antenna, num_user);
    
        for idx = 1 : num_user
            w_opt(:,idx) = (channel_her(idx,:) * W(:,:,idx) * channel(:,idx))^(-1/2) * W(:,:,idx) * channel(:,idx);
            W_opt(:,:,idx) = w_opt(:,idx) * ctranspose(w_opt(:,idx));
        end
    
        R_opt = zeros(num_antenna, num_antenna);
        for idx = 1 : num_user
            R_opt = R_opt + W(:,:,idx) - W_opt(:,:,idx);
        end
        R_opt = R_opt + R;
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
        [sum_rate_current, ~] = get_test(W_opt, W_t, R_opt, R_t, p_max, sensing_threshold_scaling, num_target, num_antenna, channel, channel_her, noise_power_scaling, steering_target, steering_target_her, distance_target);
    
        if abs(sum(sum_rate_current) - sum(sum_rate_prev)) < 1e-2
            break;
        end
    
        R_t = R_opt;
        W_t = W_opt;
        sum_rate_prev = sum_rate_current;
    end

    % fig4 = get_received_BEAM_GAIN_no_distance(W_opt, R_opt, num_user, num_antenna, sensing_threshold_scaling, scaling, distance_user, distance_target, uav_z);
end