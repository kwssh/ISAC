function sum_rate_final = ISAC_paper_Quasi_UAV()
    clear
    format long
    rng(123);

    scaling = 1000;
    
    num_target = 0;
    % num_target = 1;

    % num_user = 8;
    num_user = 1;

    num_antenna = 12;
    noise_power = 10^-14 * scaling^2;
    p_max = 0.5;
    channel_gain = 10^(-6);

    % sensing_th = 10^(-4.3) * scaling^2;
    sensing_th = 0 * scaling^2;

    num_episdoe = 10^6;
    
    % user = [370 400;
    %         380 340; 
    %         420 300; 
    %         470 270; 
    %         530 270; 
    %         580 300; 
    %         620 340; 
    %         630 400];
    user = [370 400];
    
    uav_t = [500 525];
    % uav_t = [520 595];
    % uav_t = [470 550];
    
    target_x_1 = randi([450, 550], num_target / 2, 1);
    target_x_2 = 1000 - target_x_1;
   
    target_y = randi([590, 610], num_target / 2, 1);
    
    target_1 = [target_x_1 target_y];
    target_2 = [target_x_2 target_y];
   
    target = [target_1; target_2];

    % target = [randi([450, 550], num_target, 1) randi([590, 610], num_target, 1)];
    
    channel_t = zeros(num_antenna, num_user);
    channel_her_t = zeros(num_user, num_antenna);

    channel = zeros(num_antenna, num_user);
    channel_her = zeros(num_user, num_antenna);

    steering_target_t = zeros(num_antenna, num_target);
    steering_target_her_t = zeros(num_target, num_antenna);
    steering_target = zeros(num_antenna, num_target);
    steering_target_her = zeros(num_target, num_antenna);

    distance_target_t = zeros(num_target, 1);
    distance_user_t = zeros(num_user, 1);
    distance_target = zeros(num_target, 1);
    distance_user = zeros(num_user, 1);

    sum_rate_episode = zeros(num_user, num_episdoe);

    for episode = 1:num_episdoe

        for k = 1:num_user
            distance_user_t(k) = get_distance(uav_t, user(k,:));
            channel_t(:,k) = get_channel(uav_t, user(k,:), scaling);
            channel_her_t(k,:) = transpose(conj(channel_t(:,k)));
        end
    
        for j = 1:num_target
            distance_target_t(j) = get_distance(uav_t, target(j,:));
            steering_target_t(:,j) = get_steering(distance_target_t(j), scaling);
            steering_target_her_t(j,:) = transpose(conj(steering_target_t(:,j)));
        end

        if episode == 1
            [W_t, R_t] = get_init(num_antenna, num_user, num_target, sensing_th, p_max, steering_target_t, steering_target_her_t, distance_target_t);

            [sum_rate_current, ~, ~] = get_test(W_t, W_t, R_t, R_t, p_max, sensing_th, num_target, num_antenna, channel_t, channel_her_t, noise_power, steering_target_t, steering_target_her_t, distance_target_t);
            sum_rate_episode(:,1) = sum_rate_current;
        end

        sum_rate_prev = sum_rate_current;
    
        alpha = zeros(num_user, 1);
        alpha_tmp = zeros(num_user, 1);
        sensing_constraint_tmp = 0;
        power_constraint_tmp = 0;
        B = zeros(num_antenna, num_antenna, num_user);
    
        for k = 1:num_user
    
            for i = 1:num_user
                if i == k
                    continue;
                end
    
                alpha_tmp(k) = alpha_tmp(k) + real(trace(channel_t(:,k) * channel_her_t(k,:) * W_t(:,:,i)));
            end
    
            alpha(k) = alpha_tmp(k) + real(trace(channel_t(:,k) * channel_her_t(k,:) * R_t)) + noise_power;
            alpha(k) = log(alpha(k)) / log(2);
    
            B(:,:,k) = channel_t(:,k) * channel_her_t(k,:);
            B(:,:,k) = B(:,:,k) / (alpha_tmp(k) + real(trace(channel_t(:,k) * channel_her_t(k,:) * R_t)) + noise_power);
            B(:,:,k) = B(:,:,k) / log(2);
    
        end
    
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
        
            for k = 1:num_user
        
                for i = 1:num_user
                    objective_1_tmp(k) = objective_1_tmp(k) + real(trace(channel_t(:,k) * channel_her_t(k,:) * W(:,:,i)));
    
                    if i == k
                        continue;
                    end
    
                    objective_2_tmp(k) = objective_2_tmp(k) + real(trace(B(:,:,k) * (W(:,:,i) - W_t(:,:,i))));
                end
    
                tmp(k) = objective_1_tmp(k) + real(trace(channel_t(:,k) * channel_her_t(k,:) * R)) + noise_power;
                tmp_tmp(k) = -rel_entr(1, tmp(k));
                objective_1(k) = tmp_tmp(k) / log(2);
    
                objective_2(k) = alpha(k) + objective_2_tmp(k) + real(trace(B(:,:,k) * (R - R_t)));
    
                sum_rate(k) = objective_1(k) - objective_2(k);
        
                sensing_constraint_tmp = sensing_constraint_tmp + W(:,:,k);
                power_constraint_tmp = power_constraint_tmp + real(trace(W(:,:,k)));
            end

            power_constraint = power_constraint_tmp + real(trace(R));

            for j = 1:num_target
                sensing_constraint(j) = real(steering_target_her_t(j,:) * (sensing_constraint_tmp + R) * steering_target_t(:,j));
            end
        
            maximize(sum(sum_rate));
        
            subject to
        
                for k = 1:num_user
                    W(:,:,k) == hermitian_semidefinite(num_antenna);
                end

                R == hermitian_semidefinite(num_antenna);
                
                power_constraint <= p_max;

                for j = 1 : num_target
                    sensing_constraint(j) >= sensing_th * distance_target_t(j)^2;
                end
        
        cvx_end

        [W_opt, R_opt] = get_precoder_opt(channel_t, channel_her_t, W, R, num_user, num_antenna);

%-----------------------------------------------------------------------------------------------------------------------------%

        % get_search(W_opt, W_t, R_opt, R_t, p_max, sensing_th, num_target, num_antenna, noise_power, user, target, num_user, scaling);

        [sum_rate_prev_UAV, ~, error_prev_UAV] = get_test(W_opt, W_t, R_opt, R_t, p_max, sensing_th, num_target, num_antenna, channel_t, channel_her_t, noise_power, steering_target_t, steering_target_her_t, distance_target_t);
        
        trust_region = 0.001;

        while(1)

            uav = get_UAV_position(uav_t, W_opt, R_opt, user, num_user, channel_gain, noise_power / scaling^2, sensing_th / scaling^2, num_target, target, trust_region);
            % q = get_UAV_trajectory_2(uav_t * scaling, W_opt, R_opt, user, num_user, channel_gain, noise_power, sensing_th, num_target, target, trust_region, 10000000, 1, scaling);

            disp(episode);
            disp(trust_region);
            disp(sum(sum_rate_prev_UAV));
            
            

            for k = 1:num_user
                distance_user(k) = get_distance(uav, user(k,:));
                channel(:,k) = get_channel(uav, user(k,:), scaling);
                channel_her(k,:) = transpose(conj(channel(:,k)));
            end

            for j = 1:num_target
                distance_target(j) = get_distance(uav, target(j,:));
                steering_target(:,j) = get_steering(distance_target(j), scaling);
                steering_target_her(j,:) = transpose(conj(steering_target(:,j)));
            end

            [sum_rate_current_UAV, ~, error_current_UAV] = get_test(W_opt, W_t, R_opt, R_t, p_max, sensing_th, num_target, num_antenna, channel, channel_her, noise_power, steering_target, steering_target_her, distance_target);

            disp(sum(sum_rate_current_UAV));
            
            disp(sum(sum_rate_current_UAV) - sum(sum_rate_prev_UAV));
            disp(uav_t);
            if sum(sum_rate_current_UAV) > sum(sum_rate_prev_UAV)
                uav_t = uav;
                trust_region = 0.001;

                sum_rate_prev_UAV = sum_rate_current_UAV;
            else
                trust_region = trust_region / 2;
            end

            if trust_region < 10^(-6)
                break
            end
        end

%--------------------------------------------------------------------------------------------%
        for k = 1:num_user
            distance_user_t(k) = get_distance(uav_t, user(k,:));
            channel_t(:,k) = get_channel(uav_t, user(k,:), scaling);
            channel_her_t(k,:) = transpose(conj(channel_t(:,k)));
        end
    
        for j = 1:num_target
            distance_target_t(j) = get_distance(uav_t, target(j,:));
            steering_target_t(:,j) = get_steering(distance_target_t(j), scaling);
            steering_target_her_t(j,:) = transpose(conj(steering_target_t(:,j)));
        end

        [sum_rate_current, ~, error] = get_test(W_opt, W_t, R_opt, R_t, p_max, sensing_th, num_target, num_antenna, channel_t, channel_her_t, noise_power, steering_target_t, steering_target_her_t, distance_target_t);

        disp(error);
            
        sum_rate_episode(:,episode) = sum_rate_current;

        if sum(sum_rate_current) - sum(sum_rate_prev) < 1e-6
            break;
        end

        R_t = R_opt;
        W_t = W_opt;
    end

    get_received_BEAM_GAIN(W_opt, R_opt, num_user, num_antenna, sensing_th, scaling, distance_user_t, distance_target_t);

end