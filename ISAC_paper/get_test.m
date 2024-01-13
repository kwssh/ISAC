function [rate, sum_rate, error] = get_test(W, W_t, R, R_t, p_max, sensing_th, num_target, num_antenna, channel, channel_her, noise_power, steering_target, steering_target_her, distance_target)

    error = 0;
    power_constraint_tmp = 0;
    sensing_constraint_tmp = 0;

    num_user = size(W,3);
    rate = zeros(num_user, 1);

    for k = 1:num_user

        objective_1_tmp = 0;
        objective_2_tmp = 0;

        is_psd = all(eig(W(:,:,k)) >= 0);

        if is_psd

        else
            error = error + 1;
        end

        power_constraint_tmp = power_constraint_tmp + real(trace(W(:,:,k)));
        sensing_constraint_tmp = sensing_constraint_tmp + W(:,:,k);

        for i = 1:num_user

            objective_1_tmp = objective_1_tmp + real(trace(channel(:,k) * channel_her(k,:) * W(:,:,i)));

            if i == k
                continue;
            end

            objective_2_tmp = objective_2_tmp + real(trace(channel(:,k) * channel_her(k,:) * W(:,:,i)));
        end

        objective_1 = objective_1_tmp + real(trace(channel(:,k) * channel_her(k,:) * R)) + noise_power;
        objective_1 = log2(objective_1);

        objective_2 = objective_2_tmp + real(trace(channel(:,k) * channel_her(k,:) * R)) + noise_power;
        objective_2 = log2(objective_2);

        rate(k) = objective_1 - objective_2;
    end

    is_psd = all(eig(R) >= 0);

    if is_psd

    else
        error = error + 1;
    end

    power_constraint = power_constraint_tmp + real(trace(R));

    if power_constraint <= p_max

    else
        error = error + 1;
    end

    for j = 1:num_target

        sensing_constraint = real(steering_target_her(j,:) * (sensing_constraint_tmp + R) * steering_target(:,j));
        
        if sensing_constraint >= sensing_th * distance_target(j)^2

        else
            error = error + 1;
        end
               
    end
    

    sum_rate = zeros(num_user, 1);

    objective_1 = zeros(num_user, 1);
    objective_1_tmp = zeros(num_user, 1);
    tmp = zeros(num_user, 1);
    tmp_tmp = zeros(num_user, 1);

    objective_2 = zeros(num_user, 1);
    objective_2_tmp = zeros(num_user, 1);

    alpha = zeros(num_user, 1);
    alpha_tmp = zeros(num_user, 1);

    B = zeros(num_antenna, num_antenna, num_user);

    for k = 1:num_user
        for i = 1:num_user
            objective_1_tmp(k) = objective_1_tmp(k) + real(trace(channel(:,k) * channel_her(k,:) * W(:,:,i)));

            if i == k
                continue;
            end

            alpha_tmp(k) = alpha_tmp(k) + real(trace(channel(:,k) * channel_her(k,:) * W_t(:,:,i)));
        end

        alpha(k) = alpha_tmp(k) + real(trace(channel(:,k) * channel_her(k,:) * R_t)) + noise_power;
        alpha(k) = log(alpha(k)) / log(2);

        B(:,:,k) = channel(:,k) * channel_her(k,:);
        B(:,:,k) = B(:,:,k) / (alpha_tmp(k) + real(trace(channel(:,k) * channel_her(k,:) * R_t)) + noise_power);
        B(:,:,k) = B(:,:,k) / log(2);

        for i = 1:num_user

            if i == k
                continue;
            end

            objective_2_tmp(k) = objective_2_tmp(k) + real(trace(B(:,:,k) * (W(:,:,i) - W_t(:,:,i))));
        end

        tmp(k) = objective_1_tmp(k) + real(trace(channel(:,k) * channel_her(k,:) * R)) + noise_power;
        tmp_tmp(k) = log(tmp(k));
        objective_1(k) = tmp_tmp(k) / log(2);

        objective_2(k) = alpha(k) + objective_2_tmp(k) + real(trace(B(:,:,k) * (R - R_t)));

        sum_rate(k) = objective_1(k) - objective_2(k);
    end
end