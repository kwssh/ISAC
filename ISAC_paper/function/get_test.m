function [rate, error] = get_test(W, W_t, R, R_t, p_max, sensing_th, num_target, num_antenna, channel, channel_her, noise_power, steering_target, steering_target_her, distance_target)

    error = 0;
    power_constraint_tmp = 0;
    sensing_constraint_tmp = 0;

    num_user = size(W,3);
    rate = zeros(num_user, 1);

    for k = 1:num_user

        objective_1_tmp = 0;
        objective_2_tmp = 0;

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

    for j = 1:num_target

        sensing_constraint = real(steering_target_her(j,:) * (sensing_constraint_tmp + R) * steering_target(:,j));
        
        if sensing_constraint >= sensing_th * distance_target(j)^2

        else
            error = error + 1;
        end
               
    end
end