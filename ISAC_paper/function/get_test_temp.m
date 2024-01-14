function [rate, error] = get_test_temp (W, W_t, R, R_t, uav, target, p_max, sensing_th, num_target, num_antenna, channel, channel_her, noise_power)

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

        distance_target = get_distance(uav, target(j,:));
        sterring_target = get_steering(distance_target);
        sterring_target_her = transpose(conj(sterring_target));

        sensing_constraint = real(sterring_target_her * (sensing_constraint_tmp + R) * sterring_target);
        
        if sensing_constraint >= sensing_th * distance_target^2

        else
            error = error + 1;
        end
               
    end

end