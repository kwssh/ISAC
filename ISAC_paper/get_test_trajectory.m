function [rate, error] = get_test_trajectory(W, R, p_max, sensing_th, num_target, channel, channel_her, noise_power, steering_target, steering_target_her, distance_target, N)

    error = 0;
    num_user = size(W,3);
    rate = zeros(num_user, N);

    for n = 1:N
        power_constraint_tmp = 0;
        sensing_constraint_tmp = 0;

        for k = 1:num_user
    
            objective_1_tmp = 0;
            objective_2_tmp = 0;
    
            is_psd = all(eig(W(:,:,k,n)) >= 0);
    
            if is_psd
    
            else
                error = error + 1;
            end
    
            power_constraint_tmp = power_constraint_tmp + real(trace(W(:,:,k,n)));
            sensing_constraint_tmp = sensing_constraint_tmp + W(:,:,k,n);
    
            for i = 1:num_user
    
                objective_1_tmp = objective_1_tmp + real(trace(channel(:,k,n) * channel_her(k,:,n) * W(:,:,i,n)));
    
                if i == k
                    continue;
                end
    
                objective_2_tmp = objective_2_tmp + real(trace(channel(:,k,n) * channel_her(k,:,n) * W(:,:,i,n)));
            end
    
            objective_1 = objective_1_tmp + real(trace(channel(:,k,n) * channel_her(k,:,n) * R(:,:,n))) + noise_power;
            objective_1 = log2(objective_1);
    
            objective_2 = objective_2_tmp + real(trace(channel(:,k,n) * channel_her(k,:,n) * R(:,:,n))) + noise_power;
            objective_2 = log2(objective_2);
    
            rate(k,n) = objective_1 - objective_2;
        end


        is_psd = all(eig(R(:,:,n)) >= 0);
    
        if is_psd
    
        else
            error = error + 1;
        end
    
        power_constraint = power_constraint_tmp + real(trace(R(:,:,n)));
    
        if power_constraint <= p_max
    
        else
            error = error + 1;
        end
    
        for j = 1:num_target
    
            sensing_constraint = real(steering_target_her(j,:,n) * (sensing_constraint_tmp + R(:,:,n)) * steering_target(:,j,n));
            
            if sensing_constraint >= sensing_th * distance_target(j,n)^2
    
            else
                error = error + 1;
            end
                   
        end

    end
end