function [rate, error] = get_test_trajectory_no_interference(W, R, p_max, sensing_th, num_target, channel, channel_her, noise_power, steering_target, steering_target_her, distance_target, N)

    error = 0;
    num_user = size(W,3);
    rate = zeros(num_user, N);

    for n = 1:N
        power_constraint_tmp = 0;
        sensing_constraint_tmp = 0;

        for k = 1:num_user
    
            objective_1_tmp = 0;
            objective_2_tmp = 0;
    
     
    
            power_constraint_tmp = power_constraint_tmp + real(trace(W(:,:,k,n)));
            sensing_constraint_tmp = sensing_constraint_tmp + W(:,:,k,n);
    
            
            objective_1 = real(trace(channel(:,k,n) * channel_her(k,:,n) * W(:,:,k,n))) + real(trace(channel(:,k,n) * channel_her(k,:,n) * R(:,:,n))) + noise_power;
            objective_1 = log2(objective_1);
    
            objective_2 = real(trace(channel(:,k,n) * channel_her(k,:,n) * R(:,:,n))) + noise_power;
            objective_2 = log2(objective_2);
    
            rate(k,n) = objective_1 - objective_2;
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