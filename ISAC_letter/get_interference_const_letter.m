function interference_constraint = get_interference_const_letter(interference_constraint, interference, non_interference, num_user, channel, channel_her, precoder, sinr_threshold, noise_power)
    
    for k = 1:num_user

        non_interference(k) = real(channel_her(k,:) * precoder(:,:,k) * channel(:,k));

        for i = 1:num_user + 1

            if i == k
                continue
            end

            interference(k) = interference(k) + real(channel_her(k,:) * precoder(:,:,i) * channel(:,k));
        end
        
        interference_constraint(k) = non_interference(k) - sinr_threshold * (interference(k) + noise_power);

    end
end