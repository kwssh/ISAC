function sensing_power = get_sensing_power_const_letter(sensing_power, num_target, sterring_target, sterring_target_her, precoder_total, sensing_threshold)
   
    for m = 1:num_target

        sensing_power(m) = real(sterring_target_her(m,:) * precoder_total * sterring_target(:,m));
        sensing_power(m) = sensing_power(m) - sensing_threshold;
    end

end