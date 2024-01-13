function sum_rate = get_sum_rate_const_letter(sum_rate_right, sum_rate_left, num_user, precoder, precoder_l, channel, channel_her, noise_power)

    alpha = zeros(num_user,1);

    for k = 1:num_user

        for i = 1:num_user + 1
            
            sum_rate_left(k) = sum_rate_left(k) + real(channel_her(k,:) * precoder(:,:,i) * channel(:,k));

            if i == k
                continue
            end
            
            alpha(k) = alpha(k) + real(channel_her(k,:) * precoder_l(:,:,i) * channel(:,k));

            sum_rate_right(k) = sum_rate_right(k) + real(channel_her(k,:) * (precoder(:,:,i) - precoder_l(:,:,i)) * channel(:,k));
            
        end

        alpha(k) = alpha(k) + noise_power;
    end

    sum_rate_left_tmp = log(sum_rate_left + noise_power) / log(2);
    sum_rate_right_tmp = log(alpha) / log(2) + sum_rate_right ./ (log(2) * alpha);
    rate = sum_rate_left_tmp - sum_rate_right_tmp;

    sum_rate = sum(rate);
end