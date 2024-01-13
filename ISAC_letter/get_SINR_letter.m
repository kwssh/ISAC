function SINR = get_SINR_letter(precoder, channel, channel_her, num_user, noise_power)

    sinr_left = zeros(num_user, 1);
    sinr_right = zeros(num_user, 1);

    for k = 1:num_user

        sinr_left(k) = real(channel_her(k,:) * precoder(:,:,k) * channel(:,k));

        for i = 1:num_user + 1

            if i == k
                continue
            end

            sinr_right(k) = sinr_right(k) + real(channel_her(k,:) * precoder(:,:,i) * channel(:,k));

        end

        sinr_right(k) = sinr_right(k) + noise_power;

    end

    SINR = sinr_left ./ sinr_right;
end