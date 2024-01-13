function V_opt = get_precoder_opt_letter(channel, channel_her, precoder, num_user, num_antenna)

    v_opt = zeros(num_antenna, num_user);
    V_opt = zeros(num_antenna, num_antenna, num_user + 1);
    V_opt_total = 0;

    for k = 1:num_user

        v_opt(:,k) = (channel_her(k,:) * precoder(:,:,k) * channel(:,k))^(-1/2) * precoder(:,:,k) * channel(:,k);
        V_opt(:,:,k) = v_opt(:,k) * transpose(conj(v_opt(:,k)));

        V_opt_total = V_opt_total + V_opt(:,:,k);
    end

    V_opt(:,:,end) = get_precoder_total_letter(num_user + 1, precoder) - V_opt_total;

end