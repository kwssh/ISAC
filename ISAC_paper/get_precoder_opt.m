function [W_opt, R_opt] = get_precoder_opt(channel, channel_her, W, R, num_user, num_antenna)

    w_opt = zeros(num_antenna, num_user);
    W_opt = zeros(num_antenna, num_antenna, num_user);
    R_opt = zeros(num_antenna, num_antenna);

    for k = 1 : num_user
        w_opt(:,k) = (channel_her(k,:) * W(:,:,k) * channel(:,k))^(-1/2) * W(:,:,k) * channel(:,k);
        W_opt(:,:,k) = w_opt(:,k) * ctranspose(w_opt(:,k));
    end

    for idx = 1 : num_user
            R_opt = R_opt + W(:,:,idx) - W_opt(:,:,idx);
    end

    R_opt = R_opt + R;


end