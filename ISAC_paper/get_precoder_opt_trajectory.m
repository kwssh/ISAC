function [W_opt, R_opt] = get_precoder_opt_trajectory(channel, channel_her, W, R, num_user, num_antenna, N)

    w_opt = zeros(num_antenna, num_user, N);
    W_opt = zeros(num_antenna, num_antenna, num_user, N);
    R_opt = zeros(num_antenna, num_antenna, N);
    
    for n = 1:N
        for k = 1 : num_user
            w_opt(:,k,n) = (channel_her(k,:,n) * W(:,:,k,n) * channel(:,k,n))^(-1/2) * W(:,:,k,n) * channel(:,k,n);
            W_opt(:,:,k,n) = w_opt(:,k,n) * ctranspose(w_opt(:,k,n));
        end
    
        for idx = 1 : num_user
                R_opt(:,:,n) = R_opt(:,:,n) + W(:,:,idx,n) - W_opt(:,:,idx,n);
        end
    
        R_opt(:,:,n) = R_opt(:,:,n) + R(:,:,n);
    end
end