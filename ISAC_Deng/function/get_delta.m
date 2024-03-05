function delta = get_delta(channel_user, channel_target, W, R, V, PSI, noise_power, PEAK, ZETA)
    
    N = size(channel_user, 4);
    num_user = size(channel_user, 3);
    num_target = size(R, 3);

    delta_DL = zeros(num_user, N);
    delta_UL = zeros(num_user, N);
    
    for n = 1 : N
        for k = 1 : num_user

            interference_user_tmp_DL = 0;
            interference_target_tmp_DL = 0;

            interference_user_tmp_UL = 0;
            interference_target_tmp_UL = 0;
            
            for i = 1 : num_user
                if i == k
                    continue
                end
                interference_user_tmp_DL = interference_user_tmp_DL + real(trace(channel_user(:,:,k,n) * W(:,:,i,n)));
                interference_user_tmp_UL = interference_user_tmp_UL + PEAK * real(trace(channel_user(:,:,i,n) * V(:,:,k,n)));
            end

            for j = 1 : num_target
                interference_target_tmp_DL = interference_target_tmp_DL + PSI(n) * real(trace(channel_user(:,:,k,n) * R(:,:,j,n)));
            end

            delta_DL_tmp = real(trace(channel_user(:,:,k,n) * W(:,:,k,n)));
            delta_UL_tmp = real(trace(channel_user(:,:,k,n) * V(:,:,k,n)));

            delta_DL(k, n) = delta_DL_tmp / (interference_user_tmp_DL + interference_target_tmp_DL);

        end

    end
end