function [V, delta_DL, delta_UL] = get_receive_precoder_com(channel_user_DL, channel_user_UL, channel_target, W, R, V_old, PSI, noise_power, PEAK, DURATION, RATE_TH_UL)
    
    N = size(channel_user_DL, 4);
    num_user = size(channel_user_DL, 3);
    num_target = size(R, 3);
    num_antenna = size(R, 1);
    num_episode = 10^6;
    objective_val = zeros(num_episode, 1);

    for episode = 1 : num_episode

        [epsilon_DL, delta_DL, C_DL, D_DL, epsilon_UL, delta_UL, C_UL, D_UL] = deal(zeros(num_user, N));
    
        cvx_begin
            cvx_solver Mosek
    
            variable V(num_antenna, num_antenna, num_user, N) complex
            expressions C_UL_new(num_user, N)
            expressions D_UL_new(num_user, N)
            expressions objective(num_user, N)
    
            for n = 1 : N
        
                W_sum = sum(W(:,:,1:num_user,n), 3);
                R_sum = sum(R(:,:,1:num_target,n), 3);
        
                for k = 1 : num_user
        
                    interference_user_tmp_DL = 0;
                    interference_target_tmp_DL = 0;
        
                    interference_user_tmp_UL = 0;
                    interference_user_tmp_UL_new = 0;
    
                    interference_target_tmp_UL = 0;
                    interference_target_tmp_UL_new = 0;
                    
                    for i = 1 : num_user
                        if i == k
                            continue
                        end
                        interference_user_tmp_DL = interference_user_tmp_DL + real(trace(channel_user_DL(:,:,k,n) * W(:,:,i,n)));
    
                        interference_user_tmp_UL = interference_user_tmp_UL + PEAK * real(trace(channel_user_UL(:,:,i,n) * V_old(:,:,k,n)));
                        interference_user_tmp_UL_new = interference_user_tmp_UL_new + PEAK * real(trace(channel_user_UL(:,:,i,n) * V(:,:,k,n)));
                    end
        
                    for j = 1 : num_target
                        interference_target_tmp_DL = interference_target_tmp_DL + PSI(n) * real(trace(channel_user_DL(:,:,k,n) * R(:,:,j,n)));
                        
                        interference_target_tmp_UL = interference_target_tmp_UL + PSI(n) * real(trace(channel_target(:,:,j,n)' * V_old(:,:,k,n) * channel_target(:,:,j,n) * (W_sum + R_sum)));
                        interference_target_tmp_UL_new = interference_target_tmp_UL_new + PSI(n) * real(trace(channel_target(:,:,j,n)' * V(:,:,k,n) * channel_target(:,:,j,n) * (W_sum + R_sum)));
                       
                    end
    
                    delta_DL_tmp = real(trace(channel_user_DL(:,:,k,n) * W(:,:,k,n)));
    
                    delta_UL_tmp = PEAK * real(trace(channel_user_UL(:,:,k,n) * V_old(:,:,k,n)));
                    delta_UL_tmp_new = PEAK * real(trace(channel_user_UL(:,:,k,n) * V(:,:,k,n)));
    
                    delta_DL(k, n) = delta_DL_tmp / (interference_user_tmp_DL + interference_target_tmp_DL + noise_power);
                    delta_UL(k, n) = delta_UL_tmp / (interference_user_tmp_UL + interference_target_tmp_UL + noise_power);
        
                    C_DL(k, n) = DURATION * (1 + delta_DL(k, n)) * delta_DL_tmp;
                    D_DL(k, n) = delta_DL_tmp + interference_user_tmp_DL + interference_target_tmp_DL + noise_power;
        
                    C_UL(k, n) = DURATION * (1 + delta_UL(k, n)) * delta_UL_tmp;
                    C_UL_new(k, n) = DURATION * (1 + delta_UL(k, n)) * delta_UL_tmp_new;
    
                    D_UL(k, n) = delta_UL_tmp + interference_user_tmp_UL + interference_target_tmp_UL + noise_power;
                    D_UL_new(k, n) = delta_UL_tmp_new + interference_user_tmp_UL_new + interference_target_tmp_UL_new + noise_power;
        
                    epsilon_DL(k, n) = sqrt(C_DL(k, n)) / D_DL(k, n);
                    epsilon_UL(k, n) = sqrt(C_UL(k, n)) / D_UL(k, n);
        
                    objective_DL = DURATION * (log(1 + delta_DL(k, n)) - delta_DL(k, n)) + 2 * epsilon_DL(k, n) * sqrt(C_DL(k, n)) - epsilon_DL(k, n)^2 * D_DL(k, n);
                    objective_UL = DURATION * (log(1 + delta_UL(k, n)) - delta_UL(k, n)) + 2 * epsilon_UL(k, n) * sqrt(C_UL_new(k, n)) - epsilon_UL(k, n)^2 * D_UL_new(k, n);
        
                    objective(k, n) = objective_DL + objective_UL;
    
                    subject to
                        
                        delta_UL_tmp_new >= RATE_TH_UL * (interference_user_tmp_UL_new + interference_target_tmp_UL_new + noise_power);
    
                        V(:,:,k,n) == hermitian_semidefinite(num_antenna);
                        % real(trace(V(:,:,k,n))) == 1;
                end
            end
    
            maximize(sum(sum(objective)));
    
        cvx_end

        objective_val(episode) = sum(sum(objective));

        if episode >= 2
            if objective_val(episode) - objective_val(episode - 1) <= 0.01
                break
            end
        end

        V_old = V;

    end
end