function W = get_transmit_precoder_com(channel_user_DL, channel_user_UL, channel_target, W_old, R, V, PSI, noise_power, PEAK, DURATION, RATE_TH_DL, RATE_TH_UL, P_MAX, uav_old, P_UAV, P_0, U_TIP, P_1, C_0, V_0, G_0, scaling)
    
    N = size(channel_user_DL, 4);
    num_user = size(channel_user_DL, 3);
    num_target = size(R, 3);
    num_antenna = size(R, 1);
    num_episode = 10^6;
    objective_val = zeros(num_episode, 1);
    
    for episode = 1 : num_episode

        E_UAV = 0;
        [epsilon_DL, delta_DL, C_DL, D_DL, epsilon_UL, delta_UL, C_UL, D_UL] = deal(zeros(num_user, N));

        cvx_begin
            cvx_solver Mosek
    
            variable W(num_antenna, num_antenna, num_user, N) complex
            expressions C_DL_new(num_user, N)
            expressions D_DL_new(num_user, N)
            expressions D_UL_new(num_user, N)
            expressions objective(num_user, N)
            expressions W_sum_new(num_antenna, num_antenna, N) 
    
            for n = 1 : N
        
                W_sum_old = sum(W_old(:,:,1:num_user,n), 3);
                W_sum_new(:,:,n) = sum(W(:,:,1:num_user,n), 3);
    
                R_sum = sum(R(:,:,1:num_target,n), 3);
        
                for k = 1 : num_user
        
                    interference_user_tmp_DL = 0;
                    interference_user_tmp_DL_new = 0;
    
                    interference_target_tmp_DL = 0;
        
                    interference_user_tmp_UL = 0;
    
                    interference_target_tmp_UL = 0;
                    interference_target_tmp_UL_new = 0;
                    
                    for i = 1 : num_user
                        if i == k
                            continue
                        end
                        interference_user_tmp_DL = interference_user_tmp_DL + real(trace(channel_user_DL(:,:,k,n) * W_old(:,:,i,n)));
                        interference_user_tmp_DL_new = interference_user_tmp_DL_new + real(trace(channel_user_DL(:,:,k,n) * W(:,:,i,n)));
    
                        interference_user_tmp_UL = interference_user_tmp_UL + PEAK * real(trace(channel_user_UL(:,:,i,n) * V(:,:,k,n)));
                    end
        
                    for j = 1 : num_target
                        interference_target_tmp_DL = interference_target_tmp_DL + PSI(n) * real(trace(channel_user_DL(:,:,k,n) * R(:,:,j,n)));
                        
                        interference_target_tmp_UL = interference_target_tmp_UL + PSI(n) * real(trace(channel_target(:,:,j,n)' * V(:,:,k,n) * channel_target(:,:,j,n) * (W_sum_old + R_sum)));
                        interference_target_tmp_UL_new = interference_target_tmp_UL_new + PSI(n) * real(trace(channel_target(:,:,j,n)' * V(:,:,k,n) * channel_target(:,:,j,n) * (W_sum_new(:,:,n) + R_sum)));
                    end
    
                    delta_DL_tmp = real(trace(channel_user_DL(:,:,k,n) * W_old(:,:,k,n)));
                    delta_DL_tmp_new = real(trace(channel_user_DL(:,:,k,n) * W(:,:,k,n)));
    
                    delta_UL_tmp = PEAK * real(trace(channel_user_UL(:,:,k,n) * V(:,:,k,n)));
    
                    delta_DL(k, n) = delta_DL_tmp / (interference_user_tmp_DL + interference_target_tmp_DL + noise_power);
                    delta_UL(k, n) = delta_UL_tmp / (interference_user_tmp_UL + interference_target_tmp_UL + noise_power);
        
                    C_DL(k, n) = DURATION * (1 + delta_DL(k, n)) * delta_DL_tmp;
                    C_DL_new(k, n) = DURATION * (1 + delta_DL(k, n)) * delta_DL_tmp_new;
    
                    D_DL(k, n) = delta_DL_tmp + interference_user_tmp_DL + interference_target_tmp_DL + noise_power;
                    D_DL_new(k, n) = delta_DL_tmp_new + interference_user_tmp_DL_new + interference_target_tmp_DL + noise_power;
        
                    C_UL(k, n) = DURATION * (1 + delta_UL(k, n)) * delta_UL_tmp;
    
                    D_UL(k, n) = delta_UL_tmp + interference_user_tmp_UL + interference_target_tmp_UL + noise_power;
                    D_UL_new(k, n) = delta_UL_tmp + interference_user_tmp_UL + interference_target_tmp_UL_new + noise_power;
        
                    epsilon_DL(k, n) = sqrt(C_DL(k, n)) / D_DL(k, n);
                    epsilon_UL(k, n) = sqrt(C_UL(k, n)) / D_UL(k, n);
        
                    objective_DL = DURATION * (log(1 + delta_DL(k, n)) - delta_DL(k, n)) + 2 * epsilon_DL(k, n) * sqrt(C_DL_new(k, n)) - epsilon_DL(k, n)^2 * D_DL_new(k, n);
                    objective_UL = DURATION * (log(1 + delta_UL(k, n)) - delta_UL(k, n)) + 2 * epsilon_UL(k, n) * sqrt(C_UL(k, n)) - epsilon_UL(k, n)^2 * D_UL_new(k, n);
        
                    objective(k, n) = objective_DL + objective_UL;
    
                    subject to
                        
                        (delta_DL_tmp_new) >= (RATE_TH_DL * (interference_user_tmp_DL_new + interference_target_tmp_DL + noise_power));
                        (delta_UL_tmp) >= (RATE_TH_UL * (interference_user_tmp_UL + interference_target_tmp_UL_new + noise_power));

                        % scaling * (delta_DL_tmp_new) >= scaling * (RATE_TH_DL * (interference_user_tmp_DL_new + interference_target_tmp_DL + noise_power));
                        % scaling * (delta_UL_tmp) >= scaling * (RATE_TH_UL * (interference_user_tmp_UL + interference_target_tmp_UL_new + noise_power));
    
                        W(:,:,k,n) == hermitian_semidefinite(num_antenna);
                end
    
                real(trace(W_sum_new(:,:,n) + PSI(n) * R_sum)) <= P_MAX;
    
                E_UAV_beam_tmp = DURATION * real(trace(W_sum_new(:,:,n) + PSI(n) * R_sum));
    
                if n < N
                    v_xy = (norm([uav_old(n,1) - uav_old(n+1,1), uav_old(n,2) - uav_old(n+1,2)])) / (DURATION * 2);
                    v_z = (norm(uav_old(n,3) - uav_old(n+1,3))) / (DURATION * 2);
    
                    E_UAV_tmp = (P_0 * (1 + 3 * v_xy^2 / U_TIP^2) + P_1 * sqrt(sqrt(1 + v_xy^4 / (4 * V_0^4)) - v_xy^2 / (2 * V_0^2)) + C_0 * (v_xy)^3 + G_0 * v_z) * 2 * DURATION;
                end
    
                E_UAV = E_UAV + E_UAV_beam_tmp + E_UAV_tmp;
            end
    
            E_UAV <= P_UAV;
    
            maximize(sum(sum(objective)));
    
        cvx_end

        objective_val(episode) = sum(sum(objective));

        if episode >= 2
            if abs(objective_val(episode) - objective_val(episode - 1)) <= 0.1
                break
            end
        end

        W_old = W;
    end
end