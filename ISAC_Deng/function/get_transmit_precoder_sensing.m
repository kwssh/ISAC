function [R, X_DL, X_UL] = get_transmit_precoder_sensing(channel_user_DL, channel_user_UL, channel_target, W, R_old, V, PSI, noise_power, PEAK, DURATION, RATE_TH_DL, RATE_TH_UL, P_MAX, uav_old, P_UAV, P_0, U_TIP, P_1, C_0, V_0, G_0, channel_target_diff, SENSING_TH, RCS, distance_target, X_DL_old, X_UL_old)
    
    N = size(channel_user_DL, 4);
    num_user = size(channel_user_DL, 3);
    num_target = size(R_old, 3);
    num_antenna = size(R_old, 1);
    num_episode = 10^6;
    objective_val = zeros(num_episode, 1);

    for episode = 1 : num_episode

        E_UAV = 0;
    
        cvx_begin
            cvx_solver Mosek
    
            variable R(num_antenna, num_antenna, num_target, N) complex
            variable X_DL(num_user, N)
            variable X_UL(num_user, N)

            expressions objective(num_user, N)
            expressions R_sum_new(num_antenna, num_antenna, N)
    
            for n = 1 : N
        
                W_sum = sum(W(:,:,1:num_user,n), 3);

                R_sum_old = sum(R_old(:,:,1:num_target,n), 3);
                R_sum_new(:,:,n) = sum(R(:,:,1:num_target,n), 3);
        
                for k = 1 : num_user
        
                    interference_user_tmp_DL = 0;

                    interference_target_tmp_DL = 0;
                    interference_target_tmp_DL_new = 0;
        
                    interference_user_tmp_UL = 0;

                    interference_target_tmp_UL = 0;
                    interference_target_tmp_UL_new = 0;
                    
                    for i = 1 : num_user
                        if i == k
                            continue
                        end
                        interference_user_tmp_DL = interference_user_tmp_DL + real(trace(channel_user_DL(:,:,k,n) * W(:,:,i,n)));
                        interference_user_tmp_UL = interference_user_tmp_UL + PEAK * real(trace(channel_user_UL(:,:,i,n) * V(:,:,k,n)));
                    end
        
                    for j = 1 : num_target
                        interference_target_tmp_DL = interference_target_tmp_DL + PSI(n) * real(trace(channel_user_DL(:,:,k,n) * R_old(:,:,j,n)));
                        interference_target_tmp_DL_new = interference_target_tmp_DL_new + PSI(n) * real(trace(channel_user_DL(:,:,k,n) * R(:,:,j,n)));

                        interference_target_tmp_UL = interference_target_tmp_UL + PSI(n) * real(trace(channel_target(:,:,j,n)' * V(:,:,k,n) * channel_target(:,:,j,n) * (W_sum + R_sum_old)));
                        interference_target_tmp_UL_new = interference_target_tmp_UL_new + PSI(n) * real(trace(channel_target(:,:,j,n)' * V(:,:,k,n) * channel_target(:,:,j,n) * (W_sum + R_sum_new(:,:,n))));
                    end

                    delta_DL_tmp = real(trace(channel_user_DL(:,:,k,n) * W(:,:,k,n)));
                    delta_UL_tmp = PEAK * real(trace(channel_user_UL(:,:,k,n) * V(:,:,k,n)));
    
                    theta_DL = (interference_user_tmp_DL + interference_target_tmp_DL + noise_power) / X_DL_old(k, n);
                    theta_UL = (interference_user_tmp_UL + interference_target_tmp_UL + noise_power) / X_UL_old(k, n);
        
                    objective(k, n) = DURATION * log(1 + X_DL(k, n)) + DURATION * log(1 + X_UL(k, n));
    
                    subject to
                        
                        delta_DL_tmp >= RATE_TH_DL * (interference_user_tmp_DL + interference_target_tmp_DL_new + noise_power);
                        delta_UL_tmp >= RATE_TH_UL * (interference_user_tmp_UL + interference_target_tmp_UL_new + noise_power);

                        delta_DL_tmp >= (interference_user_tmp_DL + interference_target_tmp_DL_new + noise_power)^2 / (2 * theta_DL) + theta_DL * X_DL_old(k,n)^2 / 2;
                        delta_UL_tmp >= (interference_user_tmp_UL + interference_target_tmp_UL_new + noise_power)^2 / (2 * theta_UL) + theta_UL * X_UL_old(k,n)^2 / 2;

                end

                for j = 1 : num_target
                    R(:,:,j,n) == hermitian_semidefinite(num_antenna);

                    PSI(n) * real(trace(channel_target_diff(:,:,j,n)' * channel_target_diff(:,:,j,n) * R(:,:,j,n))) >= PSI(n) * noise_power / (2 * SENSING_TH * (RCS / (2 * distance_target(j,n))^2));
                end

                real(trace(W_sum + PSI(n) * R_sum_new(:,:,n))) <= P_MAX;

                E_UAV_beam_tmp = DURATION * real(trace(W_sum + PSI(n) * R_sum_new(:,:,n)));

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
            if objective_val(episode) - objective_val(episode - 1) <= 0.01
                break
            end
        end

        R_old = R;
        X_DL_old = X_DL;
        X_UL_old = X_UL;
    end
end