function uav = get_trajectory(channel_user_DL_hat, channel_user_UL_hat, channel_target_hat, W, R, V, PSI, noise_power, PEAK, DURATION, RATE_TH_DL, RATE_TH_UL, P_MAX, uav_old, P_UAV, P_0, U_TIP, P_1, C_0, V_0, G_0, distance_user_old, distance_target_old, USER, TARGET, theta_old, channel_target_diff, RCS, sensing_th)
    
    N = size(channel_user_DL_hat, 4);
    num_user = size(channel_user_DL_hat, 3);
    num_target = size(R, 3);
    num_episode = 10^6;
    objective_val = zeros(num_episode, 1);

    for episode = 1 : num_episode

        [E_DL, F_DL] = deal(zeros(num_user, N));
        E_UAV = 0;
        
        cvx_begin
            cvx_solver Mosek
    
            variable uav(N, 3)
            variable eta_user(num_user, N)
            variable eta_target(num_target, N)
            variable theta(N-1, 1)
            variable v_xy_slack(N-1, 1)
            
            expressions objectve(num_user, N)
            expressions distance_user(num_user, N)
            expressions distance_target(num_target, N)
    
            expressions gamma_tilda_DL(num_user, N)
            expressions gamma_tilda_UL(num_user, N)
            expressions gamma_tilda_UL_right_1(num_user, N)
            expressions gamma_tilda_UL_right_2(num_user, N)
            expressions gamma_low_DL(num_user, N)
            expressions gamma_low_UL(num_user, N)
            expressions gamma_low_UL_E(num_user, N)
            expressions gamma_low_UL_F(num_user, N)

            expressions M_UL_1(num_user, N)
            expressions M_UL_2(num_user, N)
    
            distance_user = get_distance_cvx(USER, uav, distance_user);
            distance_target = get_distance_cvx(TARGET, uav, distance_target);
    
            for n = 1 : N
    
                W_sum = sum(W(:,:,1:num_user,n), 3);
                R_sum = sum(R(:,:,1:num_target,n), 3);
    
                for k = 1 : num_user
                    interference_user_tmp_DL = 0;
                    interference_target_tmp_DL = 0;
        
                    interference_user_tmp_UL = 0;
                    interference_target_tmp_UL = 0;
    
                    for i = 1 : num_user
    
                        E_UL_tmp = PEAK * real(trace(channel_user_UL_hat(:,:,i,n) * V(:,:,k,n)));
                        gamma_low_UL_E(k, n) = gamma_low_UL_E(k, n) + E_UL_tmp / distance_user_old(i, n);
                        M_UL_1(k, n) = M_UL_1(k, n) + -log2(exp(1)) * E_UL_tmp / (distance_user_old(i,n)^2) * (distance_user_old(i,n) - distance_user(i,n));
    
                        if i == k
                            continue;
                        end
    
                        interference_user_tmp_DL = interference_user_tmp_DL + real(trace(channel_user_DL_hat(:,:,k,n) * W(:,:,i,n)));
                        gamma_tilda_UL_right_1(k,n) = gamma_tilda_UL_right_1(k,n) + log_sum_exp([log(E_UL_tmp) eta_user(i, n)]);
                    end
    
                    for j = 1 : num_target
                        interference_target_tmp_DL = interference_target_tmp_DL + PSI(n) * real(trace(channel_user_DL_hat(:,:,k,n) * R(:,:,j,n)));
    
                        F_UL_tmp = real(trace(channel_target_hat(:,:,j,n)' * V(:,:,k,n) * channel_target_hat(:,:,j,n) * (W_sum + R_sum)));
                        gamma_low_UL_F(k, n) = PSI(n) * F_UL_tmp / distance_target_old(j, n);
                        M_UL_2(k, n) = M_UL_2(k, n) + -log2(exp(1)) * PSI(n) * F_UL_tmp / distance_target_old(j,n)^2 * (distance_target_old(j,n) - distance_target(j,n));
                        gamma_tilda_UL_right_2(k,n) = gamma_tilda_UL_right_2(k,n) + PSI(n) * log_sum_exp([log(F_UL_tmp) eta_target(j, n)]);
                    end
    
                    E_DL(k,n) = real(trace(channel_user_DL_hat(:,:,k,n) * W(:,:,k,n)));
                    F_DL(k,n) = interference_user_tmp_DL + interference_target_tmp_DL;
    
                    first_val_DL = log2((E_DL(k,n) + F_DL(k,n)) / distance_user_old(k,n) + noise_power);
                    first_diff_DL_tmp1 = -log2(exp(1)) * (E_DL(k,n) + F_DL(k,n)) / distance_user_old(k,n)^2;
                    first_diff_DL_tmp2 = (E_DL(k,n) + F_DL(k,n)) / distance_user_old(k,n) + noise_power;
                    first_diff = first_diff_DL_tmp1 / first_diff_DL_tmp2;
    
                    gamma_low_DL(k,n) = first_val_DL + first_diff * (distance_user(k,n) - distance_user_old(k,n));
                    gamma_tilda_DL(k,n) = gamma_low_DL(k,n) - log_sum_exp([log(F_DL(k,n)) + eta_user(k,n) log(noise_power)]);
    
                    first_val_UL = log2(gamma_low_UL_E(k, n) + gamma_low_UL_F(k, n) + noise_power);
                    
                    gamma_low_UL(k,n) = first_val_UL + (M_UL_1(k,n) + M_UL_2(k,n)) / (gamma_low_UL_E(k, n) + gamma_low_UL_F(k, n) + noise_power);
                    gamma_tilda_UL(k,n) = gamma_low_UL(k,n) - log2(gamma_tilda_UL_right_1(k,n) + gamma_tilda_UL_right_2(k,n) + noise_power);
    
                    objective(k,n) = DURATION * (gamma_tilda_DL(k,n) + gamma_tilda_UL(k,n));
    
                    subject to
                        1 / exp(eta_user(k,n)) <= uav_old(n, 3)^2 + norm([uav_old(n,1) - USER(k,1), uav_old(n,2) - USER(k,2)])^2 + 2 * uav_old(n,3) * (uav(n,3) - uav_old(n,3)) + 2 * (uav_old(n,1:2) - USER(k,1:2)) .* (uav(n,1:2) - uav_old(n,1:2));
    
                        gamma_tilda_DL(k,n) >= log2(1 + RATE_TH_DL);
                        gamma_tilda_UL(k,n) >= log2(1 + RATE_TH_UL);
                end
    
                for j = 1 : num_target
                    1 / exp(eta_target(j,n)) <= uav_old(n, 3)^2 + norm([uav_old(n,1) - TARGET(j,1), uav_old(n,2) - TARGET(j,2)])^2 + 2 * uav_old(n,3) * (uav(n,3) - uav_old(n,3)) + 2 * (uav_old(n,1:2) - TARGET(j,1:2)) .* (uav(n,1:2) - uav_old(n,1:2));
                
                    RCS_abs = abs(RCS);
                    PSI(n) * (noise_power / (2 * (RCS_abs / distance_target(j,n)^2) * channel_target_diff(:,:,j,n)' * channel_target_diff(:,:,j,n) * R(:,:,j,n))) <= PSI(n) * sensing_th;

                end

                E_UAV_beam_tmp = DURATION * real(trace(W_sum + PSI(n) * R_sum(:,:,n)));

                if n < N
                    v_xy_old = (norm([uav_old(n,1) - uav_old(n+1,1), uav_old(n,2) - uav_old(n+1,2)])) / (DURATION * 2);

                    v_xy = (norm([uav(n,1) - uav(n+1,1), uav(n,2) - uav(n+1,2)])) / (DURATION * 2);
                    v_z = (norm(uav(n,3) - uav(n+1,3))) / (DURATION * 2);

                    v_xy <= V_MAX_xy;
                    v_z <= V_MAX_z;
    
                    E_UAV_tmp = (P_0 * (1 + 3 * v_xy^2 / U_TIP^2) + P_1 * theta(n, 1) + C_0 * (v_xy)^3 + G_0 * v_z) * 2 * DURATION;
                    
                    v_xy_slack(n, 1) >= v_xy;
                    pow_p(theta(n, 1), -2) <= theta_old(n, 1)^2 + v_xy_old^2 / V_0^2 + 2 * theta_old(n, 1) * (theta(n, 1) - theta_old(n, 1)) + 2 * v_xy_old / V_0^2 * (v_xy_slack(n, 1) - v_xy_old(n, 1));
                end
    
                E_UAV = E_UAV + E_UAV_beam_tmp + E_UAV_tmp;

                Z_MIN <= uav(n, 3) <= Z_MAX;

                if n == 1
                    uav(1,:) = UAV_INIT;
                    uav(N,:) = UAV_FIN;
                end
            end

            E_UAV <= P_UAV;
    
            maximize(sum(sum(objective)))

        cvx_end

        objective_val(episode) = sum(sum(objective));

        if episode >= 2
            if objective_val(episode) - objective_val(episode - 1) <= 0.01
                break
            end
        end

        uav_old = uav;
        theta_old = theta;

        distance_user_old = get_distance(USER, uav_old);
        distance_target_old = get_distance(TARGET, uav_old);
    end
end