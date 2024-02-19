function [A, E] = get_period_tmp(A_bar, E_bar, num_antenna, p_max, distance_user, num_user, distance_target, num_target, sensing_th, N, eta, N_L, L, rate_th, gamma_0)

    rate_tmp = zeros(1 + num_target);
    rate_th_tmp = zeros(1 + num_target);
    
    cvx_begin

        cvx_solver Mosek

        variable X(num_user + (num_user * num_target), N)

        expressions A(num_user , N)
        expressions E(num_user * num_target, N)

        expressions E_sum_user(num_target, N)
        expressions user_rate_ISAC_sum(num_user, N)
        expressions beam_pattern_gain
        expressions sum_rate
        expressions user_rate_th(num_user, N)

        A = X(1 : num_user, :);
        E = X(num_user + 1 : end, :);

        user_rate_th = get_user_rate(gamma_0, num_antenna, p_max, distance_user, num_user, distance_target, num_target, E, sensing_th, A, user_rate_ISAC_sum);
        
        user_rate_comm = log2(1 + (gamma_0 * num_antenna * p_max ./ (distance_user .* distance_user)));
    
        for k = 1 : num_user
            
            user_rate_ISAC_tmp = log2(1 + ((gamma_0 * num_antenna * p_max - (distance_target .* distance_target * sensing_th)) ./ (distance_user(k,:) .* distance_user(k,:))));
    
            for n = 1 : N

                x = [A(k,n) ; E(num_target * k - num_target + 1 : num_target * k,n)];

                rate_diff = (user_rate_ISAC_tmp(:,n) - user_rate_comm(k,n)) / 2;
                
                rate_tmp(1,1) = user_rate_comm(k,n) / N - ((1-A_bar(k,n))^2 + 1) / (2 * eta);
                rate_tmp(2:end,1) = rate_diff;
                rate_tmp(1,2:end) = rate_diff';
                rate_tmp(2:end,2:end) = diag(-((1 - E_bar(num_target * k - num_target + 1 : num_target * k,n)).^2 + 1) / (2 * eta));
                
                user_rate_quad = x' * rate_tmp * x;
                user_rate_linear = [A_bar(k,n) E_bar(num_target * k - num_target + 1 : num_target * k,n)'] * x / eta;
                user_rate_const = -sum([A_bar(k,n) E_bar(num_target * k - num_target + 1 : num_target * k,n)'].^2) / (2 * eta);
                user_rate = user_rate_quad + user_rate_linear + user_rate_const;

                sum_rate = sum_rate + user_rate;

                % rate_th_tmp(1,1) = user_rate_comm(k,n) / N;
                % rate_th_tmp(2:end,1) = rate_diff;
                % rate_th_tmp(1,2:end) = rate_diff';
                % 
                % [eig_vec, eig_val] = eig(rate_th_tmp);
                % eig_val(eig_val < 0) = 0;
                % rate_th_tmp_new = eig_vec * eig_val / (eig_vec);

                % user_rate_th(k, n) = x' * rate_th_tmp_new * x;
            end
        end

        minimize(-sum_rate);

        subject to

            sum(A,1) <= 1;

            for k = 1 : num_user
                repmat(A(k, :), num_target, 1) >= E(num_target * k - num_target + 1 : num_target * k,:);

                beam_pattern_gain = beam_pattern_gain + E(num_target * k - num_target + 1 : num_target * k,:) .* (num_antenna * p_max - (distance_target .* distance_target) * sensing_th);

                E_sum_user = E_sum_user + E(num_target * k - num_target + 1 : num_target * k,:);
            end

            for l = 1 : L
                sum(E_sum_user(:,(l-1) * N_L + 1 : l * N_L), 2) == 1;

                sum(user_rate_th(:,(l-1) * N_L + 1 : l * N_L), 2) / N_L >= rate_th;
            end

            sum(E_sum_user) <= 1;

            beam_pattern_gain >= 0;

    cvx_end
end