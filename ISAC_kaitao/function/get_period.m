function [A, E] = get_period(A_bar, E_bar, channel_gain, noise_power, num_antenna, p_max, distance_user, num_user, distance_target, num_target, sensing_th, N, eta, N_L, L, rate_th)

    cvx_begin

        cvx_solver Mosek

        variable A(num_user , N)
        variable E(num_user * num_target, N)

        expressions E_sum_user(num_target, N)
        expressions user_rate_ISAC_sum(num_user, N)
        expressions beam_pattern_gain

        user_rate_left_tmp = get_user_rate(channel_gain, noise_power, num_antenna, p_max, distance_user, num_user, distance_target, num_target, E, sensing_th, A, user_rate_ISAC_sum);
        user_rate_left = sum(sum(user_rate_left_tmp)) / N;

        user_rate_A_tmp = pow_pos(abs(A .* (1 - A_bar)),2) + pow_pos(abs(A - A_bar),2);
        user_rate_E_tmp = pow_pos(abs(E .* (1 - E_bar)),2) + pow_pos(abs(E - E_bar),2);

        user_rate_A = sum(sum(user_rate_A_tmp)) / (-2 * eta);
        user_rate_E = sum(sum(user_rate_E_tmp)) / (-2 * eta);

        user_rate = user_rate_left + user_rate_A + user_rate_E;

        minimize(-user_rate);

        subject to

            sum(A,1) <= 1;

            for k = 1 : num_user
                repmat(A(k, :), num_target, 1) >= E(num_target * k - num_target + 1 : num_target * k,:);

                beam_pattern_gain = beam_pattern_gain + E(num_target * k - num_target + 1 : num_target * k,:) .* (num_antenna * p_max - (distance_target .* distance_target) * sensing_th);

                E_sum_user = E_sum_user + E(num_target * k - num_target + 1 : num_target * k,:);
            end

            for l = 1 : L
                sum(E_sum_user(:,(l-1) * N_L + 1 : l * N_L), 2) == 1;

                sum(user_rate_left_tmp(:,(l-1) * N_L + 1 : l * N_L), 2) / N_L >= rate_th;
            end

            sum(E_sum_user) <= 1;

            beam_pattern_gain >= 0;

    cvx_end
end