function objective_val = get_objective_val(distance_user, distance_target, num_user, N, gamma_0, p_max, num_antenna, sensing_th, A, E, A_bar, E_bar, eta)

    objective_val = 0;

    user_rate_comm = log2(1 + (gamma_0 * num_antenna * p_max ./ (distance_user.^2)));

        for k = 1 : num_user

            user_rate_ISAC_user = log2(1 + gamma_0 * (num_antenna * p_max - distance_target.^2 * sensing_th) ./ (distance_user(k,:).^2));

            for n = 1 : N

                x = [A(k,n) ; E(k,:,n)'];

                rate_diff = (user_rate_ISAC_user(:,n) - user_rate_comm(k,n)) / 2;

                rate_tmp(1,1) = user_rate_comm(k,n) / N - ((1-A_bar(k,n))^2 + 1) / (2 * eta);
                rate_tmp(2:end,1) = rate_diff;
                rate_tmp(1,2:end) = rate_diff';
                rate_tmp(2:end,2:end) = diag(-((1 - E_bar(k,:,n)).^2 + 1) / (2 * eta));
                
                user_rate_quad = x' * rate_tmp * x;
                user_rate_linear = [A_bar(k,n) E_bar(k,:,n)] * x / eta;
                user_rate_const = -sum([A_bar(k,n) E_bar(k,:,n)].^2) / (2 * eta);
                user_rate = user_rate_quad + user_rate_linear + user_rate_const;

                objective_val = objective_val + user_rate;
            end
        end
end