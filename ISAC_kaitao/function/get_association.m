function [A, E] = get_association(A_bar, E_bar, num_antenna, p_max, distance_user, num_user, distance_target, num_target, sensing_th, N, eta, gamma_0, isac_duration, rate_th)
    
    rate_tmp = zeros(1 + num_target);

    cvx_begin
        cvx_solver Mosek

        variable A(num_user, N)
        variable E(num_user, num_target, N)

        expressions sum_rate
        expressions sensing_constraint(1, num_target, N)
        expressions user_rate_th(num_user, N)
        expressions user_rate_ISAC_sum(num_user, N)

        user_rate_th = get_user_rate(gamma_0, num_antenna, p_max, distance_user, num_user, distance_target, num_target, E, sensing_th, A, user_rate_ISAC_sum);
    
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

                sum_rate = sum_rate + user_rate;
            end
        end

        minimize(-sum_rate);

        subject to

            sum(A,1) <= 1;

            for j = 1 : num_target
                E_reshpae = reshape(E(:,j,:), [num_user, N]);
                A >= E_reshpae;

                for k = 1 : num_user
                    sensing_constraint_tmp = num_antenna * p_max - distance_target(j,:).^2 * sensing_th;
                    sensing_constraint_tmp_reshape = reshape(sensing_constraint_tmp, [1 size(sensing_constraint_tmp)]);
                    sensing_constraint(:,j,:) = sensing_constraint(:,j,:) + E(k,j,:) .* sensing_constraint_tmp_reshape;
                end
            end

            sensing_constraint >= 0;

            sum(sum(E,1), 2) <= 1;

            E_user_sum = sum(E,1);

            for i = 1:isac_duration:N
                sum(E_user_sum(:,:,i:i+isac_duration-1), 3) == 1;
                sum(user_rate_th(:,i:i+isac_duration-1),2) >= rate_th;
            end

    cvx_end

end