function [A, E] = get_association_CVX(num_antenna, p_max, distance_user, num_user, distance_target, num_target, sensing_th, N, gamma_0, isac_duration, rate_th)

    cvx_begin
            cvx_solver Mosek
    
            variable A(num_user, N)
            variable E(num_user, num_target, N)
    
            expressions sum_rate(num_user, N)
            expressions sensing_constraint(1, num_target, N)
            expressions user_rate_th(num_user, N)
            expressions user_rate_ISAC_sum(num_user, N)
            expressions user_rate_quad(num_user, N)
            expressions user_rate_linear(num_user, N)
            expressions x(1+num_target, num_user, N)

            user_rate_th = get_user_rate(gamma_0, num_antenna, p_max, distance_user, num_user, distance_target, num_target, E, sensing_th, A, user_rate_ISAC_sum, N);
    
            minimize(1);
    
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
                    sum(user_rate_th(:,i:i+isac_duration-1),2) / isac_duration >= rate_th;
                end
    
                E >= 0;
                A >= 0;
    
        cvx_end

        

end