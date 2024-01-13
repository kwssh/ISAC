function [V_l, lambda_l] = get_init_letter(num_user, num_antenna, num_target, channel, channel_her, sinr_threshold, noise_power, sterring_target, sterring_target_her, sensing_threshold, p_max, p_circuit, pho, scaling)
    

    cvx_begin

        cvx_solver Mosek
    
        variable V(num_antenna, num_antenna, num_user + 1) complex;

        expressions interference(num_user, 1)
        expressions non_interference(num_user, 1)
        expressions sensing_power(num_target, 1)
        expressions interference_constraint(num_user, 1)
    
        minimize(1)
    
        subject to

            V_total = get_precoder_total_letter(num_user + 1, V);
            interference_const = get_interference_const_letter(interference_constraint, interference, non_interference, num_user, channel, channel_her, V, sinr_threshold, noise_power);
            sensing_power_const = get_sensing_power_const_letter(sensing_power, num_target, sterring_target, sterring_target_her, V_total, sensing_threshold);
            
            interference_const >= 0;
            real(trace(V_total)) <= p_max;
            sensing_power_const >= 0;

            for k = 1:num_user + 1
                    V(:,:,k) == hermitian_semidefinite(num_antenna);
            end

    cvx_end

   V_opt = get_precoder_opt_letter(channel, channel_her, V, num_user, num_antenna);

   V_opt_total = get_precoder_total_letter(num_user + 1, V_opt);

   u = real(trace(V_opt_total)) / pho + p_circuit;
   sinr = get_SINR_letter(V_opt, channel, channel_her, num_user, noise_power / scaling^2);
   t = sum(log2(1+sinr));

   lambda_l = u / t;
   V_l = V_opt;

end