function get_IPSAC(W, R, p_max, sensing_th, num_target, channel, channel_her, noise_power, steering_target, steering_target_her, distance_target, N)

    c_total_target = 0;
    

    cvx_begin

            cvx_solver Mosek

            variable c(N, num_target)

            [user_rate, ~] = get_test_trajectory_IPSAC(W, R, p_max, sensing_th, num_target, channel, channel_her, noise_power, steering_target, steering_target_her, distance_target, N, c);

            subject to
                
                for n = 1 : N
                    
                    c_total_n = 0;

                    for j = 1 : num_target

                        c_total_n = c_total_n + c()




end