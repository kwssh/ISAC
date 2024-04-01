function [W_opt, R_opt] = get_precoder_no_interference(PARAM, channel_t, channel_her_t, steering_target_t, steering_target_her_t, distance_target_t)

    %-----------------------------precoder CVX start-----------------------------------------------------------------------------------------------------------------------------%
    cvx_begin

        cvx_solver Mosek

        expressions sum_rate(PARAM.NUM_USER, PARAM.N)
        expressions objective_1(PARAM.NUM_USER, PARAM.N)
        expressions objective_1_tmp(PARAM.NUM_USER, PARAM.N)
        expressions tmp(PARAM.NUM_USER, PARAM.N)
        expressions tmp_tmp(PARAM.NUM_USER, PARAM.N)

        expressions objective_2(PARAM.NUM_USER, PARAM.N)
        expressions objective_2_tmp(PARAM.NUM_USER, PARAM.N)

        expressions sensing_constraint(PARAM.NUM_TARGET, PARAM.N)

        variable W(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N) complex
        variable R(PARAM.NUM_ANTENNA, PARAM.NUM_ANTENNA, PARAM.N) complex

        for n = 1 :PARAM.N

            sensing_constraint_tmp = 0;
            power_constraint_tmp = 0;

            for k = 1 : PARAM.NUM_USER

                tmp = 1 + (real(trace(channel_t(:,k,n) * channel_her_t(k,:,n) * W(:,:,k,n))) / PARAM.NOISE_POWER_SCALING);
                sum_rate(k,n) = -rel_entr(1, tmp) / log(2);

                sensing_constraint_tmp = sensing_constraint_tmp + W(:,:,k,n);
                power_constraint_tmp = power_constraint_tmp + real(trace(W(:,:,k,n)));
            end

            power_constraint = power_constraint_tmp + real(trace(R(:,:,n)));

            for j = 1 : PARAM.NUM_TARGET
                sensing_constraint(j,n) = real(steering_target_her_t(j,:,n) * (sensing_constraint_tmp + R(:,:,n)) * steering_target_t(:,j,n));
            end

            subject to

                for k = 1 : PARAM.NUM_USER
                    W(:,:,k,n) == hermitian_semidefinite(PARAM.NUM_ANTENNA);
                end

                R(:,:,n) == hermitian_semidefinite(PARAM.NUM_ANTENNA);

                power_constraint <= PARAM.P_MAX;

                for j = 1 : PARAM.NUM_TARGET
                    1000 * (sensing_constraint(j,n)) >= 1000 * (PARAM.SENSING_TH_SCALING * distance_target_t(j,n)^2);
                end
        end

        maximize(sum(sum(sum_rate)));

    cvx_end
    %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    
    [W_opt, R_opt] = get_precoder_opt_trajectory(channel_t, channel_her_t, W, R, PARAM.NUM_USER, PARAM.NUM_ANTENNA, PARAM.N);
end