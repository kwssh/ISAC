function [W_opt, R_opt, V] = get_init(channel_user_DL, channel_user_UL, channel_target, channel_target_diff, PSI, distance_target, noise_power, sensing_th, RCS, rate_th_DL, rate_th_UL, PEAK, P_MAX, scaling)
    
    rng(123)
    num_antenna = size(channel_user_DL, 1);
    num_user = size(channel_user_DL, 3);
    num_target = size(channel_target, 3);
    N = size(channel_user_DL, 4);

    W_opt = zeros(num_antenna, num_antenna, num_user, N);
    R_opt = zeros(num_antenna, num_antenna, num_user, N);

    % W = zeros(num_antenna, num_antenna, num_user, N);
    V = zeros(num_antenna, num_antenna, num_user, N);

    cvx_begin
        cvx_solver Mosek

        variable W(num_antenna, num_antenna, num_user, N) complex
        variable R(num_antenna, num_antenna, num_target, N) complex
        % variable V(num_antenna, num_antenna, num_user, N) complex

        expressions R_sum(num_antenna, num_antenna, N)
        expressions W_sum(num_antenna, num_antenna, N)
        expressions interference_DL(num_user, N)
        expressions interference_UL(num_user, N)

        for n = 1 : N
           
            W_sum(:,:,n) = sum(W(:,:,1:num_user,n), 3);
            R_sum(:,:,n) = sum(R(:,:,1:num_target,n), 3);
            
            for k = 1 : num_user

                interference_user_tmp_DL = 0;
                interference_target_tmp_DL = 0;

                interference_user_tmp_UL = 0;
                interference_target_tmp_UL = 0;

                V_tmp_tmp = ones(num_antenna, 1);
                V_tmp = V_tmp_tmp ./ abs(V_tmp_tmp);
                V(:,:,k,n) = (V_tmp * V_tmp') / num_antenna;

                for i = 1 : num_user
                    if i == k
                        continue
                    end
                    interference_user_tmp_DL = interference_user_tmp_DL + real(trace(channel_user_DL(:,:,k,n) * W(:,:,i,n)));
                    interference_user_tmp_UL = interference_user_tmp_UL + PEAK * real(trace(channel_user_UL(:,:,i,n) * V(:,:,k,n)));
                end

                for j = 1 : num_target
                    alpha_sensing = abs(RCS / (2 * distance_target(j,n)));
                    PSI(n) * real(trace(channel_target_diff(:,:,j,n)' * channel_target_diff(:,:,j,n) * R(:,:,j,n))) >= PSI(n) * (noise_power / scaling) / (2 * sensing_th * alpha_sensing^2);
    
                    interference_target_tmp_DL = interference_target_tmp_DL + PSI(n) * real(trace(channel_user_DL(:,:,k,n) * R(:,:,j,n)));
                    interference_target_tmp_UL = interference_target_tmp_UL + PSI(n) * real(trace(channel_target(:,:,j,n)' * V(:,:,k,n) * channel_target(:,:,j,n) * (W_sum(:,:,n) + R_sum(:,:,n))));

                    R(:,:,j,n) == hermitian_semidefinite(num_antenna);
                end

                interference_DL(k,n) = interference_user_tmp_DL + interference_target_tmp_DL + noise_power;
                interference_UL(k,n) = interference_user_tmp_UL + interference_target_tmp_UL + noise_power;

                % scaling * (real(trace(channel_user_DL(:,:,k,n) * W(:,:,k,n)))) >= scaling * (rate_th_DL * interference_DL(k,n));
                % scaling * (PEAK * real(trace(channel_user_UL(:,:,k,n) * V(:,:,k,n)))) >= scaling * (rate_th_UL * interference_UL(k,n));

                real(trace(channel_user_DL(:,:,k,n) * W(:,:,k,n))) >= rate_th_DL * interference_DL(k,n);
                PEAK * real(trace(channel_user_UL(:,:,k,n) * V(:,:,k,n))) >= rate_th_UL * interference_UL(k,n);

                W(:,:,k,n) == hermitian_semidefinite(num_antenna);
            end

            real(trace(W_sum(:,:,n))) + PSI(n) * real(trace(R_sum(:,:,n))) <= P_MAX;
        end

        maximize(1)

    cvx_end

    for n = 1 : N
        for k = 1 : num_user
            [V_W, D_W] = eig(W(:,:,k,n));
            W_opt(:,:,k,n) = D_W(num_antenna, num_antenna) * V_W(:, num_antenna) * V_W(:, num_antenna)';

            [V_R, D_R] = eig(R(:,:,k,n));
            R_opt(:,:,k,n) = D_R(num_antenna, num_antenna) * V_R(:, num_antenna) * V_R(:, num_antenna)';
        end
    end
end