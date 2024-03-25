function [uav_init, W_opt, R_opt] = get_init_trajectory_SF(target, num_antenna, num_user, num_target, sensing_th, p_max, scaling, V_max, N, start_x, end_x, uav_y, uav_z)

    uav_init_tmp = linspace(start_x, end_x, N);
    uav_init = [uav_init_tmp' ones(N, 1) * uav_y];
    % uav_init = [-100 100; 0 0; 100 100];

    R_opt = zeros(num_antenna, num_antenna, N);

    distance_target = zeros(num_target, N);
    steering_target = zeros(num_antenna, num_target, N);
    steering_target_her = zeros(num_target, num_antenna, N);

    for n = 1:N
        for j = 1:num_target
            distance_target(j, n) = get_distance(uav_init(n,:), target(j,:), uav_z);
            steering_target(:, j, n) = get_steering(distance_target(j, n), scaling, uav_z, num_antenna);
            steering_target_her(j, :, n) = transpose(conj(steering_target(:, j, n)));
        end
    end

    cvx_begin

        cvx_solver Mosek

        variable R_init(num_antenna, num_antenna, N) complex;

        expressions power_constraint(1, N)
        expressions sensing_constraint(num_target, N)

        minimize(1)

        subject to

            for n = 1 : N

                R_init(:, :, n) == hermitian_semidefinite(num_antenna);
    
                power_constraint(:, n) = real(trace(R_init(:, :, n)));
                power_constraint(:, n) <= p_max;
    
                for j = 1:num_target
                    sensing_constraint(j, n) = real(steering_target_her(j,:,n) * (R_init(:, :, n)) * steering_target(:,j,n));
                    sensing_constraint(j, n) >= sensing_th * distance_target(j,n)^2;
                end
            end

    cvx_end

    W_opt = zeros(num_antenna, num_antenna, num_user, N);
    R_opt = R_init;
end