function [W_opt, R_opt] = get_init(num_antenna, num_user, num_target, sensing_th, p_max, steering_target, steering_target_her, distance_target)

    cvx_begin

        cvx_solver Mosek

        variable R_init(num_antenna, num_antenna) complex;

        minimize(1)

        subject to

            R_init == hermitian_semidefinite(num_antenna);

            power_constraint = real(trace(R_init));
            power_constraint <= p_max;

            for j = 1:num_target
                sensing_constraint = real(steering_target_her(j,:) * (R_init) * steering_target(:,j));
                sensing_constraint >= sensing_th * distance_target(j)^2;
            end

    cvx_end

    W_opt = zeros(num_antenna, num_antenna, num_user);
    R_opt = R_init;

end