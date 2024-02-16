function [A, E, A_bar, E_bar, uav_init] = get_init(start_x, end_x, uav_y, N, N_L, num_user, num_target)

    uav_init_tmp = linspace(start_x, end_x, N);
    uav_init = [uav_init_tmp' ones(N, 1) * uav_y];

    % A = zeros(num_user, N);
    A = [eye(num_user) eye(num_user)];
    C = zeros(num_target, N);
    E = zeros(num_user * num_target, N);
    
    for j = 1 : num_target
        C(j, (j-1) * N_L + 1) = 1;
    end

    for k = 1 : num_user
        E(num_target * k - num_target + 1 : num_target * k,:) = A(k,:) .* C;
    end

    A_bar = A;
    E_bar = E;
end