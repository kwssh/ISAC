function [A, E, A_bar, E_bar, uav_init] = get_init(start_x, end_x, uav_y, N, num_user, num_target, isac_duration)

    uav_init_tmp = linspace(start_x, end_x, N);
    uav_init = [uav_init_tmp' ones(N, 1) * uav_y];

    A = zeros(num_user, N);
    C = zeros(num_target, N);
    E = zeros(num_user, num_target, N);

    eye_matrix = eye(num_user);
    A(:, 1:num_user) = eye_matrix;

    for i = num_user+1:N
        col_index = mod(i-num_user, num_user);
        if col_index == 0 
            col_index = num_user;
        end
        A(:, i) = eye_matrix(:, col_index);
    end

    for i = 1:isac_duration:N
        for j = 1 : num_target
            C(j,i+j-1) = 1;
        end
    end

    C(1,1) = 0;
    C(1,4) = 1;

    for n = 1 : N
        for i = 1 : num_user
            for j = 1: num_target
                E(i,j,n) = A(i,n) * C(j,n);
            end
        end
    end
    
    A_bar = A;
    E_bar = E;
end