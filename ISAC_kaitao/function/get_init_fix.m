function [A, E, A_bar, E_bar, uav_init] = get_init_fix(start_x, end_x, uav_y, N, num_user, num_target, isac_duration, PARAM)

    rng(124);

    uav_init_tmp = linspace(start_x, end_x, N);
    uav_init = [uav_init_tmp' ones(N, 1) * uav_y];

    [A, A_bar, E, E_bar] = get_association_rule(PARAM, uav_init);

    % A = zeros(num_user, N) * 0.1;
    % C = zeros(num_target, N) * 0.1;
    % E = zeros(num_user, num_target, N);
    % 
    % cols = N / num_user;
    % 
    % A(1, 1:cols) = 0.9;
    % A(2, cols+1:2*cols) = 0.9;
    % A(3, 2*cols+1:3*cols) = 0.9;
    % A(4, 3*cols+1:4*cols) = 0.9;

    % eye_matrix = eye(num_user);
    % A(:, 1:num_user) = eye_matrix;
    % 
    % for i = num_user+1:N
    %     col_index = mod(i-num_user, num_user);
    %     if col_index == 0 
    %         col_index = num_user;
    %     end
    %     A(:, i) = eye_matrix(:, col_index);
    % end
    % 
    % for i = 1:isac_duration:N
    %     for j = 1 : num_target
    %         C(j,i+j-1) = 1;
    %     end
    % end
    % 
    % C(1,1) = 0;
    % C(1,4) = 1;

    % A_tmp = rand(num_user, N);
    % A = A_tmp ./ repmat(sum(A_tmp), num_user, 1);

    % A = ones(num_user, N) / num_user;
    % C = ones(num_target, N);

    % C_tmp = rand(num_target, N);
    % C = C_tmp ./ repmat(sum(C_tmp), num_target, 1);

    % for i = 1:isac_duration:N
    %     period_sum = sum(C(:, i:min(i+isac_duration-1, N)), 2);
    %     scaling_factor = 1 ./ period_sum;
    %     C(:, i:min(i+isac_duration-1, N)) = C(:, i:min(i+isac_duration-1, N)) .* scaling_factor;
    % end
    % 
    for n = 1 : N
        for i = 1 : num_user
            for j = 1: num_target
                E(i,j,n) = A(i,n) * C(j,n);
            end
        end
    end

    % [A_bar, E_bar] = get_slack_variable(A, E);
end