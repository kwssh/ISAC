function [A, E, A_bar, E_bar, uav_init] = get_init(start_x, end_x, uav_y, N, num_user, num_target, isac_duration, PARAM)

    rng(124);

    uav_init_tmp = linspace(start_x, end_x, N);
    uav_init = [uav_init_tmp' ones(N, 1) * uav_y];

    distance_user = get_distance(PARAM.USER, uav_init, PARAM.UAV_Z);
    distance_target = get_distance(PARAM.TARGET, uav_init, PARAM.UAV_Z);

    % [A, E] = get_association_rule(PARAM, uav_init);
    [A, E] = get_association_CVX(PARAM.NUM_ANTENNA, PARAM.P_MAX, distance_user, PARAM.NUM_USER, distance_target, PARAM.NUM_TARGET, PARAM.SENSING_TH, PARAM.TOTAL_TIME_SLOT, PARAM.GAMMA, PARAM.ISAC_DURATION, PARAM.RATE_TH);
    [A, E] = get_association_rule(PARAM, uav_init);

    



    % E_bar = zeros(num_user, num_target, N);

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

    % A_bar = ones(num_user, N) / num_user;
    % C_bar = ones(num_target, N);
    % 
    % C_tmp = rand(num_target, N);
    % C = C_tmp ./ repmat(sum(C_tmp), num_target, 1);
    % 
    % for i = 1:isac_duration:N
    %     period_sum = sum(C_bar(:, i:min(i+isac_duration-1, N)), 2);
    %     scaling_factor = 1 ./ period_sum;
    %     C_bar(:, i:min(i+isac_duration-1, N)) = C_bar(:, i:min(i+isac_duration-1, N)) .* scaling_factor;
    % end
    % 
    % for n = 1 : N
    %     for i = 1 : num_user
    %         for j = 1: num_target
    %             E_bar(i,j,n) = A_bar(i,n) * C_bar(j,n);
    %         end
    %     end
    % end

    [A_bar, E_bar] = get_slack_variable(A, E);
end