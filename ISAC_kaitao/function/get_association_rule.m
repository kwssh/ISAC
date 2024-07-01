function [A, E] = get_association_rule(PARAM, uav)

    N_L = PARAM.ISAC_DURATION;
    N = PARAM.TOTAL_TIME_SLOT;

    A = zeros(PARAM.NUM_USER, N);
    C = zeros(PARAM.NUM_TARGET, N);
    E = zeros(PARAM.NUM_USER, PARAM.NUM_TARGET, N);

    distance_target = get_distance(PARAM.TARGET, uav, PARAM.UAV_Z);
    distance_user = get_distance(PARAM.USER, uav, PARAM.UAV_Z);

    sensing_feasible = PARAM.NUM_ANTENNA * PARAM.P_MAX - distance_target.^2 * PARAM.SENSING_TH;

    for col_start = 1 : N_L : N

        col_end = min(col_start + N_L - 1, N);

        sensing_feasible_block = sensing_feasible(:, col_start : col_end);
        % distance_user_block = distance_user(:, col_start : col_end);

        sensing_feasible_block_positive_idx = sensing_feasible_block > 0 ;
        sensing_feasible_block(~sensing_feasible_block_positive_idx) = inf;

        [~, idx_min_target] = min(sensing_feasible_block, [], 2);

        idx_C = sub2ind(size(C), (1 : PARAM.NUM_TARGET)', col_start - 1 + idx_min_target);
        C(idx_C) = 1;

        % [~, idx_min_user] = min(distance_user_block, [], 1);
        % 
        % [~, idx_min_each_user] = min(distance_user_block, [], 2);
        % [~, unique_idx, ~] = unique(idx_min_each_user, 'stable');
        % duplicate_idx = setdiff(1:length(idx_min_each_user), unique_idx);
        % idx_min_each_user(duplicate_idx) = idx_min_each_user(duplicate_idx) - 1;
        % 
        % idx_min_user(idx_min_each_user) = 1:length(idx_min_each_user);
        % 
        % idx_A = sub2ind(size(A(:,1:N_L)), idx_min_user, 1:size(A(:,1:N_L), 2));
        % 
        % A(idx_A) = 1;
    end

    eye_matrix = eye(PARAM.NUM_USER);
    A(:, 1:PARAM.NUM_USER) = eye_matrix;

    for i = PARAM.NUM_USER+1:N
        col_index = mod(i-PARAM.NUM_USER, PARAM.NUM_USER);
        if col_index == 0 
            col_index = PARAM.NUM_USER;
        end
        A(:, i) = eye_matrix(:, col_index);
    end

    for n = 1 : N
        for i = 1 : PARAM.NUM_USER
            for j = 1 : PARAM.NUM_TARGET
                E(i,j,n) = A(i,n) * C(j,n);
            end
        end
    end

end