function [uav_init, W_opt, R_opt] = get_init_trajectory(target, num_antenna, num_user, num_target, sensing_th, p_max, scaling, V_max, N, start_position, end_position)

    uav_set = [];

    for x = 450 : 25 : 550
        for y = 450 : 25 : 550

            distance_target = zeros(num_target);
            steering_target = zeros(num_antenna, num_target);
            steering_target_her = zeros(num_target, num_antenna);
            
            for j = 1:num_target
                distance_target(j) = get_distance([x y], target(j,:));
                steering_target(:, j) = get_steering(distance_target(j), scaling);
                steering_target_her(j, :) = transpose(conj(steering_target(:, j)));
            end

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

        if strcmp(cvx_status, 'Solved')
            uav_set = [uav_set; [x y]];
        end

        end
    end

    set_size = size(uav_set, 1);
    adjList = cell(set_size, 1);
    path = {};

    for i = 1:set_size

        for j = i+1:set_size

            dist = sqrt((uav_set(i,1) - uav_set(j,1))^2 + (uav_set(i,2) - uav_set(j,2))^2);
            
            if dist < V_max
                adjList{i} = [adjList{i}, j];
                adjList{j} = [adjList{j}, i];
            end
        end
    end

    visited = false(1, set_size);

    start_idx = find(ismember(uav_set, start_position, 'rows'));
    end_idx = find(ismember(uav_set, end_position, 'rows'));

    path = depth_first_search(start_idx, visited, adjList, [], path, end_idx);

    uav_init = uav_set(path{1}, :);
    uav_init = [repmat(uav_init(1,:), N - length(path{1}), 1); uav_init];

    %-----------------------------------------------------------------------------------------------------------------------------------%
    R_opt = zeros(num_antenna, num_antenna, N);

    for n = 1:N

        distance_target = zeros(num_target);
        steering_target = zeros(num_antenna, num_target);
        steering_target_her = zeros(num_target, num_antenna);
            
        for j = 1:num_target
            distance_target(j) = get_distance(uav_init(n,:), target(j,:));
            steering_target(:, j) = get_steering(distance_target(j), scaling);
            steering_target_her(j, :) = transpose(conj(steering_target(:, j)));
        end

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

        R_opt(:,:,n) = R_init;
    end

    W_opt = zeros(num_antenna, num_antenna, num_user, N);
    %-----------------------------------------------------------------------------------------------------------------------------------%

end