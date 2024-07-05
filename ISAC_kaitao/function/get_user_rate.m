function user_rate = get_user_rate(gamma_0, num_antenna, p_max, distance_user, num_user, distance_target, num_target, E, sensing_th, A, user_rate_ISAC_sum, N)

    E_zero = ones(num_user, num_target, N);

    feasible_region = num_antenna * p_max - distance_target.^2 * sensing_th;
    feasible_negative = feasible_region < 0;
    [negative_row, negative_col] = find(feasible_negative);
    feasible_negative_idx = [negative_row, negative_col];

    for k = 1:size(feasible_negative_idx, 1)
        idx = feasible_negative_idx(k, :);
        if ~isempty(idx)
            E_zero(:, idx(1), idx(2)) = 0;
        end
    end

    user_rate_comm = log2(1 + (gamma_0 * num_antenna * p_max ./ (distance_user .* distance_user)));
    
    for k = 1 : num_user
        
        user_rate_ISAC_tmp = log2(1 + gamma_0 * (num_antenna * p_max - distance_target.^2 * sensing_th) ./ (distance_user(k,:).^2));

        E_tmp = reshape(E(k,:,:), [size(E(k,:,:),2) size(E(k,:,:),3)]);
        E_zero_tmp = reshape(E_zero(k,:,:), [size(E_zero(k,:,:),2) size(E_zero(k,:,:),3)]);

        user_rate_ISAC_sum(k,:) = sum(E_zero_tmp .* E_tmp .* (user_rate_ISAC_tmp - user_rate_comm(k,:)));
    end
    
    user_rate = A .* user_rate_comm + user_rate_ISAC_sum;
end