function [channel_t, channel_her_t, steering_target_t, steering_target_her_t, distance_target_t, distance_user_t] = get_channel_steering(PARAM, uav_t)

    channel_t = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_USER, PARAM.N);
    channel_her_t = zeros(PARAM.NUM_USER, PARAM.NUM_ANTENNA, PARAM.N);

    steering_target_t = zeros(PARAM.NUM_ANTENNA, PARAM.NUM_TARGET, PARAM.N);
    steering_target_her_t = zeros(PARAM.NUM_TARGET, PARAM.NUM_ANTENNA, PARAM.N);

    distance_target_t = zeros(PARAM.NUM_TARGET, PARAM.N);
    distance_user_t = zeros(PARAM.NUM_USER, PARAM.N);

    for n = 1 : PARAM.N
        for k = 1 : PARAM.NUM_USER
            distance_user_t(k, n) = get_distance(uav_t(n, :), PARAM.USER(k,:), PARAM.UAV_Z);
            channel_t(:, k, n) = get_channel(uav_t(n, :), PARAM.USER(k,:), PARAM.SCALING, PARAM.UAV_Z, PARAM.NUM_ANTENNA);
            channel_her_t(k, :, n) = transpose(conj(channel_t(:, k, n)));
        end
    
        for j = 1 : PARAM.NUM_TARGET
            distance_target_t(j, n) = get_distance(uav_t(n, :), PARAM.TARGET(j,:), PARAM.UAV_Z);
            steering_target_t(:, j, n) = get_steering(distance_target_t(j, n), PARAM.SCALING, PARAM.UAV_Z, PARAM.NUM_ANTENNA);
            steering_target_her_t(j, :, n) = transpose(conj(steering_target_t(:, j, n)));
        end
    end
end