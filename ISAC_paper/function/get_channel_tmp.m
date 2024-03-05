function channel = get_channel_tmp(x, y, scaling, uav_z, num_antenna)

    distance = get_distance(x, y, uav_z);
    sterring_vector = get_steering(distance, scaling, uav_z, num_antenna);

    channel_gain = 10^-6;

    channel = sqrt(channel_gain) * sterring_vector;
end