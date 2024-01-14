function channel = get_channel(x, y, scaling, uav_z, num_antenna)

    distance = get_distance(x, y, uav_z);
    sterring_vector = get_steering(distance, scaling, uav_z, num_antenna);

    channel_gain = 10^-6;

    channel = sqrt(channel_gain / (distance^2)) * sterring_vector;
end