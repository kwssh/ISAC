function channel = get_channel(x, y, scaling)

    distance = get_distance(x, y);
    sterring_vector = get_steering(distance, scaling);

    channel_gain = 10^-6;

    channel = sqrt(channel_gain / (distance^2)) * sterring_vector;
end