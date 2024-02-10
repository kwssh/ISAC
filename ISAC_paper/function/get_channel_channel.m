function channel = get_channel_channel(x, y, scaling, uav_z, num_antenna)

    distance = get_distance(x, y, uav_z);
    channel_tmp = ones(num_antenna, 1);

    channel_gain = 10^-6;

    channel = scaling * sqrt(channel_gain / (distance^2)) * channel_tmp;
end