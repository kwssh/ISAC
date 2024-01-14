function beam_gain = get_beam_gain_diff_UAV(precoder, distance, uav_t, user)

    num_antenna = size(precoder,1);
    uav_z = 100;

    tmp = 0;

    for p = 1:num_antenna

        for q = p+1:num_antenna

            tmp = tmp + abs(precoder(p,q)) * sin(angle(precoder(p,q)) + pi * (q-p) * uav_z / distance) * ...
                        uav_z * (q-p) * (uav_t - user) / (2 *  distance^3);
        end
    end

    beam_gain = 4 * pi * tmp;
end