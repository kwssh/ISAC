function beam_gain = get_beam_gain_UAV(precoder, distance)

    num_antenna = size(precoder,1);
    uav_z = 100;

    tmp_1 = real(trace(precoder));
    tmp_2 = 0;

    for p = 1:num_antenna

        for q = p+1:num_antenna

            tmp_2 = tmp_2 + abs(precoder(p,q)) * cos(angle(precoder(p,q)) + pi * (q-p) * uav_z / distance);
        end
    end

    beam_gain = tmp_1 + 2 * tmp_2;
end