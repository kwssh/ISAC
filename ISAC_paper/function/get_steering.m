function steering_vector = get_steering(distance, scaling, uav_z, num_antenna)

    steering_vector = zeros(num_antenna, 1);

    cos = uav_z / distance;

    for i = 1:num_antenna
        steering_vector(i) = exp(1j * pi * (i-1) * cos) * scaling;
    end
end
