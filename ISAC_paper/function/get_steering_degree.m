function steering_vector = get_steering_degree(degree, scaling, num_antenna)

    steering_vector = zeros(num_antenna, 1);

    for i = 1:num_antenna
        steering_vector(i) = exp(1j * pi * (i-1) * cos(deg2rad(degree))) * scaling;
    end
end
