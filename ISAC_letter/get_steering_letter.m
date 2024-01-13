function steering_vector = get_steering_letter(num_antenna, degree, scaling)

    steering_vector = zeros(num_antenna, 1);

    radian = degree * pi / 180;

    for i = 1:num_antenna
        steering_vector(i) = exp(1j * pi * (i-1) * sin(radian)) * scaling / sqrt(num_antenna);
    end
end
