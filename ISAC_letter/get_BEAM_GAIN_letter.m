function get_BEAM_GAIN_letter(precoder_total, num_antenna, user_direction, target_direction, sensing_threshold, scaling)

    degree = -90 : 1 : 90;
    num_degree = size(degree, 2);
    beam_gain = zeros(num_degree, 1);
    sensing_threshold = sensing_threshold / scaling^2;
    sensing = zeros(num_degree, 1) + sensing_threshold;

    for i = 1:num_degree

        steering = get_steering_letter(num_antenna, degree(i), 1);
        steering_her = transpose(conj(steering));

        beam_gain(i) = real(steering_her * precoder_total * steering);

    end

    plot(degree, beam_gain);
    xlabel('Degree');
    ylabel('Beam Gain');
    title('Beam Gain vs. Degree');
    grid on;

    hold on;
    plot(user_direction, zeros(size(user_direction,2),1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(target_direction, zeros(size(target_direction,2),1), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    plot(degree, sensing, ':', 'LineWidth', 2);
   
    
end