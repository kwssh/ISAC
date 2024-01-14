function get_received_BEAM_GAIN(W, R, num_user, num_antenna, sensing_th, scaling, distance_user, distance_target, uav_z)

    precoder_total = zeros(num_antenna, num_antenna);

    distance = 100 : 1 : 1000;
    num_distance = size(distance, 2);
    beam_gain = zeros(num_distance, 1);
    sensing_th = sensing_th / scaling^2;
    sensing = zeros(num_distance, 1) + sensing_th;

    for k = 1:num_user
        precoder_total = precoder_total + W(:,:,k);
    end

    precoder_total = precoder_total + R;

    for idx = 1:num_distance

        steering_distance = get_steering(distance(idx), 1, uav_z, num_antenna);
        steering_distance_her = transpose(conj(steering_distance));

        beam_gain(idx) = (steering_distance_her * precoder_total * steering_distance) / distance(idx)^2;
    end

    plot(distance, beam_gain);
    xlabel('Distance');
    ylabel('Beam Gain');
    title('Beam Gain vs. Distance');
    grid on;

    hold on;
    plot(distance, sensing, ':', 'LineWidth', 2);
    plot(distance_user, zeros(size(distance_user,1),1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(distance_target, zeros(size(distance_target,1),1), 'rx', 'MarkerSize', 10, 'LineWidth', 2);

end