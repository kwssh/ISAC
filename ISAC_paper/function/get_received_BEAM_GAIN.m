function get_received_BEAM_GAIN(W, R, num_user, num_antenna, sensing_th, scaling, distance_user, distance_target, uav_z)

    precoder_total = zeros(num_antenna, num_antenna);
    [num_target, ~] = size(distance_target);
    degree_user = zeros(num_user, 1);
    degree_target = zeros(num_target, 1);

    
    for i = 1 : num_user
        degree_user(i) = rad2deg(acos(uav_z / distance_user(i)));
    end

    for j = 1 : num_target
        degree_target(j) = rad2deg(acos(uav_z / distance_target(j)));
    end

    
    degree = 0 : 1 : 360;
    num_distance = size(degree, 2);
    beam_gain = zeros(num_distance, 1);
    sensing_th = sensing_th / scaling^2;
    sensing = zeros(num_distance, 1) + sensing_th;

    for k = 1:num_user
        precoder_total = precoder_total + W(:,:,k);
    end

    precoder_total = precoder_total + R;

    for idx = 1:num_distance

        steering_distance = get_steering_degree(degree(idx), 1, num_antenna);
        steering_distance_her = transpose(conj(steering_distance));

        beam_gain(idx) = (steering_distance_her * precoder_total * steering_distance);
    end
    
    plot(degree, beam_gain);
    xlabel('Degree');
    ylabel('Beam Gain');
    title('Beam Gain vs. Degree');
    grid on;

    hold on;
    plot(degree, sensing, ':', 'LineWidth', 2);
    plot(degree_user, zeros(size(distance_user,1),1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(degree_target, zeros(size(distance_target,1),1), 'rx', 'MarkerSize', 10, 'LineWidth', 2);

    % polarplot(deg2rad(degree), beam_gain);
    % title('Beam Gain vs. Angle');
end