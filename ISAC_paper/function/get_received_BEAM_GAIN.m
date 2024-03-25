function fig = get_received_BEAM_GAIN(W, R, num_user, num_antenna, sensing_th, scaling, distance_user, distance_target, uav_z)

    precoder_total = zeros(num_antenna, num_antenna);
    [num_target, ~] = size(distance_target);
    degree_user = zeros(num_user, 1);
    degree_target = zeros(num_target, 1);

    degree = 0 : 1 : 360;
    num_distance = size(degree, 2);
    beam_gain = zeros(num_distance, 1);
    beam_gain_user = zeros(num_user, 1);
    beam_gain_target = zeros(num_target, 1);

    sensing_th = sensing_th / scaling^2;
    sensing = zeros(num_distance, 1) + sensing_th;

    for k = 1:num_user
        precoder_total = precoder_total + W(:,:,k);
    end

    precoder_total = precoder_total + R;

    
    for i = 1 : num_user
        degree_user(i) = rad2deg(acos(uav_z / distance_user(i)));

        steering_distance = get_steering_degree(degree_user(i), 1, num_antenna);
        steering_distance_her = transpose(conj(steering_distance));

        distance = uav_z / (cos(deg2rad(degree_user(i))));

        beam_gain_user(i) = (steering_distance_her * precoder_total * steering_distance) / distance^2;
    end

    for j = 1 : num_target
        degree_target(j) = rad2deg(acos(uav_z / distance_target(j)));

        steering_distance = get_steering_degree(degree_target(j), 1, num_antenna);
        steering_distance_her = transpose(conj(steering_distance));

        distance = uav_z / (cos(deg2rad(degree_target(j))));

        beam_gain_target(j) = (steering_distance_her * precoder_total * steering_distance) / distance^2;
    end

    for idx = 1:num_distance

        steering_distance = get_steering_degree(degree(idx), 1, num_antenna);
        steering_distance_her = transpose(conj(steering_distance));

        distance = uav_z / (cos(deg2rad(degree(idx))));

        beam_gain(idx) = (steering_distance_her * precoder_total * steering_distance) / distance^2;
        % beam_gain(idx) = (steering_distance_her * precoder_total * steering_distance);
    end

    plot(degree, beam_gain);

    % beam_gain_user_tmp = linspace(0, beam_gain_user, 100);
    % plot(ones(size(beam_gain_user_tmp)) * degree_user, beam_gain_user_tmp, '--');
    % plot(degree_user, beam_gain_user, 'go', 'MarkerSize', 10, 'LineWidth', 2);
    % 
    % plot(degree_target, beam_gain_target, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    xlim([0 361]);
    xlabel('Degree');
    ylabel('Beam Gain');
    title('Beam Gain vs. Degree');
    grid on;

    hold on;
    plot(degree, sensing, ':', 'LineWidth', 2);
    plot(degree_user, zeros(size(distance_user,1),1), 'x', 'MarkerSize', 10, 'LineWidth', 3, 'Color', 'red');
    plot(degree_target, zeros(size(distance_target,1),1), 'o', 'MarkerSize', 10, 'LineWidth', 3, 'Color', 'magenta');
    
    legend('', '', 'User position', 'Target position');
    fig = gcf;

    % polarplot(deg2rad(degree), beam_gain);
    % title('Beam Gain vs. Angle');
end