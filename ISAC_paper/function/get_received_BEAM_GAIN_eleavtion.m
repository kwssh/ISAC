function get_received_BEAM_GAIN_eleavtion(W, R, user, uav, target, num_antenna, num_user, uav_z)
    
    precoder_total = zeros(num_antenna, num_antenna);

    for k = 1:num_user
        precoder_total = precoder_total + W(:,:,k);
    end

    precoder_total = precoder_total + R;

    axis_min = 250;
    axis_max = 750;

    x = linspace(axis_min, axis_max, axis_max - axis_min + 1);
    y = linspace(axis_min, axis_max, axis_max - axis_min + 1);

    beam_gain = zeros(size(x,2), size(y,2));

    for i = 1 : size(x,2)
        for j = 1 : size(y,2)
            
            steering = get_steering(get_distance(uav, [x(i) y(j)], uav_z), 1, uav_z, num_antenna);
            steering_her = transpose(conj(steering));

            beam_gain(i,j) = real(steering_her * precoder_total * steering);
        end
    end

    figure;
    imagesc(x, y, beam_gain');
    set(gca, 'YDir','normal');
    hold on;
    plot(uav(:, 1), uav(:, 2), "+", 'MarkerSize', 10, 'LineWidth', 3, 'Color', 'green');
    plot(user(:, 1), user(:, 2), 'x', 'MarkerSize', 10, 'LineWidth', 3, 'Color', 'red');
    scatter(target(:, 1), target(:, 2), 'ro', 'filled');
    
    colorbar;
    xlabel('X (m)');
    ylabel('Y (m)');

    legend('UAV position', 'User position');
end