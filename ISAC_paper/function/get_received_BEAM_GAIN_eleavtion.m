function get_received_BEAM_GAIN_eleavtion(W, R, user, uav, target, num_antenna, num_user, uav_z)
    
    precoder_total = zeros(num_antenna, num_antenna);

    for k = 1:num_user
        precoder_total = precoder_total + W(:,:,k);
    end

    precoder_total = precoder_total + R;

    axis_min = -350;
    axis_max = 700;

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
    q = [450 525];
    plot(q(:, 1), q(:, 2), "+", 'MarkerSize', 10, 'LineWidth', 3, 'Color', 'black');
    plot(uav(:, 1), uav(:, 2), "+", 'MarkerSize', 10, 'LineWidth', 3, 'Color', 'green');
    plot(user(:, 1), user(:, 2), 'x', 'MarkerSize', 10, 'LineWidth', 3, 'Color', 'red');
    plot(target(:, 1), target(:, 2), 'o', 'MarkerSize', 10, 'LineWidth', 3, 'Color', 'magenta');
    
    colorbar;
    % caxis([0, 6]);
    
    xlabel('X (m)');
    ylabel('Y (m)');

    legend('UAV position(initial point)', 'UAV position(final point)', 'User position', 'Target position');
    % legend('UAV position(final point)', 'User position', 'Target position');

end