function fig = get_received_BEAM_GAIN_eleavtion(W, R, user, uav, target, num_antenna, num_user, uav_z, num_target)
    
    min_x = min([min(uav(:,1)) min(user(:,1)) min(target(:,1))]);
    max_x = max([max(uav(:,1)) max(user(:,1)) max(target(:,1))]);
    min_y = min([min(uav(:,2)) min(user(:,2)) min(target(:,2))]);
    max_y = max([max(uav(:,2)) max(user(:,2)) max(target(:,2))]);

    precoder_total = zeros(num_antenna, num_antenna);

    for k = 1:num_user
        precoder_total = precoder_total + W(:,:,k);
    end

    precoder_total = precoder_total + R;

    x_axis_min = min_x - 50;
    x_axis_max = max_x + 50;

    y_axis_min = min_y - 50;
    y_axis_max = max_y + 50;

    x = linspace(x_axis_min, x_axis_max, x_axis_max - x_axis_min + 1);
    y = linspace(y_axis_min, y_axis_max, y_axis_max - y_axis_min + 1);

    beam_gain = zeros(size(x,2), size(y,2));

    for i = 1 : size(x,2)
        for j = 1 : size(y,2)
            
            distance_tmp = get_distance(uav, [x(i) y(j)], uav_z);

            steering = get_steering(distance_tmp, 1, uav_z, num_antenna);
            steering_her = transpose(conj(steering));

            beam_gain(i,j) = real(steering_her * precoder_total * steering) / distance_tmp^2;
            % beam_gain(i,j) = real(steering_her * precoder_total * steering);
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
    
    if num_target
        plot(target(:, 1), target(:, 2), 'o', 'MarkerSize', 10, 'LineWidth', 3, 'Color', 'magenta');
    end
    
    colorbar;
    % caxis([0, 6]);
    
    xlabel('X (m)');
    ylabel('Y (m)');

    legend('UAV position(initial point)', 'UAV position(final point)', 'User position', 'Target position');
    % legend('UAV position(final point)', 'User position', 'Target position');

    fig = gcf;

end