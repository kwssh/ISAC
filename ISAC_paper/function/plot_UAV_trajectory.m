function fig = plot_UAV_trajectory(uav_t, PARAM)

    hold on
    grid on
    box on
    x = uav_t(:, 1);
    y = uav_t(:, 2);

    plot(x, y, '-^','LineWidth', 3, 'Color', [0.8 0.8 0.8], 'MarkerSize', 7);
    grid on;
    xlabel('X (m)');
    ylabel('Y (m)');
    title('UAV trajectory');

    for i = 1:length(x)
        text(x(i), y(i), [' ' num2str(i)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
    
    plot(PARAM.USER(:,1), PARAM.USER(:,2), 'x', 'MarkerSize', 10, 'LineWidth', 3, 'Color', 'red');

    if PARAM.NUM_TARGET
        plot(PARAM.TARGET(:,1), PARAM.TARGET(:,2), 'o', 'MarkerSize', 10, 'LineWidth', 3, 'Color', 'magenta');
        legend('uav trajectory', 'user position', 'target position');
    else
        legend('uav trajectory', 'user position');
    end
    
    
    
    fig = gcf;
end