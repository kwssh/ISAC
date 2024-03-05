function [channel_matrix_user, channel_matrix_target] = get_channel(C, D, K, channel_gain, num_antenna, distance_user, uav, distance_target, RCS)
    rng(123)
    
    num_user = size(distance_user, 1);
    num_target = size(distance_target, 1);
    N = size(distance_user, 2);

    theta_user = zeros(num_user, N);
    theta_target = zeros(num_target, N);

    prob_LoS_user = zeros(num_user, N);
    prob_n_LoS_user = zeros(num_user, N);

    prob_LoS_target = zeros(num_target, N);
    prob_n_LoS_target = zeros(num_target, N);

    channel_matrix_user = zeros(num_antenna, num_antenna, num_user, N);
    channel_user = zeros(num_antenna, num_user, N);
    channel_her_user = zeros(num_user, num_antenna, N);
    channel_n_LoS_user = (randn(num_antenna, num_user, N) + 1i * randn(num_antenna, num_user, N)) / sqrt(2);

    channel_matrix_target = zeros(num_antenna, num_antenna, num_target, N);
    channel_n_LoS_target = (randn(num_antenna, num_antenna, num_target, N) + 1i * randn(num_antenna, num_antenna, num_target, N)) / sqrt(2);

    for n = 1 : N
        theta_user(:, n) = asin(uav(n, 3) ./ distance_user(:, n));
        theta_target(:, n) = asin(uav(n, 3) ./ distance_target(:, n));

        prob_LoS_user(:, n) = 1 ./ (1 + C * exp(-D * ((180 / pi)* theta_user(:, n) - C)));
        prob_n_LoS_user(:, n) = 1 - prob_LoS_user(:, n);

        prob_LoS_target(:, n) = 1 ./ (1 + C * exp(-D * ((180 / pi)* theta_target(:, n) - C)));
        prob_n_LoS_target(:, n) = 1 - prob_LoS_target(:, n);

        for k = 1 : num_user
            steering_tmp_user = [0, -1i * pi * sin(theta_user(k, n)) * (1 : num_antenna - 1)];

            channel_user(:, k, n) = sqrt(prob_LoS_user(k, n) * channel_gain) * (transpose(exp(steering_tmp_user))) + sqrt(prob_n_LoS_user(k, n) * channel_gain * K) * channel_n_LoS_user(:, k, n);
            channel_user(:, k, n) = channel_user(:, k, n) ./ distance_user(k, n);
            channel_her_user(k, :, n) = channel_user(:, k, n)';
            channel_matrix_user(:, :, k, n) = channel_user(:, k, n) * channel_her_user(k, :, n);
        end

        for j = 1 : num_target
            steering_tmp_target = [0, -1i * pi * sin(theta_target(j, n)) * (1 : num_antenna - 1)];

            channel_LoS_target_tmp = transpose(exp(steering_tmp_target));
            channel_LoS_target = channel_LoS_target_tmp * channel_LoS_target_tmp';

            channel_matrix_target(:, :, j, n) = sqrt(prob_LoS_target(j, n)) * channel_LoS_target + sqrt(prob_n_LoS_target(j, n) * K) * channel_n_LoS_target(:,:,j,n);
            channel_matrix_target(:, :, j, n) = RCS * channel_matrix_target(:, :, j, n) / (2 * distance_target(j, n));
        end
    end
end