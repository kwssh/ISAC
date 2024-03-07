function [channel_matrix_user_DL, channel_matrix_user_UL, channel_matrix_target, channel_matrix_target_diff] = get_channel(C, D, K, channel_gain, num_antenna, distance_user, uav, distance_target, RCS, X_DL_old)
    rng(123)
    
    num_user = size(distance_user, 1);
    num_target = size(distance_target, 1);
    N = size(distance_user, 2);

    theta_user = zeros(num_user, N);
    theta_target = zeros(num_target, N);

    prob_LoS_user = zeros(num_user, N);
    prob_n_LoS_user = zeros(num_user, N);

    prob_LoS_target = zeros(num_target, N);
    prob_LoS_target_dev = zeros(num_target, N);

    prob_n_LoS_target = zeros(num_target, N);

    channel_matrix_user_DL = zeros(num_antenna, num_antenna, num_user, N);
    channel_matrix_user_UL = zeros(num_antenna, num_antenna, num_user, N);
    channel_user_DL = zeros(num_antenna, num_user, N);
    channel_her_user_DL = zeros(num_user, num_antenna, N);
    channel_user_UL = zeros(num_antenna, num_user, N);
    channel_her_user_UL = zeros(num_user, num_antenna, N);
    channel_n_LoS_user_DL = (randn(num_antenna, num_user, N) + 1i * randn(num_antenna, num_user, N)) / sqrt(2);
    channel_n_LoS_user_UL = (randn(num_antenna, num_user, N) + 1i * randn(num_antenna, num_user, N)) / sqrt(2);

    channel_matrix_target = zeros(num_antenna, num_antenna, num_target, N);
    channel_matrix_target_diff = zeros(num_antenna, num_antenna, num_target, N);

    channel_n_LoS_target = (randn(num_antenna, num_antenna, num_target, N) + 1i * randn(num_antenna, num_antenna, num_target, N)) / sqrt(2);

    for n = 1 : N
        theta_user(:, n) = asin(uav(n, 3) ./ distance_user(:, n));
        theta_target(:, n) = asin(uav(n, 3) ./ distance_target(:, n));

        prob_LoS_user(:, n) = 1 ./ (1 + C * exp(-D * ((180 / pi)* theta_user(:, n) - C)));
        prob_n_LoS_user(:, n) = 1 - prob_LoS_user(:, n);

        prob_LoS_target(:, n) = 1 ./ (1 + C * exp(-D * ((180 / pi)* theta_target(:, n) - C)));
        prob_LoS_target_dev(:, n) = 0.5 * (1 + C * exp(-D * ((180 / pi)* theta_target(:, n) - C))) .^ (-1.5) .* (C * D * exp(-D * ((180 / pi)* theta_target(:, n) - C)));

        prob_n_LoS_target(:, n) = 1 - prob_LoS_target(:, n);

        for k = 1 : num_user
            steering_tmp_user = [0, -1i * pi * sin(theta_user(k, n)) * (1 : num_antenna - 1)];

            channel_user_DL_tmp = sqrt(prob_LoS_user(k, n) * channel_gain) * (transpose(exp(steering_tmp_user))) + sqrt(prob_n_LoS_user(k, n) * channel_gain * K) * channel_n_LoS_user_DL(:, k, n);
            channel_user_DL(:, k, n) = channel_user_DL_tmp ./ distance_user(k, n);
            channel_her_user_DL(k, :, n) = channel_user_DL(:, k, n)';
            channel_matrix_user_DL(:, :, k, n) = channel_user_DL(:, k, n) * channel_her_user_DL(k, :, n);

            channel_user_UL_tmp = sqrt(prob_LoS_user(k, n) * channel_gain) * (transpose(exp(steering_tmp_user))) + sqrt(prob_n_LoS_user(k, n) * channel_gain * K) * channel_n_LoS_user_UL(:, k, n);
            channel_user_UL(:, k, n) = channel_user_UL_tmp ./ distance_user(k, n);
            channel_her_user_UL(k, :, n) = channel_user_UL(:, k, n)';
            channel_matrix_user_UL(:, :, k, n) = channel_user_UL(:, k, n) * channel_her_user_UL(k, :, n);
        end

        for j = 1 : num_target
            steering_tmp_target = [0, -1i * pi * sin(theta_target(j, n)) * (1 : num_antenna - 1)];
            steering_tmp_target_diff = [0, 1i * pi * cos(theta_target(j, n)) * (1 : num_antenna - 1)];
            steering_target_diff = toeplitz(steering_tmp_target_diff);

            channel_LoS_target_tmp = transpose(exp(steering_tmp_target));
            channel_LoS_target = channel_LoS_target_tmp * channel_LoS_target_tmp';

            channel_matrix_target(:, :, j, n) = sqrt(prob_LoS_target(j, n)) * channel_LoS_target + sqrt(prob_n_LoS_target(j, n) * K) * channel_n_LoS_target(:,:,j,n);
            channel_matrix_target(:, :, j, n) = RCS * channel_matrix_target(:, :, j, n) / (2 * distance_target(j, n));

            channel_matrix_target_diff(:, :, j, n) = prob_LoS_target_dev(j, n) * channel_LoS_target + sqrt(prob_LoS_target(j, n)) * channel_LoS_target .* steering_target_diff;
        end
    end
end