clear

num_antenna = 16;
num_user = 2;
num_target = 2;
scaling = 1000;

p_max = 2;
p_circuit = 10^(-0.5);
noise_power = 10^(-11) * scaling^2;

path_loss_gain = sqrt(10^(-9.9));

sinr_threshold = 10^(0.5);
sensing_threshold = 10^(-2) * scaling^2 * 0;
pho = 0.35;
dynamic_power = 10^(-5.6);
epsilon = 10^(-3);

user_direction = [-30 20];
target_direction = [-54, 18];

channel = zeros(num_antenna, num_user);
channel_her = zeros(num_user, num_antenna);

sterring_target = zeros(num_antenna, num_target);
sterring_target_her = zeros(num_target, num_antenna);

num_episode = 100;
sum_rate_episode = zeros(num_user, num_episode);

for k = 1:num_user

    channel(:,k) = path_loss_gain * get_steering_letter(num_antenna, user_direction(k), scaling);
    channel_her(k,:) = transpose(conj(channel(:,k)));
end

for m = 1:num_target

    sterring_target(:,m) = get_steering_letter(num_antenna, target_direction(m), scaling);
    sterring_target_her(m,:) = transpose(conj(sterring_target(:,m)));
end

[V_l, lambda_l] = get_init_letter(num_user, num_antenna, num_target, channel, channel_her, sinr_threshold, noise_power, sterring_target, sterring_target_her, sensing_threshold, p_max, p_circuit, pho, scaling);

sum_rate_current = get_SINR_letter(V_l, channel, channel_her, num_user, noise_power);
sum_rate_episode(:,1) = sum_rate_current;

for i = 2:num_episode

    sum_rate_prev = sum_rate_current;

    cvx_begin

        cvx_solver Mosek
    
        variable V(num_antenna, num_antenna, num_user + 1) hermitian semidefinite;
        variable t
        variable u

        expressions interference_constraint(num_user, 1)
        expressions interference(num_user, 1)
        expressions non_interference(num_user, 1)
        expressions sensing_power(num_target, 1)
        expressions sum_rate_right(num_user, 1)
        expressions sum_rate_left(num_user, 1)

        V_total = get_precoder_total_letter(num_user + 1, V);
    
        maximize(t)
    
        subject to

            interference_const = get_interference_const_letter(interference_constraint, interference, non_interference, num_user, channel, channel_her, V, sinr_threshold, noise_power);
            sensing_power_const = get_sensing_power_const_letter(sensing_power, num_target, sterring_target, sterring_target_her, V_total, sensing_threshold);
            sum_rate = get_sum_rate_const_letter(sum_rate_right, sum_rate_left, num_user, V, V_l, channel, channel_her, noise_power);

            0.5 * (lambda_l * t^2 + u^2 / lambda_l) - sum_rate <= 0;
            u >= trace(V_total) / pho + p_circuit;
            interference_const >= 0;
            real(trace(V_total)) <= p_max;
            sensing_power_const >= 0;
            t >= 0;
            u >= 0;

    cvx_end

    V_opt = get_precoder_opt_letter(channel, channel_her, V, num_user, num_antenna);
    V_total_opt = get_precoder_total_letter(num_user + 1, V_opt);

    sum_rate_current = get_SINR_letter(V_opt, channel, channel_her, num_user, noise_power);
    sum_rate_episode(:,i) = sum_rate_current;

    if abs(sum(sum_rate_current) - sum(sum_rate_prev)) < epsilon
        break;
    end

    V_l = V_opt;
    lambda_l = u / t;

end

get_BEAM_GAIN_letter(V_total_opt, num_antenna, user_direction, target_direction, sensing_threshold, scaling);


    
    
    