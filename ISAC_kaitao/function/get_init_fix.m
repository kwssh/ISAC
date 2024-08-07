function [A, E, A_bar, E_bar, uav_init] = get_init_fix(start_x, end_x, uav_y, N, num_user, num_target, isac_duration, PARAM)

    rng(124);

    uav_init_tmp = linspace(start_x, end_x, N);
    uav_init = [uav_init_tmp' ones(N, 1) * uav_y];

    [A, E] = get_association_rule(PARAM, uav_init);

    [A_bar, E_bar] = get_slack_variable(A, E);
end