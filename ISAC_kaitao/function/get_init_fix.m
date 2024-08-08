function [A, E, A_bar, E_bar, uav_init] = get_init_fix(start_x, end_x, uav_y, N, num_user, num_target, isac_duration, PARAM)

    rng(124);

    mid_x = (start_x + end_x) / 2;

    distance_to_mid = abs(mid_x - start_x);
    distance_to_end = abs(end_x - mid_x);

    time_slot_to_mid = ceil(distance_to_mid / (PARAM.V_MAX * PARAM.TOTAL_DURATION)) + 1;
    time_slot_to_end = ceil(distance_to_end / (PARAM.V_MAX * PARAM.TOTAL_DURATION)) + 1;
    time_slot_at_mid = PARAM.TOTAL_TIME_SLOT - time_slot_to_mid -time_slot_to_end;

    uav_init_tmp = [linspace(start_x, mid_x, time_slot_to_mid), ...
                repmat(mid_x, 1, time_slot_at_mid), ...
                linspace(mid_x, end_x, time_slot_to_end)];



    % uav_init_tmp = linspace(start_x, end_x, N);
    uav_init = [uav_init_tmp' ones(N, 1) * uav_y];

    [A, E] = get_association_rule(PARAM, uav_init);

    [A_bar, E_bar] = get_slack_variable(A, E);
end