function [uav, user_rate] = get_uav_trajectory_BCD_SCA(distance_user, distance_target, num_user, num_target, N, gamma_0, p_max, num_antenna, sensing_th, PARAM, uav_t, V_max, delta_t, A_opt, E_opt, rate_th, isac_duration)

    num_episode_BCD = 10^6;
    user_rate_episode_BCD = zeros(num_user, N, num_episode_BCD);

    distance_user_l = distance_user;
    distance_target_l = distance_target;

    % for episode_BCD = 1 : num_episode_BCD
   
        [uav, z_user, user_rate]= get_trajectory_user_SCA(distance_user_l, distance_target_l, num_user, num_target, N, gamma_0, p_max, num_antenna, sensing_th, PARAM, uav_t, V_max, delta_t, A_opt, E_opt, rate_th, isac_duration);
        % fig3 = plot_UAV_trajectory(uav, PARAM);

        % [z_traget]= get_trajectory_target_SCA(uav, z_user, gamma_0, num_antenna, p_max, A_opt, E_opt, N, num_user, num_target, sensing_th, PARAM, isac_duration, rate_th);

        % distance_target_l = get_distance(PARAM.TARGET, uav, PARAM.UAV_Z);
        % distance_user_l = get_distance(PARAM.USER, uav, PARAM.UAV_Z);
        % 
        % user_rate_episode_BCD(:,:,episode_BCD) = get_user_rate_real(distance_user_l, distance_target_l, num_user, num_target, N, gamma_0, p_max, num_antenna, sensing_th, A_opt, E_opt, uav, PARAM, PARAM.CHANNEL_GAIN);
        % 
        % % break
        % 
        % if episode_BCD > 1
        %     if abs(sum(sum(user_rate_episode_BCD(:,:,episode_BCD))) - sum(sum(user_rate_episode_BCD(:,:,episode_BCD-1)))) < 0.01
        %         break
        %     end
        % end
    % end
end