function uav = get_uav_trajectory_BCD_SCA(distance_user, distance_target, num_user, num_target, N, gamma_0, p_max, num_antenna, sensing_th, PARAM, uav_t, V_max, delta_t, A_opt, E_opt, rate_th, isac_duration)

    num_episode_BCD = 10^6;
    user_rate_episode_BCD = zeros(num_user, N, num_episode_BCD);

    distance_user_l = distance_user;
    distance_target_l = distance_target;
   
    for episode_BCD = 1 : num_episode_BCD

        uav= get_trajectory_user_SCA(distance_user_l, distance_target_l, num_user, num_target, N, gamma_0, p_max, num_antenna, sensing_th, PARAM, uav_t, V_max, delta_t, A_opt, E_opt, rate_th, isac_duration);
        % fig3 = plot_UAV_trajectory(uav, PARAM);
        distance_target_l = get_distance(PARAM.TARGET, uav, PARAM.UAV_Z);
        distance_user_l = get_distance(PARAM.USER, uav, PARAM.UAV_Z);

        user_rate = get_user_rate_real(distance_user_l, distance_target_l, num_user, num_target, N, gamma_0, p_max, num_antenna, sensing_th, A_opt, E_opt, uav, PARAM, PARAM.CHANNEL_GAIN);

        uav_t = [
            270, 200;  % 시작점 q_I
            330, 350;  % 첫 번째 꼭지점
            370, 450;  % 두 번째 꼭지점
            434, 350;  % 두 번째 꼭지점
            485.2, 270;  % 세 번째 꼭지점
            530, 200;  % 끝점 q_F
                ];
        % 
        % distance_target_l = get_distance(PARAM.TARGET, uav_t, PARAM.UAV_Z);
        % distance_user_l = get_distance(PARAM.USER, uav_t, PARAM.UAV_Z);
        % 
        % user_rate = get_user_rate_real(distance_user_l, distance_target_l, num_user, num_target, N, gamma_0, p_max, num_antenna, sensing_th, A_opt, E_opt, uav_t, PARAM, PARAM.CHANNEL_GAIN);





        user_rate_episode_BCD(:,:,episode_BCD) = get_objective(distance_user_l, distance_target_l, num_user, num_target, N, gamma_0, p_max, num_antenna, sensing_th, A_opt, E_opt);
     
        if episode_BCD > 1
            if abs(sum(sum(user_rate_episode_BCD(:,:,episode_BCD))) - sum(sum(user_rate_episode_BCD(:,:,episode_BCD-1)))) < 0.01
                break
            end
        end


    end



    

end