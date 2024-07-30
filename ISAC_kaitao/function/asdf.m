function [uav, z_user, user_rate] = asdf(distance_user, distance_target, num_user, num_target, N, gamma_0, p_max, num_antenna, sensing_th, PARAM, uav_t, V_max, delta_t, A_opt, E_opt, rate_th, isac_duration, episilon_sca)

    tic
    num_episode_SCA = 10^6;
    user_rate_episode_SCA = zeros(num_user, N, num_episode_SCA);
    uav_episode = zeros(N, 2, num_episode_SCA);

    scaling = 10000;

    z_user_l = distance_user.^2;
    z_target_l = distance_target.^2;

    for episode_SCA = 1 : num_episode_SCA

        cvx_begin
    
            cvx_solver Mosek
    
            variable z_user(num_user, N)
            variable uav(N, 2)
    
            expressions user_rate(num_user, N)
            expressions user_rate_comm(num_user, N)
            expressions user_rate_ISAC(num_user, num_target, N)
            expressions user_rate_ISAC_target(num_user, num_target, N)
            expressions sensing_constraint(num_target, N)
            expressions distance_target_new(num_target, N)
            expressions distance_user_new(num_user, N)
            expressions user_rate_th(num_user, N)
            expressions user_rate_ISAC_sum(num_user, N)
            expressions user_rate_comm_tmp(num_user, N)
            expressions user_rate_ISAC_tmp(num_user, num_target, N)
            expressions user_rate_comm_tmp_tmp(num_user, N)

            % user_rate = zeros(num_user, N);
            % user_rate_comm = zeros(num_user, N);
            % user_rate_ISAC = zeros(num_user, N);
            % user_rate_ISAC_target = zeros(num_user, num_target, N);
            % sensing_constraint = zeros(num_target, N);
            % distance_target_new = zeros(num_target, N);
            % distance_user_new = zeros(num_user, N);
            % user_rate_th = zeros(num_user, N);
            % user_rate_ISAC_sum = zeros(num_user, N);
            % user_rate_comm_tmp = zeros(num_user, N);
            % user_rate_ISAC_tmp = zeros(num_user, num_target, N);
    
            distance_target_new = get_distance_cvx(PARAM.TARGET, uav, PARAM.UAV_Z, distance_target_new);

            user_rate_comm_tmp = (-gamma_0 * num_antenna * p_max * (z_user - z_user_l) ./ (z_user_l.^2 + gamma_0 * num_antenna * p_max * z_user_l) / log(2)) * scaling;
            user_rate_comm = user_rate_comm_tmp + (log2(1 + gamma_0 * num_antenna * p_max ./ z_user_l)) * scaling;

            user_rate_comm_tmp_tmp = ((rel_entr(z_user / (gamma_0* num_antenna* p_max), z_user / (gamma_0* num_antenna* p_max) + 1) + rel_entr(z_user / (gamma_0* num_antenna* p_max) + 1, z_user / (gamma_0* num_antenna* p_max))) / log(2)) * scaling;

            for k = 1 : num_user

                z_user_repmat = repmat(z_user(k,:), num_target, 1);
                z_user_l_repmat = repmat(z_user_l(k,:), num_target, 1);
                user_rate_comm_tmp_tmp_repmat = repmat(user_rate_comm_tmp_tmp(k,:), num_target, 1);
                user_rate_comm_tmp_tmp_repmat_reshpae = reshape(user_rate_comm_tmp_tmp_repmat, [1 size(user_rate_comm_tmp_tmp_repmat)]);

                user_rate_ISAC_tmp(k,:,:) = (-gamma_0 * (num_antenna * p_max - z_target_l * sensing_th) .* (z_user_repmat - z_user_l_repmat) ./ (z_user_l_repmat.^2 + gamma_0 * (num_antenna * p_max - z_target_l * sensing_th) .* z_user_l_repmat) / log(2)) *scaling;
                
                user_rate_ISAC_const = (log2(1 + gamma_0 * (num_antenna * p_max - z_target_l * sensing_th) ./ z_user_l_repmat)) * scaling;
                user_rate_ISAC_const_reshpae = reshape(user_rate_ISAC_const, [1 size(user_rate_ISAC_const)]);
                
                user_rate_ISAC_target(k,:,:) = user_rate_ISAC_tmp(k,:,:) + user_rate_ISAC_const_reshpae;

                user_rate_ISAC(k,:,:) = E_opt(k,:,:) .* (real(user_rate_ISAC_target(k,:,:)) - user_rate_comm_tmp_tmp_repmat_reshpae);
            end

            user_rate_ISAC_target_sum = squeeze(sum(user_rate_ISAC, 2));

            user_rate = A_opt .* user_rate_comm + user_rate_ISAC_target_sum;

            %--------------------------------------------------------------------------------------------------------------------------------------------------------------------%

            user_expand = reshape(PARAM.USER, num_user, 1, 2);
            user_expand = repmat(user_expand, 1, N, 1);

            uav_expand = reshape(uav, 1, N, 2);
            uav_expand = repmat(uav_expand, num_user, 1, 1);
            
            distance_user_uav_2D = user_expand - uav_expand;
            distance_user_uav_2D_square = sum(distance_user_uav_2D.^2, 3);
            distance_user_uav_square = distance_user_uav_2D_square + PARAM.Z^2;

            z_user * scaling >= distance_user_uav_square * scaling;
            
            sensing_constraint_repmat = (num_antenna * p_max - pow_pos(distance_target_new, 2) * sensing_th) * scaling;
            sensing_constraint_repmat = repmat(reshape(sensing_constraint_repmat, 1 ,size(sensing_constraint_repmat, 1), size(sensing_constraint_repmat, 2)), num_user, 1, 1);
            
            sensing_constraint_tmp = E_opt .* sensing_constraint_repmat;
            sensing_constraint = sum(sensing_constraint_tmp, 1);

            sensing_constraint >= 0;
                  
            uav(1,1) == uav_t(1,1);
            uav(1,2) == uav_t(1,2);

            uav(N,1) == uav_t(N,1);
            uav(N,2) == uav_t(N,2);

            uav_diff = uav(2:N, :) - uav(1:N-1, :);

            for n = 1 : N - 1
                norm(uav_diff(n,:)) <= V_max * delta_t;
            end

            maximize(sum(sum(user_rate)))

            % for i = 1:isac_duration:N
            %     sum(user_rate(:,i:i+isac_duration-1),2)/ isac_duration >= rate_th * scaling;
            % end
    
        cvx_end

        toc

        if strcmp('Infeasible', cvx_status)
            disp("qweqwew")
        end

        if strcmp('Failed', cvx_status)
            disp("qweqwew")
        end

        user_rate_episode_SCA(:,:,episode_SCA) = user_rate;
        uav_episode(:,:,episode_SCA) = uav;

        % break

        if episode_SCA > 1
            if sum(sum(user_rate_episode_SCA(:,:,episode_SCA))) - sum(sum(user_rate_episode_SCA(:,:,episode_SCA-1))) <= episilon_sca * scaling
                break
            end
        end
        
        z_user_l = z_user;
        % z_user_l = get_distance(PARAM.USER, uav, PARAM.UAV_Z).^2;
    end
end