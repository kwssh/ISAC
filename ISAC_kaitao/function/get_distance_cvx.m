function distance = get_distance_cvx(x, uav, uav_z, distance)

    num_x = size(x,1);
    num_time_slot = size(uav,1);

    for k = 1 : num_x
        for j = 1 : num_time_slot
            distance(k,j) = norm([uav_z, x(k,1)-uav(j,1), x(k,2)-uav(j,2)]);
        end
    end

end