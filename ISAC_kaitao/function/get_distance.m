function distance = get_distance(x, uav, uav_z)

    num_x = size(x,1);
    num_uav = size(uav,1);

    distance = zeros(num_x, num_uav);

    for k = 1 : num_x
        for j = 1 : num_uav
            distance(k,j) = norm([uav_z, x(k,1)-uav(j,1), x(k,2)-uav(j,2)]);
        end
    end

end