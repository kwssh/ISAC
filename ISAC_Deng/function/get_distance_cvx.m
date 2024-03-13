function distance = get_distance_cvx(x, uav, distance)

    num_x = size(x,1);
    N = size(uav,1);

    for n = 1 : N
        for k = 1 : num_x
            distance(k,n) = norm([uav(n,3), x(k,1)-uav(n,1), x(k,2)-uav(n,2)]);
        end
    end

end