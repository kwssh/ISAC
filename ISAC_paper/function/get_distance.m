function distance = get_distance(x, y, uav_z)

    distance = norm([uav_z, x(1,1)-y(1,1), x(1,2)-y(1,2)]);

end