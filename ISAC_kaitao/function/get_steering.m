function steering = get_steering(distance, num_antenna, uav, user)

    num_antenna_x = sqrt(num_antenna);
    num_antenna_y = sqrt(num_antenna);

    num_user = size(user, 1);
    num_time_slot = size(uav, 1);

    steering = zeros(num_user, num_antenna, num_time_slot);
    steering_x = zeros(num_antenna_x, 1);
    steering_y = zeros(num_antenna_y, 1);

    for n = 1 : num_time_slot
        for k = 1 : num_user
            
            psi = (uav(n,1) - user(k,1)) / distance(k,n);
            omega = (uav(n,2) - user(k,2)) / distance(k,n);

            for x = 1 : num_antenna_x
                steering_x(x) = exp(-1j * pi * (x-1) * psi);
            end

            for y = 1 : num_antenna_y
                steering_y(y) = exp(-1j * pi * (y-1) * omega);
            end
                
            steering(k,:,n) = kron(steering_x, steering_y);

        end
    end
end