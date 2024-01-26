function target = get_target(num_target)

    rng(123);

    if num_target == 1
        target = [randi([450, 550], num_target, 1) randi([590, 610], num_target, 1)];

    else
        target_x_1 = randi([450, 550], num_target / 2, 1);
        target_x_2 = 1000 - target_x_1;
       
        target_y = randi([590, 610], num_target / 2, 1);
        
        target_1 = [target_x_1 target_y];
        target_2 = [target_x_2 target_y];

        target = [target_1; target_2];
    end
end