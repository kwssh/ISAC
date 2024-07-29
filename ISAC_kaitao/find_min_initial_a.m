function min_initial_a = find_min_initial_a()
    % 초기값 설정
    a0_lower = 0;
    a0_upper = 1;
    tolerance = 1e-6;
    max_iterations = 1000;

    % 이분법을 사용하여 초기값 탐색
    while (a0_upper - a0_lower) > tolerance
        a0 = (a0_lower + a0_upper) / 2;
        if iterate(a0, tolerance, max_iterations)
            a0_upper = a0;
        else
            a0_lower = a0;
        end
    end

    min_initial_a = a0;
end