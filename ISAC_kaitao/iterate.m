function converged = iterate(a0, tolerance, max_iterations)
    a = a0;
    target_a = 1;
    for i = 1:max_iterations
        b = (a + a^2) / (1 + a^2);
        a = b / (b^2 - 2*b + 2);
        if abs(a - target_a) < tolerance
            converged = true;
            return;
        end
    end
    converged = false;
end