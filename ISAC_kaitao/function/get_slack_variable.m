function [A_bar_opt, E_bar_opt] = get_slack_variable(A_opt, E_opt)

    A_bar_opt = (A_opt + A_opt.^2) ./ (1 + A_opt.^2);
    E_bar_opt = (E_opt + E_opt.^2) ./ (1 + E_opt.^2);

end