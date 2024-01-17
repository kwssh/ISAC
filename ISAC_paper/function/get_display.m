function get_display(data, data_name)

    [row, col] = size(data);

    fprintf(data_name);
    for i = 1:row
        for j = 1:col
            fprintf('%.4f ', data(i, j));
        end
        fprintf('\n');

        if i < row
            fprintf(repmat(' ', 1, length(data_name)));
        end
    end
end