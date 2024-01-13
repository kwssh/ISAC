function all_path = depth_first_search(node, visited, adjList, path_tmp, all_path, end_node)

    visited(node) = true;
    path_tmp(end + 1) = node;
    is_leaf = true;

    for i = 1:length(adjList{node})

        next_node = adjList{node}(i);

        if ~visited(next_node)
            is_leaf = false;
            all_path = depth_first_search(next_node, visited, adjList, path_tmp, all_path, end_node);

        end
    end

    if ~isempty(all_path)
        return;
    end

    if is_leaf && node == end_node
        all_path{end+1} = path_tmp;
    end

    

end