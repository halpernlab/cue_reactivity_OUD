function p_val = perm_mdiff(vec1, vec2, nperm)
    % ensure row vectors
    if size(vec1,1) ~= 1, vec1 = vec1'; end
    if size(vec2,1) ~= 1, vec2 = vec2'; end

    both = [vec1, vec2];
    true_diff = mean(vec2) - mean(vec1);

    n1 = length(vec1);
    n2 = length(vec2);
    n_total = n1 + n2;

    perm_diffs = zeros(1, nperm);
    for i = 1:nperm
        idx = randperm(n_total);
        perm_group1 = both(idx(1:n1));
        perm_group2 = both(idx(n1+1:end));
        perm_diffs(i) = mean(perm_group2) - mean(perm_group1);
    end

    zvalue = (true_diff - mean(perm_diffs)) / std(perm_diffs);
    p_val = (1 - normcdf(abs(zvalue)))*2; 
end
