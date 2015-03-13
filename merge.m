function [ match ] = merge( match, sep_match, col1, col2 )
% merge two match tables in to one
temp = [];
for i = 1:size(sep_match,1)
    row = find(match(:, col1) == sep_match(i, 1));
    if isempty(row)
        temp(end+1, col1) = sep_match(i, 1);
        temp(end, col2) = sep_match(i, 2);
    else
        match(row, col2) = sep_match(i, 2);
    end
end
if ~isempty(temp)
    match = [match; temp];
end
end

