function [ distance ] = dist( match_up_src, mapped_ref )
distance = sum(sqrt(sum((match_up_src - mapped_ref).^2)));

end

