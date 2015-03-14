function [ distance ] = wdist( match_up_src, mapped_ref, W_star )
if(isempty(W_star))
    distance = 0;
else
    distance = (sqrt(sum((match_up_src - mapped_ref).^2))) * W_star' ;
end
end

