function werr = mdlt_cost(center, h, match_up_src, match_up_ref, sigma, gamma)
    mapped_ref = regularize(h * match_up_ref);
    Gki_star = exp(-pdist2(center',match_up_ref(1:2,:)')./sigma^2);  
    W_star = max(gamma,Gki_star); 
    werr = wdist(match_up_src,mapped_ref,W_star);