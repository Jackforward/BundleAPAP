function [ err ] = bundle_cost( match_table, pt_src, center, NUM_PIC, sigma, gamma )
err = 0; 
pt_ref = match_table(1:3, 3*NUM_PIC+1:end);
Gki_star = exp(-pdist2(center',pt_ref(1:2,:)')./sigma^2);  
W_star = max(gamma,Gki_star); 

for i = 1 : NUM_PIC
    pt_src_ = pt_src{i}';
    mapped_ref = regularize(match_table(1:3, 3*i-2:3*i) * pt_ref); 
    mapped_ref(1,:) = mapped_ref(1,:) .* pt_src_(3,:);
    mapped_ref(2,:) = mapped_ref(2,:) .* pt_src_(3,:);
    mapped_ref(3,:) = mapped_ref(3,:) .* pt_src_(3,:);
    err = err + wdist(pt_src_, mapped_ref, W_star);
end
end

