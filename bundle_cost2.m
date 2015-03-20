function err = bundle_cost2(match_table, pt_src1, pt_src3, pt_ref_orig, center, sigma, gamma)

h1 = match_table(1:3,1:3);
h3 = match_table(1:3,4:6);
pt_ref = match_table(1:3,7:end);

pt_src1 = pt_src1';
pt_src3 = pt_src3';
pt_ref_orig = pt_ref_orig';

Gki_star = exp(-pdist2(center',pt_ref(1:2,:)')./sigma^2);  
W_star = max(gamma,Gki_star); 

mapped_ref1 = regularize(h1 * pt_ref);
% mapped_ref1 = mapped_ref1(:,size(pt_src1,2));
mapped_ref1(1,:) = mapped_ref1(1,:) .* pt_src1(3,:);
mapped_ref1(2,:) = mapped_ref1(2,:) .* pt_src1(3,:);
mapped_ref1(3,:) = mapped_ref1(3,:) .* pt_src1(3,:);

mapped_ref3 = regularize(h3 * pt_ref);
% mapped_ref3 = mapped_ref3(:,size(pt_src3,2));
mapped_ref3(1,:) = mapped_ref3(1,:) .* pt_src3(3,:);
mapped_ref3(2,:) = mapped_ref3(2,:) .* pt_src3(3,:);
mapped_ref3(3,:) = mapped_ref3(3,:) .* pt_src3(3,:);

err = wdist(pt_src1, mapped_ref1, W_star) + wdist(pt_src3, mapped_ref3, W_star) + wdist(pt_ref_orig, pt_ref, W_star);
end