function [matchab_inliar, H, A, T1, T2, D1, D2] = ransac(valid_ptsa, valid_ptsb, matchab)
X1 = valid_ptsa(matchab(:,1)).Location;X1(:,3) = 1;X1 = double(X1)';
X2 = valid_ptsb(matchab(:,2)).Location;X2(:,3) = 1;X2 = double(X2)';
numMatches = size(matchab,1);
for t = 1:200
    %estimate homography(need more detail)
    subset = vl_colsubset(1:numMatches,4);
    A = [];
    for i = subset
        A = cat(1,A,kron(X1(:,i)',vl_hat(X2(:,i)))); %???
    end
    [U,S,V] = svd(A);
    H{t} = reshape(V(:,9),3,3);
    %score homography
    X2_ = H{t} * X1;
    du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:);
    dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:);
    ok{t} = (du.*du + dv.*dv) < 6*6;
    score(t) = sum(ok{t});
end
[score,best] = max(score);
%H = double(H{best});
ok = ok{best};
confidence = score / (8 + 0.3 * numMatches);
if confidence < 2
    matchab_inliar = [];
else
    matchab_inliar = matchab(ok,:);
end
% Refine homography using DLT on inliers.
X1 = X1(:,ok);X2 = X2(:,ok);
[ dat_norm_img1,T1 ] = normalise2dpts(X1);
[ dat_norm_img2,T2 ] = normalise2dpts(X2);
data_norm = [ dat_norm_img1 ; dat_norm_img2 ];
[ h,A,D1,D2 ] = feval('homography_fit',data_norm(:,:));
H = T2\(reshape(h,3,3)*T1);

end