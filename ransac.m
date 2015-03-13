function [matchab_inliar, H] = ransac(valid_ptsa, valid_ptsb, matchab)
X1 = valid_ptsa(matchab(:,1)).Location;X1(:,3) = 1;X1 = X1';
X2 = valid_ptsb(matchab(:,2)).Location;X2(:,3) = 1;X2 = X2';
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
H = double(H{best});
ok = ok{best};
confidence = score / (8 + 0.3 * numMatches);
if confidence < 2
    matchab_inliar = [];
else
    matchab_inliar = matchab(ok,:);
end
end