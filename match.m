function [ H,inliar,confidence,f1,f2 ] = match( I1,I2 )
%%
%params
showmatch = false;

%read image
% I1 = imread('E:\stitching\source\low_res\G9PQ0282.jpg');
% I2 = imread('E:\stitching\source\low_res\G9PQ0283.jpg');

%make single and grayscale
I1g = im2single(rgb2gray(I1));
I2g = im2single(rgb2gray(I2));

%extract features
[f1,d1] = vl_sift(I1g);
[f2,d2] = vl_sift(I2g);

%number of features
for i = 1: size(f1,2)
    f1(5,i) = i;
end
for i = 1: size(f2,2)
    f2(5,i) = i;
end

%match features
[matches,scores] = vl_ubcmatch(d1,d2);
numMatches = size(matches,2);

X1 = f1(1:2,matches(1,:)); X1(3,:)=1;  %keypoint position
X2 = f2(1:2,matches(2,:)); X2(3,:)=1;

%RANSAC
clear H score ok;
for t = 1:100
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
H = H{best};
ok = ok{best};

X1(4,:)=f1(5,matches(1,:));
X2(4,:)=f2(5,matches(2,:));

inliar = [X1(1:2,ok);X1(4,ok);X2(1:2,ok);X2(4,ok)];
confidence = score / (8 + 0.3 * numMatches);


%(opt)show matches
if showmatch
    dh1 = max(size(I2,1)-size(I1,1),0) ;
    dh2 = max(size(I1,1)-size(I2,1),0) ;
    
    figure(1) ; clf ;
    subplot(2,1,1) ;
    imagesc([padarray(I1,dh1,'post') padarray(I2,dh2,'post')]) ;
    o = size(I1,2) ;
    line([f1(1,matches(1,:));f2(1,matches(2,:))+o], ...
        [f1(2,matches(1,:));f2(2,matches(2,:))]) ;
    title(sprintf('%d tentative matches', numMatches)) ;
    axis image off ;
    
    subplot(2,1,2) ;
    imagesc([padarray(I1,dh1,'post') padarray(I2,dh2,'post')]) ;
    o = size(I1,2) ;
    line([f1(1,matches(1,ok));f2(1,matches(2,ok))+o], ...
        [f1(2,matches(1,ok));f2(2,matches(2,ok))]) ;
    title(sprintf('%d (%.2f%%) inliner matches out of %d', ...
        sum(ok), ...
        100*sum(ok)/numMatches, ...
        numMatches)) ;
    axis image off ;
    
    drawnow ;
    
end



end

