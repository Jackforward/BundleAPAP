function [linear_mdlt,  Hmdlt]= stitching(img1, img2, data_orig, off, cw, ch)

data_orig = double(data_orig);
gamma = 0.01; % Normalizer for Moving DLT. (0.0015-0.1 are usually good numbers).
sigma = 8.5;  % Bandwidth for Moving DLT. (Between 8-12 are good numbers).   

C1 = 100; % Resolution/grid-size for the mapping function in MDLT (C1 x C2).
C2 = 100;

[ dat_norm_img1,T1 ] = normalise2dpts(data_orig(1:3,:));
[ dat_norm_img2,T2 ] = normalise2dpts(data_orig(4:6,:));
data_norm = [ dat_norm_img1 ; dat_norm_img2 ];

%-----------------------
% Global homography (H).
%-----------------------
% Refine homography using DLT on inliers.
[ ~,A,D1,D2 ] = feval('homography_fit',data_norm(:,:));

%-------------------------
% Moving DLT (projective).
%-------------------------
% Image keypoints coordinates.
Kp = [data_orig(1,:)' data_orig(2,:)'];
% Generating mesh for MDLT.
[ X,Y ] = meshgrid(linspace(1,cw,C1),linspace(1,ch,C2));
% Mesh (cells) vertices' coordinates.
Mv = [X(:)-off(1), Y(:)-off(2)];
% Perform Moving DLT
Hmdlt = zeros(size(Mv,1),9);
parfor i=1:size(Mv,1)    
    % Obtain kernel    
    Gki = exp(-pdist2(Mv(i,:),Kp)./sigma^2);   
    % Capping/offsetting kernel
    Wi = max(gamma,Gki);     
    % This function receives W and A and obtains the least significant 
    % right singular vector of W*A by means of SVD on WA (Weighted SVD).
    v = wsvd(Wi,A);
    h = reshape(v,3,3)';            
    % De-condition
    h = D2\h*D1;
    % De-normalize
    h = T2\h*T1;
    
    Hmdlt(i,:) = h(:);
end
%---------------------------------
% Image stitching with Moving DLT.
%---------------------------------
% Warping images with Moving DLT.
warped_img1 = uint8(zeros(ch,cw,3));
warped_img1(off(2):(off(2)+size(img1,1)-1),off(1):(off(1)+size(img1,2)-1),:) = img1;
[warped_img2] = imagewarping(double(ch),double(cw),double(img2),Hmdlt,double(off),X(1,:),Y(:,1)');
warped_img2 = reshape(uint8(warped_img2),size(warped_img2,1),size(warped_img2,2)/3,3);
% Blending images by averaging (linear blending)
linear_mdlt = imageblending(warped_img1,warped_img2);
end