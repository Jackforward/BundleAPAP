clear;clc;close all;
NUM_PIC = 3;

%-----------------------------
% Set Path and Global Varience
%-----------------------------
gamma = 0.01; % Normalizer for Moving DLT. (0.0015-0.1 are usually good numbers).
sigma = 8.5;  % Bandwidth for Moving DLT. (Between 8-12 are good numbers).   

C1 = 100; % Resolution/grid-size for the mapping function in MDLT (C1 x C2).
C2 = 100;

% addpath('../mdlt/modelspecific');
% addpath('../mdlt/mexfiles');
% addpath('../mdlt/multigs');
% addpath('../mdlt');

I1c = imread('E:\stitching\source\low_res\G9PQ0282.jpg');
I2c = imread('E:\stitching\source\low_res\G9PQ0283.jpg');
I3c = imread('E:\stitching\source\low_res\G9PQ0284.jpg');
I1 = rgb2gray(I1c);I2 = rgb2gray(I2c);I3 = rgb2gray(I3c);

%--------------------------------------
% Find Matching Features Between Images
%--------------------------------------
% Detect features in images.
ptsI1 = detectSURFFeatures(I1, 'MetricThreshold', 500);
ptsI2 = detectSURFFeatures(I2, 'MetricThreshold', 500);
ptsI3 = detectSURFFeatures(I3, 'MetricThreshold', 500);
% Extract feature descriptors.
[featuresI1, valid_ptsI1] = extractFeatures(I1, ptsI1);
[featuresI2, valid_ptsI2] = extractFeatures(I2, ptsI2);
[featuresI3, valid_ptsI3] = extractFeatures(I3, ptsI3);
% Match features by using their descriptors.
match12 = matchFeatures(featuresI1, featuresI2);
match13 = matchFeatures(featuresI1, featuresI3);
match23 = matchFeatures(featuresI2, featuresI3);
% Ransac.
[match12_inliar,H12,A{1,2},T1{1,2},T2{1,2},D1{1,2},D2{1,2}] = ransac(valid_ptsI1, valid_ptsI2, match12);
[match13_inliar,H13,A{1,3},T1{1,3},T2{1,3},D1{1,3},D2{1,3}] = ransac(valid_ptsI1, valid_ptsI3, match13);
[match23_inliar,H23,A{2,3},T1{2,3},T2{2,3},D1{2,3},D2{2,3}] = ransac(valid_ptsI2, valid_ptsI3, match23);
% Merge separate match tables into one.
match = match12_inliar;
match(:, 3) = 0;
match = merge(match, match13_inliar, 1, 3);
match = merge(match, match23_inliar, 2, 3);

%-----------------------------------------
% Map all keypoints to the reference frame (For now reference frame is 2nd pic)
%-----------------------------------------
% Find inliars' coordinates.
valid_pts{1} = valid_ptsI1.Location;
valid_pts{2} = valid_ptsI2.Location;
valid_pts{3} = valid_ptsI3.Location;
for row = 1:size(match,1)
    for col = 1:size(match,2)
        if match(row, col) ~= 0
            keypoint(row, 2 * col - 1:2 * col) = valid_pts{col}(match(row, col), :);
        else
            keypoint(row, 2 * col - 1:2 * col) = 0;
        end
    end
end
% Map keypoints & average their coordinates
point1 = keypoint(:,1:2)';point1(3, :) =(point1(1, :) + point1(2, :)) ~= 0;
point2 = keypoint(:,3:4)';point2(3, :) =(point2(1, :) + point2(2, :)) ~= 0;
point3 = keypoint(:,5:6)';point3(3, :) =(point3(1, :) + point3(2, :)) ~= 0;
point1_map = H12 * point1;point1_map = regularize(point1_map);
point2_map = point2;
point3_map = H23 \ point3;point3_map = regularize(point3_map);
point = point1_map + point2_map + point3_map;point = regularize(point)';
% Obtaining size of canvas (ref frame does not change)
TL1 = regularize(H12 * [1;1;1]);
BL1 = regularize(H12 * [1;size(I1,1);1]);
TR1 = regularize(H12 * [size(I1,2);1;1]);
BR1 = regularize(H12 * [size(I1,2);size(I1,1);1]);
TL2 = [1;1;1];
BL2 = [1;size(I2,1);1];
TR2 = [size(I2,2);1;1];
BR2 = [size(I2,2);size(I2,1);1];
TL3 = regularize(H23 \ [1;1;1]);
BL3 = regularize(H23 \ [1;size(I3,1);1]);
TR3 = regularize(H23 \ [size(I3,2);1;1]);
BR3 = regularize(H23 \ [size(I3,2);size(I3,1);1]);
corners = [TL1 BL1 TR1 BR1 TL2 BL2 TR2 BR2 TL3 BL3 TR3 BR3]; %TL3 BL3 TR3 BR3
cw = round(max(corners(1,:)) - min(corners(1,:)) + 1);
ch = round(max(corners(2,:)) - min(corners(2,:)) + 1);
off = round([ 1 - min(corners(1,:)) + 1 ; 1 - min(corners(2,:)) + 1 ]);

%------------------------------------
% As-Projective-As-Possible Stitching
%------------------------------------
% % Stitching 1st and 2nd pictures.
% clear data1 data2 data_orig pos
% pos = keypoint(:,1) + keypoint(:,2) ~= 0;
% data1 = point(pos, 1:2);
% data2 = keypoint(pos, 1:2);
% data_orig(1:2, :) = data1';data_orig(3, :) = 1;
% data_orig(4:5, :) = data2';data_orig(6, :) = 1;
% data_orig = double(data_orig);
% [result1, Hmdlt12] = stitching(I2c, I1c, data_orig, off, cw, ch);
% imshow(result1);
% 
% % Add 3rd picthre to pano.
% clear data1 data2 data_orig pos
% pos = keypoint(:,5) + keypoint(:,6) ~= 0;
% data1 = point(pos, 1:2);
% % data1(:,1) = data1(:,1) + prev_off(1);
% % data1(:,2) = data1(:,2) + prev_off(2);
% data2 = keypoint(pos, 5:6);
% data_orig(1:2, :) = data1';data_orig(3, :) = 1;
% data_orig(4:5, :) = data2';data_orig(6, :) = 1;
% data_orig = double(data_orig);
% [result3, Hmdlt32] = stitching(I2c, I3c, data_orig, off, cw, ch);
% figure;imshow(result3);


% %------------------
% % Bundle Adjustment
% %------------------
img1 = I2c;
img2 = I1c;

% Convert all point data to double
point = double(point);
keypoint = double(keypoint);
data_orig =  point';
% Prepare keypoint match data
clear pos                                              % pair between 1st and ref
pos = keypoint(:,1) + keypoint(:,2) ~= 0;
match_pair_1 = [point(pos,:)';keypoint(pos, 1:2)'];match_pair_1(6,:) = 1;
% Image keypoints coordinates.
Kp = [data_orig(1,:)' data_orig(2,:)'];
% Generating mesh for MDLT.
[ X,Y ] = meshgrid(linspace(1,cw,C1),linspace(1,ch,C2));
BW = X(1,2) - X(1,1); BH = Y(2,1) - Y(1,1); % Block width & height
% Mesh (cells) vertices' coordinates.
Mv = [X(:)-off(1), Y(:)-off(2)];
% Perform Moving DLT
Hmdlt = zeros(size(Mv,1),9);
%err = 0;werr = 0;
%test = [];
parfor i=1:size(Mv,1)    
    % Obtain kernel    
    Gki = exp(-pdist2(Mv(i,:),Kp)./sigma^2);   
    % Capping/offsetting kernel
    Wi = max(gamma,Gki);     
    % This function receives W and A and obtains the least significant 
    % right singular vector of W*A by means of SVD on WA (Weighted SVD).
    
    % picture wise mdlt
    % 1st
    v = wsvd(Wi,A{1,2});
    h = reshape(v,3,3)';            
    % De-condition
    d1 = D1{1,2}; 
    d2 = D2{1,2}; 
    t1 = T1{1,2}; 
    t2 = T2{1,2};
    h = d2\h*d1;
    h = t2\h*t1;
    h = inv(h);

    
    % Compute err.
    pos = Mv(i,1) < match_pair_1(1,:)  & match_pair_1(1,:)  < Mv(i,1) + BW & Mv(i,2) < match_pair_1(2,:)  & match_pair_1(2,:) < Mv(i,2) + BH;
    center = [Mv(i,1) + BW /2; Mv(i,2) + BH /2];
    match_up_ref = match_pair_1(1:3,pos); match_up_src = match_pair_1(4:6,pos);
%     mapped_ref = regularize(h * match_up_ref);
%     err = err + dist(match_up_src,mapped_ref);
    %test = [test;sum(pos)];
    
%     Gki_star = exp(-pdist2(center',match_up_ref(1:2,:)')./sigma^2);  
%     W_star = max(gamma,Gki_star); 
%     werr = werr + wdist(match_up_src,mapped_ref,W_star);
    option = optimset('MaxFunEvals',500);
    if (mdlt_cost(center, h, match_up_src, match_up_ref,sigma,gamma) ~= 0)
        h = fminsearch(@(h)mdlt_cost(center, h, match_up_src, match_up_ref,sigma,gamma),h,option);
    end
    Hmdlt(i,:) = h(:);

end

% ---------------------------------
% Image stitching with Moving DLT.
% ---------------------------------
% Warping images with Moving DLT.
warped_img1 = uint8(zeros(ch,cw,3));
warped_img1(off(2):(off(2)+size(img1,1)-1),off(1):(off(1)+size(img1,2)-1),:) = img1;
[warped_img2] = imagewarping(double(ch),double(cw),double(img2),Hmdlt,double(off),X(1,:),Y(:,1)');
warped_img2 = reshape(uint8(warped_img2),size(warped_img2,1),size(warped_img2,2)/3,3);
% Blending images by averaging (linear blending)
linear_mdlt = imageblending(warped_img1,warped_img2);
imshow(linear_mdlt);