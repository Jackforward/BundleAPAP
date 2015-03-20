clear;clc;close all;
NUM_PIC = 4;
REF = 2;
varargin = {'E:\stitching\source\low_res\G9PQ0283.jpg', ...
    'E:\stitching\source\low_res\G9PQ0284.jpg', 'E:\stitching\source\low_res\G9PQ0285.jpg', 'E:\stitching\source\low_res\G9PQ0304.jpg'};
tic;
%-----------------------------
% Set Path and Global Varience
%-----------------------------
gamma = 0.01; % Normalizer for Moving DLT. (0.0015-0.1 are usually good numbers).
sigma = 8.5;  % Bandwidth for Moving DLT. (Between 8-12 are good numbers).   
C1 = 100; % Resolution/grid-size for the mapping function in MDLT (C1 x C2).
C2 = 100;
% 
% addpath('../mdlt/modelspecific');
% addpath('../mdlt/mexfiles');
% addpath('../mdlt/multigs');
% addpath('../mdlt');
% run 'E:\stitching\mdlt\vlfeat-0.9.14\toolbox\vl_setup.m'

Ic = cell(NUM_PIC, 1);
I = cell(NUM_PIC, 1);

for i = 1:NUM_PIC
    Ic{i} = imread(varargin{i});
    I{i} = rgb2gray(Ic{i});
end

%--------------------------------------
% Find Matching Features Between Images
%--------------------------------------
ptsI = cell(NUM_PIC, 1);featuresI = cell(NUM_PIC, 1);valid_ptsI = cell(NUM_PIC, 1);
for i = 1:NUM_PIC
    % Detect features in images.
    ptsI{i} = detectSURFFeatures(I{i}, 'MetricThreshold', 500);

    % Extract feature descriptors.
    [featuresI{i}, valid_ptsI{i}] = extractFeatures(I{i}, ptsI{i});
end

match = cell(NUM_PIC); match_inliar = cell(NUM_PIC);H = cell(NUM_PIC);
A = cell(NUM_PIC); T1 = cell(NUM_PIC); T2 = cell(NUM_PIC); D1 = cell(NUM_PIC); D2 = cell(NUM_PIC);
for i = 1:NUM_PIC
    for j = i+1 : NUM_PIC
        % Match features by using their descriptors.
        match{i,j} = matchFeatures(featuresI{i}, featuresI{j});

        % Ransac.
        [match_inliar{i,j},H{i,j},A{i,j},T1{i,j},T2{i,j},D1{i,j},D2{i,j}] = ransac(valid_ptsI{i}, valid_ptsI{j}, match{i,j});
    end
end

% Merge separate match tables into one.
match = match_inliar{1,2};match(:, 3:NUM_PIC) = 0;
for i = 1:NUM_PIC
    for j = i+1 : NUM_PIC
        match = merge(match, match_inliar{i,j}, i, j);
    end
end


%-----------------------------------------
% Map all keypoints to the reference frame (For now reference frame is 2nd pic)
%-----------------------------------------
% Find inliars' coordinates.
valid_pts = cell(1,NUM_PIC);
for i = 1:NUM_PIC
    valid_pts{i} = valid_ptsI{i}.Location;
end
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
% Meanwhile obtaining size of canvas (ref frame does not change)
point = cell(NUM_PIC, 1); point_map = cell(NUM_PIC, 1);
for i = 1:NUM_PIC
    point{i} = keypoint(:,2*i-1:2*i)';point{i}(3, :) =(point{i}(1, :) + point{i}(2, :)) ~= 0;
end

for i = 1:REF-1
    point_map{i} = H{i,REF} * point{i};point_map{i} = regularize(point_map{i});
    TL{i} = regularize(H{i,REF} * [1;1;1]);
    BL{i} = regularize(H{i,REF} * [1;size(I{i},1);1]);
    TR{i} = regularize(H{i,REF} * [size(I{i},2);1;1]);
    BR{i} = regularize(H{i,REF} * [size(I{i},2);size(I{i},1);1]);
end
point_map{REF} = point{REF};
TL{REF} = [1;1;1];
BL{REF} = [1;size(I{REF},1);1];
TR{REF} = [size(I{REF},2);1;1];
BR{REF} = [size(I{REF},2);size(I{REF},1);1];
for i = REF+1 : NUM_PIC
    point_map{i} = H{REF,i} \ point{i};point_map{i} = regularize(point_map{i});
    TL{i} = regularize(H{REF,i} \ [1;1;1]);
    BL{i} = regularize(H{REF,i} \ [1;size(I{i},1);1]);
    TR{i} = regularize(H{REF,i} \ [size(I{i},2);1;1]);
    BR{i} = regularize(H{REF,i} \ [size(I{i},2);size(I{i},1);1]);
end
clear point
point = point_map{1};
for i = 2:NUM_PIC
    point = point + point_map{i};
end
point = regularize(point)';

% Obtaining size of canvas (ref frame does not change)
corners = [];
for i = 1:NUM_PIC
    corners = [corners TL{i} BL{i} TR{i} BR{i}];
end
cw = round(max(corners(1,:)) - min(corners(1,:)) + 1);
ch = round(max(corners(2,:)) - min(corners(2,:)) + 1);
off = round([ 1 - min(corners(1,:)) + 1 ; 1 - min(corners(2,:)) + 1 ]);

% %------------------
% % Bundle Adjustment
% %------------------
% Convert all point data to double
point = double(point);
keypoint = double(keypoint);
% Image keypoints coordinates.
Kp = point(:,1:2);
% Generating mesh for MDLT.
[ X,Y ] = meshgrid(linspace(1,cw,C1),linspace(1,ch,C2));
BW = X(1,2) - X(1,1); BH = Y(2,1) - Y(1,1); % Block width & height
% Mesh (cells) vertices' coordinates.
Mv = [X(:)-off(1), Y(:)-off(2)];
% Perform Moving DLT
% Hmdlt = zeros(size(Mv,1),9);
Hmdlt = cell(1, NUM_PIC);
werr_orig = 0;
werr = 0;
% test = [];
toc;tic;

clear H; H = cell(1,size(Mv,1));pt_src = cell(1,size(Mv,1));
parfor i=  1:size(Mv,1)    
    % Obtain kernel    
    Gki = exp(-pdist2(Mv(i,:),Kp)./sigma^2);   
    % Capping/offsetting kernel
    Wi = max(gamma,Gki);     
    % This function receives W and A and obtains the least significant 
    % right singular vector of W*A by means of SVD on WA (Weighted SVD).
    
    for j = 1 : NUM_PIC
        if j == REF
            H{i}(1:3,3*j-2:3*j) = diag([1,1,1]);
        elseif j < REF
                v = wsvd(Wi,A{j,REF});
                h = reshape(v,3,3)';            
                % De-condition
                d1 = D1{j,REF};d2 = D2{j,REF}; t1 = T1{j,REF}; t2 = T2{j,REF};
                h = d2\h*d1; h = t2\h*t1; h = inv(h);    
                H{i}(1:3,3*j-2:3*j) = h;
        elseif j > REF
                v = wsvd(Wi,A{REF,j});
                h = reshape(v,3,3)';            
                % De-condition
                d1 = D1{REF,j}; d2 = D2{REF,j}; t1 = T1{REF,j}; t2 = T2{REF,j};
                h = d2\h*d1; h = t2\h*t1;
                H{i}(1:3,3*j-2:3*j) = h;
        end
    end
    
    
%     % For compare intention
%     for j = 1:NUM_PIC
%         Hmdlt_orig{i}(j,:) = reshape(H{i}(1:3,3*j-2:3*j),1,[]);
%     end
%     
    
    % Bundle Adjustment
    pos = Mv(i,1) < point(:,1)  & point(:,1)  < Mv(i,1) + BW...
        & Mv(i,2) < point(:,2)  & point(:,2) < Mv(i,2) + BH;
    pt_ref = point(pos,:);
    if ~isempty(pt_ref)
        match_table = [];
        center = [Mv(i,1) + BW /2; Mv(i,2) + BH /2];
        for j = 1 : NUM_PIC
            pt_src{i}(:,j*3-2:j*3-1) = keypoint(pos,j*2-1:j*2); 
            pt_src{i}(:,j*3) = (pt_src{i}(:,j*3-2) + pt_src{i}(:,j*3-1)) ~= 0;
            match_table = [match_table H{i}(1:3,3*j-2:3*j)];
        end
        match_table = [match_table pt_ref'];
        pt_src_ = pt_src{i};
        werr_orig = werr_orig + par_bundle_cost(match_table, pt_src_, center, NUM_PIC, sigma, gamma);
        option = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
        match_table = lsqnonlin(@(match_table)par_bundle_cost(match_table, pt_src_, center, NUM_PIC, sigma, gamma),match_table,[],[],option);
        for j = 1:NUM_PIC
            H{i}(1:3,3*j-2:3*j) = match_table(1:3,3*j-2:3*j);
        end
        werr = werr + par_bundle_cost(match_table, pt_src_, center, NUM_PIC, sigma, gamma);

    end
    
    for j = 1:NUM_PIC
        Hmdlt{i}(j,:) = reshape(H{i}(1:3,3*j-2:3*j),1,[]);
    end
        
end
toc;tic;

% % For compare intention
% for j = 1:NUM_PIC
%     for i =  1:size(Mv,1)
%         Hmdlt_{j}(i,:) = Hmdlt_orig{i}(j,:);
%     end
% end


for j = 1:NUM_PIC
    for i =  1:size(Mv,1)
        Hmdlt_{j}(i,:) = Hmdlt{i}(j,:);
    end
end


% ---------------------------------
% Image stitching with Moving DLT.
% ---------------------------------
% Warping images with Moving DLT.
linear_mdlt = uint8(zeros(ch,cw,3));
linear_mdlt(off(2):(off(2)+size(Ic{REF},1)-1),off(1):(off(1)+size(Ic{REF},2)-1),:) = Ic{REF};
for i = 1:NUM_PIC
    [warped_img] = imagewarping(double(ch),double(cw),double(Ic{i}),Hmdlt_{i},double(off),X(1,:),Y(:,1)');
    warped_img = reshape(uint8(warped_img),size(warped_img,1),size(warped_img,2)/3,3);
    linear_mdlt = imageblending(linear_mdlt,warped_img);
end
imshow(linear_mdlt);