function [ data_orig ] = registrator(point, keypoint, pic_no)
pos = keypoint(:,pic_no *2 -1) + keypoint(:,pic_no *2) ~= 0;
data1 = point(pos, 1:2);
data2 = keypoint(pos, pic_no *2 -1:pic_no *2);
data_orig(1:2, :) = data1';data_orig(3, :) = 1;
data_orig(4:5, :) = data2';data_orig(6, :) = 1;
data_orig = double(data_orig);
end

