function [point] = regularize(point)
%regularize points in row coordinate
point(1, :) = point(1, :) ./ point(3, :);
point(2, :) = point(2, :) ./ point(3, :);
point(3, :) = point(3, :) ./ point(3, :);
pos_nan = isnan(point);
point(pos_nan) = 0;
end