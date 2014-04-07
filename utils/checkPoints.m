function [ sign, validInx, validPoints ] = checkPoints( points, height, width )
%CHECKPOINTS Summary of this function goes here
%   Detailed explanation goes here
sign = 0;
validInx = [];
for i = 1:size(points, 2)
    x = points(1, i);
    y = points(2, i);
    if x <= 0 || y <= 0 || x > width || y > height
       sign = 1;
    else
        validInx = [validInx, i];
    end
end
validPoints = points(:, validInx);

end

