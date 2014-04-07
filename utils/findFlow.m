function [ A, new_error, sign] = findFlow(img1, img2, lr0, scale, initA)
%FINDFLOW Summary of this function goes here
%   Detailed explanation goes here
numIter = 3;

if nargin == 4
    %Initialize the matrix
    a11 = 1;
    a12 = 0;
    a13 = 0;
    a21 = 0;
    a22 = 1;
    a23 = 0;
elseif nargin == 5
    %take the matrix from previous level
    a11 = initA(1, 1);
    a12 = initA(2, 1);
    a13 = initA(3, 1);
    a21 = initA(4, 1);
    a22 = initA(5, 1);
    a23 = initA(6, 1);
else
    return;
end

[height, width] = size(img1);
Dim = height*width;
%Get the image gradient respect to x and y axis
[Gx, Gy] = calcGradient(img2);
%Adjust the gradient to fit the scale
Gx = Gx./scale;
Gy = Gy./scale;

%%Calculate the errors:
%get p1, p1 is the point in frame 1
p1 = zeros(2, Dim);
for j = 1:height
    for i = 1:width
        p1(:, (j-1)*width+i) = [i, j];
    end
end

%get p2, p2 is the point in the transformed image
p2 = [a11, a12, a13; a21, a22, a23] * [((p1.*2-1).*scale+1)/2; ones(1, Dim)];
p2 = uint8(floor((p2-1)./scale)+1);

%check if all the point in p2 locates in the image, ValidInx comprises all
%the valid points' index
[sign, ValidInx] = checkPoints(p2, height, width);

%get errors
errors = zeros(1, Dim);
for i = 1:size(p1, 2)
    %if the point locates out of the image, error = 0
    if ~isempty(find(i == ValidInx, 1))
        errors(1, i) = calcError(img1, img2, p1(:, i), p2(:, i));
    else
        errors(1, i) = 0;
    end
end

new_error = errors * errors';
fprintf('-----------------------------------------------------------------------\n');
fprintf('iteration : 0   \n');
fprintf('matrix A : \n');
fprintf('%d \t %d \t %d \t \n', a11, a12, a13);
fprintf('%d \t %d \t %d \t \n', a21, a22, a23);
fprintf('new error is : %d \n', new_error);
fprintf('-----------------------------------------------------------------------\n');

%%Gradient Descent Starts
count = 1;
while count <= numIter
    %learning rate changes dynamically
    lr = (lr0/numIter)*(numIter-count+1);
    %The gradient repect to each factors in transformation matrix
    Ga11 = zeros(Dim, 1);
    Ga12 = zeros(Dim, 1);
    Ga13 = zeros(Dim, 1);
    Ga21 = zeros(Dim, 1);
    Ga22 = zeros(Dim, 1);
    Ga23 = zeros(Dim, 1);
    old_error = new_error;
    
    for i = 1:size(p2, 2)
        new_point = p2(:, i);
        newX = round(new_point(1));
        newY = round(new_point(2));
        %if a point trasformed out of the image, its gradient is 0
        if isempty(find(i == ValidInx, 1))
            Ga11(i, 1) = 0;
            Ga12(i, 1) = 0;
            Ga13(i, 1) = 0;
            Ga21(i, 1) = 0;
            Ga22(i, 1) = 0;
            Ga23(i, 1) = 0;
        else
            %Calculate the gradient of each pixels respect to every factors
            %in transformation matrix
            Ga11(i, 1) = -2*errors(1, i) .* ((Gx(newY, newX) .* ((double(newX)*2-1)*scale+1)/2));
            Ga12(i, 1) = -2*errors(1, i) .* ((Gx(newY, newX) .* ((double(newY)*2-1)*scale+1)/2));
            Ga13(i, 1) = -2*errors(1, i) .* (Gx(newY, newX));
            Ga21(i, 1) = -2*errors(1, i) .* ((Gy(newY, newX) .* ((double(newX)*2-1)*scale+1)/2));
            Ga22(i, 1) = -2*errors(1, i) .* ((Gy(newY, newX) .* ((double(newY)*2-1)*scale+1)/2));
            Ga23(i, 1) = -2*errors(1, i) .* (Gy(newY, newX));
        end
    end
    %sum the gradients, corresponding to the integral to the euqation
    Ga11 = sum(Ga11(:));
    Ga12 = sum(Ga12(:));
    Ga13 = sum(Ga13(:));
    Ga21 = sum(Ga21(:));
    Ga22 = sum(Ga22(:));
    Ga23 = sum(Ga23(:));
    
    %gradient descent:
    JG = [Ga11, Ga12, Ga13, Ga21, Ga22, Ga23];
    GF = JG'*old_error';
    
    tempGradient = lr*GF(1);
    a11 = a11 - tempGradient;
    
    tempGradient = lr*GF(2);
    a12 = a12 - tempGradient;
    
    tempGradient = lr*GF(3);
    a13 = a13 - tempGradient;
    
    tempGradient = lr*GF(4);
    a21 = a21 - tempGradient;
    
    tempGradient = lr*GF(5);
    a22 = a22 - tempGradient;
    
    tempGradient = lr*GF(6);
    a23 = a23 - tempGradient;
    
    %get new p2
    p2 = [a11, a12, a13; a21, a22, a23] * [((p1.*2-1).*scale+1)/2; ones(1, Dim)];
    p2 = uint8(floor((p2-1)./scale)+1);
    [sign, ValidInx] = checkPoints(p2, height, width);
    
    %get new error
    for i = 1:size(p1, 2)
        if ~isempty(find(i == ValidInx, 1))
            errors(1, i) = calcError(img1, img2, p1(:, i), p2(:, i));
        else
            errors(1, i) = 0;
        end
    end
    
    new_error = errors*errors';
    fprintf('-----------------------------------------------------------------------\n');
    fprintf('iteration : %d   \n', count);
    fprintf('matrix A : \n');
    fprintf('%d \t %d \t %d \t \n', a11, a12, a13);
    fprintf('%d \t %d \t %d \t \n', a21, a22, a23);
    fprintf('new error is : %d \n', new_error);
    fprintf('-----------------------------------------------------------------------\n');
    count = count + 1;
end
%fprintf('--------------------------------------------------------------------- \n');
A = [a11, a12, a13; a21, a22, a23];

end

