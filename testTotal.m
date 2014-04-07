clear all

minNumOfPixelsHeight = 30; %the minimum number of pixels in the highest level
alpha = 1e-13;

im1 = imread('testVideo\3\frame00000.jpg');
im2 = imread('testVideo\3\frame00005.jpg');

imgray1 = rgb2gray(im1);
imgray1 = imresize(imgray1, 0.5);
imgray2 = rgb2gray(im2);
imgray2 = imresize(imgray2, 0.5);

% imgray1 = medfilt2(imgray1, [3, 3]);
% imgray2 = medfilt2(imgray2, [3, 3]);

%Calculate the height, width and number of level that can make them fits the quadtree
%structure
[height, width] = size(imgray1);
scale = height/minNumOfPixelsHeight;
nlevel = floor(log(scale)/log(2));
scale = 2^nlevel;
minNumOfPixelsHeight = floor(height/scale);
minNumOfPixelsWidth = floor(width/scale);
newHeight = scale * minNumOfPixelsHeight;
newWidth = scale * minNumOfPixelsWidth;
rimgray1 = imresize(imgray1, [newHeight, newWidth]);
rimgray2 = imresize(imgray2, [newHeight, newWidth]);

%%level 0
currentHeight = minNumOfPixelsHeight;
currentWidth = minNumOfPixelsWidth;
%the learning rate changes related to size of the image
lr = alpha/(currentHeight*currentWidth);
img1 = imresize(rimgray1, [currentHeight, currentWidth]);
img2 = imresize(rimgray2, [currentHeight, currentWidth]);
img1 = uint8(img1);
img2 = uint8(img2);
[A, error, ~] = findFlow(uint8(img1), uint8(img2), lr, 2^nlevel);
fprintf('Level: 0 finished \n');

for i = 1:nlevel
    %double the size of image
    currentHeight = currentHeight * 2;
    currentWidth = currentWidth * 2;
    %the learning rate changes related to size of the image
    lr = alpha/(currentHeight*currentWidth);
    img1 = imresize(img1, [currentHeight, currentWidth]);
    img2 = imresize(rimgray2, [currentHeight, currentWidth]);
    img1 = uint8(img1);
    img2 = uint8(img2);
    [A, error, ~] = findFlow(double(img1), double(img2), lr, 2^(nlevel-i), [A(1, 1); A(1, 2); A(1, 3); A(2, 1); A(2, 2); A(2, 3)]);
    fprintf('Level: %d finished \n', i);
end

%construct the estimated image based on the frame 1 and the transformation matrix
newPoints = zeros(currentHeight, currentWidth, 2);
showimg1 = imresize(rimgray1, [currentHeight, currentWidth]);
showimg2 = imresize(rimgray2, [currentHeight, currentWidth]);
newimg = zeros(currentHeight, currentWidth);
for i = 1:currentHeight
    for j = 1:currentWidth
        %Get where a point in transormed image coming from frame 1
        fromxy = round(inv(A(1:2, 1:2)) * ([j; i]-A(1:2, 3)));
        if fromxy(2) > height || fromxy(2) <= 0 || fromxy(1) > width || fromxy(1) <= 0
            newimg(i, j) = 0;
        else
            newimg(i, j) = showimg1(fromxy(2), fromxy(1));
        end
    end
end

figure;imshow(showimg1, [])
figure;imshow(newimg, [])
figure;imshow(showimg2, [])
figure;imagesc(abs(newimg-double(showimg2)))
figure;imagesc(abs(newimg-double(showimg1)))
figure;imagesc(abs(double(showimg2)-double(showimg1)))

