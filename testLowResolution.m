clear all
%close all

%parameters
minNumOfCellHeight = 5; %the minimum number of cells in a column
range = 10; %the range of pixels in a cell
nHeight = range*2+1;
nWidth = range*2+1;
alpha = 1e-14;

%%read image 1
im1 = imread('testVideo\frame00000.jpg');
if size(im1, 3) ~= 3
    rimgray1 = im1;
else
    rimgray1 = rgb2gray(im1);
end
rimgray1 = cv.resize(rimgray1, 0.8);
rimgray1 = medfilt2(rimgray1, [3, 3]);
%rimgray1 = edge(rimgray1);
%rimgray1 = imgradient(rimgray1);

%%read image 2
im2 = imread('testVideo\frame00005.jpg');
if size(im2, 3) ~= 3
    rimgray2 = im2;
else
    rimgray2 = rgb2gray(im2);
end
rimgray2 = cv.resize(rimgray2, 0.8);
rimgray2 = medfilt2(rimgray2, [3, 3]);
%rimgray2 = edge(rimgray2);
%rimgray2 = imgradient(rimgray2);
rimgray1 = double(rimgray1);
rimgray2 = double(rimgray2);


%%get image information
[height, width] = size(rimgray1);
maxNumOfCellHeight = floor(height/nHeight);
maxScale = floor(maxNumOfCellHeight/minNumOfCellHeight);
nlevel = floor(log(maxScale)/log(2));
maxScale = 2^nlevel;
maxNumOfCellHeight = maxScale * minNumOfCellHeight;
minNumOfCellWidth = floor(width/(height/minNumOfCellHeight));
maxNumOfCellWidth = maxScale * minNumOfCellWidth;

currentNumOfCellHeight = minNumOfCellHeight;
currentNumOfCellWidth = minNumOfCellWidth;

AInLevel = cell(1, nlevel+1);
errors = cell(1, nlevel+1);

%%LEVEL 0
m = currentNumOfCellHeight*nHeight;
n = currentNumOfCellWidth*nWidth;
imgray1 = imresize(rimgray1, [m, n]);
imgray2 = imresize(rimgray2, [m, n]);
diffLevel = nlevel;
scale = 2^diffLevel;
[Gx, Gy] = calcGradient(imgray2);
Gx = Gx ./ scale;
Gy = Gy ./ scale;
for j = range+1:nHeight:m-range
    for i = range+1:nWidth:n-range
        [A, oldc, newc, errors{1, 1}(:, end+1)] = findFlow(imgray1, imgray2, Gx, Gy, [i; j], range, alpha, scale);
        AInLevel{1, 1}(:, end+1) = [i; j; newc(1, 1); newc(2, 1); A(1, 1); A(1, 2); A(1, 3); A(2, 1); A(2, 2); A(2, 3)];
    end
end
fprintf('0 level finished!   \n');

for clevel = 1:1:nlevel
    %
    diffLevel = nlevel - clevel;
    cell_h_prev = currentNumOfCellHeight;
    cell_w_prev = currentNumOfCellWidth;
    currentNumOfCellHeight = currentNumOfCellHeight*2;
    currentNumOfCellWidth = currentNumOfCellWidth*2;
    m = currentNumOfCellHeight*nHeight;
    n = currentNumOfCellWidth*nWidth;
    imgray1 = imresize(rimgray1, [m, n]);
    imgray2 = imresize(rimgray2, [m, n]);
    scale = 2^diffLevel;
    [Gx, Gy] = calcGradient(imgray2);
    Gx = Gx ./ scale;
    Gy = Gy ./ scale;
    AinxH = 1;
    
    for j = range+1:nHeight*2:m-range
        tempAH = AInLevel{clevel}(:, AinxH:AinxH+cell_w_prev-1);
        AinxH = AinxH + cell_w_prev;
        AinxW = 1;
        for i = range+1:nWidth*2:n-range
%             if (i == 41 && j == 17) || (i == 41-nWidth && j == 17)
%                aaaa = 0; 
%             end
            tempAW = tempAH(:, AinxW);
            AinxW = AinxW + 1;
            [A, ~, newc, errors{1, clevel+1}(:, end+1)] = findFlow(imgray1, imgray2, Gx, Gy, [i; j], range, alpha, scale, [tempAW(5, 1);tempAW(6, 1);tempAW(7, 1);tempAW(8, 1);tempAW(9, 1);tempAW(10, 1)]);
            AInLevel{1, clevel+1}(:, end+1) = [i; j; newc(1, 1); newc(2, 1); A(1, 1); A(1, 2); A(1, 3); A(2, 1); A(2, 2); A(2, 3)];
            
            [A, ~, newc, errors{1, clevel+1}(:, end+1)] = findFlow(imgray1, imgray2, Gx, Gy, [i+nWidth; j], range, alpha, scale, [tempAW(5, 1);tempAW(6, 1);tempAW(7, 1);tempAW(8, 1);tempAW(9, 1);tempAW(10, 1)]);
            AInLevel{1, clevel+1}(:, end+1) = [i+nWidth; j; newc(1, 1); newc(2, 1); A(1, 1); A(1, 2); A(1, 3); A(2, 1); A(2, 2); A(2, 3)];
        end
        
        AinxW = 1;
        for i = range+1:nWidth*2:n-range
%             if (i == 41 && j == 17-nHeight) || (i == 41-nWidth && j == 17-nHeight)
%                aaaa = 0; 
%             end
            tempAW = tempAH(:, AinxW);
            AinxW = AinxW + 1;
            [A, ~, newc, errors{1, clevel+1}(:, end+1)] = findFlow(imgray1, imgray2, Gx, Gy, [i; j+nHeight], range, alpha, scale, [tempAW(5, 1);tempAW(6, 1);tempAW(7, 1);tempAW(8, 1);tempAW(9, 1);tempAW(10, 1)]);
            AInLevel{1, clevel+1}(:, end+1) = [i; j+nHeight; newc(1, 1); newc(2, 1); A(1, 1); A(1, 2); A(1, 3); A(2, 1); A(2, 2); A(2, 3)];
            
            [A, ~, newc, errors{1, clevel+1}(:, end+1)] = findFlow(imgray1, imgray2, Gx, Gy, [i+nWidth; j+nHeight], range, alpha, scale, [tempAW(5, 1);tempAW(6, 1);tempAW(7, 1);tempAW(8, 1);tempAW(9, 1);tempAW(10, 1)]);
            AInLevel{1, clevel+1}(:, end+1) = [i+nWidth; j+nHeight; newc(1, 1); newc(2, 1); A(1, 1); A(1, 2); A(1, 3); A(2, 1); A(2, 2); A(2, 3)];
        end
    end
    
    fprintf('%d level finished!   \n', clevel);
end

flow = zeros(m, n, 2);
temp = AInLevel{1, end};
for i = 1:size(temp, 2)
    x = temp(1, i);
    y = temp(2, i);
    A = [temp(5, i), temp(6, i), temp(7, i); temp(8, i), temp(9, i), temp(10, i)];
    for j = -range:1:range
        for k = -range:1:range
            c = [x+j; y+k];
            newc = A * [c; 1];
            
            u = newc(1, 1) - c(1, 1);
            flow(c(2, 1), c(1, 1), 1) = u;
            
            v = newc(2, 1) - c(2, 1);
            flow(c(2, 1), c(1, 1), 2) = v;
        end
    end
    
end
myflowcolour = flowToColor(flow);
figure;imshow(myflowcolour);

im11 = imresize(imgray1, [minNumOfCellHeight*nHeight*2^5, minNumOfCellWidth*nWidth*2^5]);
im22 = imresize(imgray2, [minNumOfCellHeight*nHeight*2^5, minNumOfCellWidth*nWidth*2^5]);
test = AInLevel{1, 6};
flow = zeros(minNumOfCellHeight*nHeight*2^5, minNumOfCellWidth*nWidth*2^5, 2);
figure;imshow(im11+im22, []);hold on;
for i = 1:size(test, 2)
    %quiver(test(1, i),test(2, i),test(3, i)-test(1, i),test(4, i)-test(2, i), 'color', 'r');
    flow(test(2, i), test(1, i), 1) = test(4, i)-test(2, i);
    flow(test(2, i), test(1, i), 2) = test(3, i)-test(1, i);
end
myflowcolour = flowToColor(flow);
figure;imshow(myflowcolour)

im11 = imresize(imgray1, [minNumOfCellHeight*nHeight*2^0, minNumOfCellWidth*nWidth*2^0]);
im22 = imresize(imgray2, [minNumOfCellHeight*nHeight*2^0, minNumOfCellWidth*nWidth*2^0]);
test = AInLevel{1, 1};
figure;subplot(1, 2, 1);imshow(im11, []);hold on;
for i = 1:size(test, 2)
    quiver(test(1, i),test(2, i),test(3, i)-test(1, i),test(4, i)-test(2, i), 'color', 'r');
end
subplot(1, 2, 2);imshow(im22, []);

im11 = imresize(imgray1, [minNumOfCellHeight*nHeight*2^1, minNumOfCellWidth*nWidth*2^1]);
im22 = imresize(imgray2, [minNumOfCellHeight*nHeight*2^1, minNumOfCellWidth*nWidth*2^1]);
test = AInLevel{1, 2};
figure;subplot(1, 2, 1);imshow(im11, []);hold on;
for i = 1:size(test, 2)
    quiver(test(1, i),test(2, i),test(3, i)-test(1, i),test(4, i)-test(2, i), 'color', 'r');
end
subplot(1, 2, 2);imshow(im22, []);

im11 = imresize(imgray1, [minNumOfCellHeight*nHeight*2^2, minNumOfCellWidth*nWidth*2^2]);
im22 = imresize(imgray2, [minNumOfCellHeight*nHeight*2^2, minNumOfCellWidth*nWidth*2^2]);
test = AInLevel{1, 3};
figure;subplot(1, 2, 1);imshow(im11, []);hold on;
for i = 1:size(test, 2)
    quiver(test(1, i),test(2, i),test(3, i)-test(1, i),test(4, i)-test(2, i), 'color', 'r');
end
subplot(1, 2, 2);imshow(im22, []);

im11 = imresize(imgray1, [minNumOfCellHeight*nHeight*2^3, minNumOfCellWidth*nWidth*2^3]);
im22 = imresize(imgray2, [minNumOfCellHeight*nHeight*2^3, minNumOfCellWidth*nWidth*2^3]);
test = AInLevel{1, 4};
figure;subplot(1, 2, 1);imshow(im11, []);hold on;
for i = 1:size(test, 2)
    quiver(test(1, i),test(2, i),test(3, i)-test(1, i),test(4, i)-test(2, i), 'color', 'r');
end
subplot(1, 2, 2);imshow(im22, []);

im11 = imresize(imgray1, [minNumOfCellHeight*nHeight*2^4, minNumOfCellWidth*nWidth*2^4]);
im22 = imresize(imgray2, [minNumOfCellHeight*nHeight*2^4, minNumOfCellWidth*nWidth*2^4]);
test = AInLevel{1, 5};
figure;subplot(1, 2, 1);imshow(im11, []);hold on;
for i = 1:size(test, 2)
    quiver(test(1, i),test(2, i),test(3, i)-test(1, i),test(4, i)-test(2, i), 'color', 'r');
end
subplot(1, 2, 2);imshow(im22, []);

im11 = imresize(imgray1, [minNumOfCellHeight*nHeight*2^5, minNumOfCellWidth*nWidth*2^5]);
im22 = imresize(imgray2, [minNumOfCellHeight*nHeight*2^5, minNumOfCellWidth*nWidth*2^5]);
test = AInLevel{1, 6};
figure;subplot(1, 2, 1);imshow(im11, []);hold on;
for i = 1:size(test, 2)
    quiver(test(1, i),test(2, i),test(3, i)-test(1, i),test(4, i)-test(2, i), 'color', 'r');
end
subplot(1, 2, 2);imshow(im22, []);

im11 = imresize(imgray1, [minNumOfCellHeight*nHeight*2^6, minNumOfCellWidth*nWidth*2^6]);
im22 = imresize(imgray2, [minNumOfCellHeight*nHeight*2^6, minNumOfCellWidth*nWidth*2^6]);
test = AInLevel{1, 7};
figure;subplot(1, 2, 1);imshow(im11, []);hold on;
for i = 1:size(test, 2)
    quiver(test(1, i),test(2, i),test(3, i)-test(1, i),test(4, i)-test(2, i), 'color', 'r');
end
subplot(1, 2, 2);imshow(im22, []);
