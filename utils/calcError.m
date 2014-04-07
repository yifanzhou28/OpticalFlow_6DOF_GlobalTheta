function e = calcError(im1, im2, p1, p2)
[m, n] = size(im1);
if p2(2) > m || p2(1) > n || p2(2) <= 0 || p2(1) <= 0
    e = 0;
else
    e = im1(p1(2), p1(1)) - im2(p2(2), p2(1));
end