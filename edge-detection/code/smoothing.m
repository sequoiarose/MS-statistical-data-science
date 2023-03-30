I = rgb2gray(imread("C:\Users\sequo\OneDrive\Desktop\PXL_20220904_184602699.MP.jpg"));
Iblur = imgaussfilt(I,6);
montage({I,Iblur})
title('Original Image (Left) Vs. Gaussian Filtered Image (Right)')