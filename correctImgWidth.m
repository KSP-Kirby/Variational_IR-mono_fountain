function [ imgOut ] = correctImgWidth( img )
%Corrects the image width



img = imcrop(img,[64, 1, 511, 480]);

A = [640/512, 0 , 0 ; 0, 1, 0 ; 0, 0, 1];

tform = affine2d(A);

imgOut = imwarp(img, tform);

imtool(imgOut)


end

