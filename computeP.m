% left Reference Image - e.g. depth map is from the left image frame
function [p, h] = computeP(depthMap, params)

    [rows, cols] = size(depthMap);

    pixel_l = 1:1:cols;
    pizel_l_delta = (-pixel_l+cols/2);
    x_l = pizel_l_delta*.006;
    X_l = x_l.*depthMap/params.f_l;
    X_r = X_l + params.b;
    x_r = params.f_r*X_r./(depthMap-params.d);
    pixel_r = -x_r/.006 + cols/2;
    p = (params.f_r/params.f_l)*(depthMap)./(depthMap-params.d);
    h = pixel_l-pixel_r;
    

end