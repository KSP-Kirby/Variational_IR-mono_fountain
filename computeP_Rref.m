% right Reference Image - e.g. depth map is from the right image frame
function [p, h] = computeP_Rref(depthMap, params)

    [rows, cols] = size(depthMap);

    pixel_r = 1:1:cols;
    pizel_r_delta = (-pixel_r+cols/2);
    x_r = pizel_r_delta*.006;
    X_r = x_r.*depthMap/params.f_r;
    X_l = X_r - params.b;
    x_l = params.f_l*X_l./(depthMap+params.d);
    pixel_l = -x_l/.006 + cols/2;
    p = (params.f_r/params.f_l)*(depthMap+params.d)./(depthMap);
    h = pixel_l-pixel_r;
    

end

% for pixel_r = 1:640
%     
%     if pixel_r == 300
%         k = 5;
%     end
% 
%     pizel_r_delta = (-pixel_r+cols/2);       % delta from center
%     x_r = pizel_r_delta*.006;
%     X_r = x_r*(depthMap_r(pixel_r))/params.f_r;
% 
%     X_l = X_r - params.b;
%     x_l = params.f_l*X_l./(depthMap_r(pixel_r)+params.d);
% 
%     pixel_l(pixel_r) = -x_l/.006 + cols/2;
% 
%     %m = (params.f_r/params.f_l)*(depthMap)./(depthMap-params.d);
% 
%     m = (params.f_r/params.f_l)*(depthMap_r(pixel_r)+params.d)./(depthMap_r(pixel_r));
%     vl_adj(pixel_r) = vl(pixel_r)/m;
% end
