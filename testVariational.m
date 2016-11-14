
% version 14 October 2015
% the function scaleAndCrop converts the different pixel size and fl of the
% IR camera to match that of the VL camera.

% reference image is left image, equation 3.24 in dissertation
% right camera is closer to scene by d than left camera

%params.columns = 400;
tic
load('flow.mat');

params.f_l = 6.1762;        % left camera is ir
params.f_r = 6.1762;         % right camera is RGB
params.pixelDim = .006;
params.b =  3*25.4;         % stereo baseline
params.d =  300;         % dual focal length baseline
params.lambda = 40000;
params.alpha = .01;
params.minZ = 600;
params.maxZ = 3000;
params.iterations1 = 50;
startSeq = 1;

results = [];

[rows, cols, p] = size(uv_vl{startSeq});

vl = -uv_vl{startSeq}(240,:,1);
ir = -uv_ir{startSeq}(240,:,1);

centerRow = round(rows/2);
i = centerRow;
% The center estimate is for left or back camera so it includes d

%[zEst0, deltaXest] = centerEstimate_Lref(uv_vl{startSeq}(i,:,1), uv_ir{startSeq}(i,:,1), params);

[zEst0, deltaXest] = centerEstimate(uv_vl{startSeq}(i,:,1), uv_ir{startSeq}(i,:,1), params);

%deltaXest = 20/.006;
%this comes from the left depth (zEst0) found above times the left flow
disp(strcat('deltaXest:',num2str(deltaXest*.006)))

h = waitbar(0, 'Working')
zEst_l = [];
for i = 1:rows
    
    zEst_l(i,:) = (deltaXest*params.f_l)./(-uv_ir{startSeq}(i,:,1));
    %zEst0(i,:) = centerEstimate(uv_vl{1}(i,:,1), uv_ir{1}(i,:,1), params);
    waitbar(i/rows)
end
close(h)

h = waitbar(0, 'Working')
zEst_r = [];
for i = 1:rows
    
    zEst_r(i,:) = (deltaXest*params.f_r)./(-uv_vl{startSeq}(i,:,1));
    %zEst0(i,:) = centerEstimate(uv_vl{1}(i,:,1), uv_ir{1}(i,:,1), params);
    waitbar(i/rows)
end
close(h)


% plot(zEst_l(240,:)-300)
% hold all
% plot(zEst_r(240,:))


for i = 1:rows
    for j = 1:cols
        if zEst_l(i,j) > params.maxZ
            zEst_l(i,j) = params.maxZ;
        end
        
        if zEst_l(i,j) < params.minZ
            zEst_l(i,j) = params.minZ;
        end    
    end
end

for i = 1:rows
    for j = 1:cols
        if zEst_r(i,j) > params.maxZ
            zEst_r(i,j) = params.maxZ;
        end
        
        if zEst_r(i,j) < params.minZ
            zEst_r(i,j) = params.minZ;
        end    
    end
end

% do a shifted plot to make sure it works

i = 240;
depthMap = zEst_l(i,:);      % depthMap is from the perspective of the left (IR) camera
depthMap_r = zEst_r(i,:);      % depthMap is from the perspective of the left (IR) camera

%[p, h] = computeP_Rref(depthMap_r, params);

% pixel_r = 1:1:640;
% pixel_l = pixel_r + h;
% 
% vl_adj = vl.*p;

% figure
% plot(vl)
% hold all
% plot(ir)

%left image is reference image
pixel_l = 1:640;
[p, h] = computeP(depthMap, params);
% plot(p.*ir)
% plot(pixel_l+h, vl)

% adjust the left flow based on depth to match right flow
rightFlow_adj = [];
for pixel_l = 1:640;
    %plot(pixel_l,ir(pixel_l),'*r')
    pizel_l_delta = (-pixel_l+cols/2);
    x_l = pizel_l_delta*params.pixelDim;
    X_l = x_l.*depthMap(pixel_l)/params.f_l;
    X_r = X_l + params.b;
    x_r = params.f_r*X_r./(depthMap(pixel_l)-params.d);
    pixel_r = -x_r/params.pixelDim + cols/2;
    if round(pixel_r) > 0 && round(pixel_r) <=640
        %plot(pixel_r, vl(round(pixel_r)),'*r')
        p(pixel_l) = (params.f_l/params.f_r)*(depthMap(pixel_l)-params.d)./(depthMap(pixel_l));
        rightFlow_adj(pixel_l,1) = pixel_r;
        rightFlow_adj(pixel_l,2) = vl(round(pixel_r))/p(pixel_l);
        %plot(pixel_r, ir(pixel_l)/p,'*m')
    end
    h1(pixel_l) = pixel_l-pixel_r;
end

iterations = 6;
h = waitbar(0,'Interating')
for i = 1:iterations
    h = waitbar(0,strcat('Iteration:',num2str(i),' of:', num2str(iterations)))
    for j = 1:rows
        [zEst_l(i,:), zEst1_1(i,:),rms(i)] = imrStereo_withGraphics(-uv_ir{startSeq}(i,:,1), -uv_vl{startSeq}(i,:,1), zEst_l(i,:), params );
        waitbar(j/rows)
    end
    close(h)
end

[ imgOut ] = correctImgWidth(zEst_l);