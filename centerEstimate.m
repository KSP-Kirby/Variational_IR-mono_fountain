function [zEst0, deltaXest] = centerEstimate(w0_left, w0_right, params )
% version 15 October 2015
% this only deals with the first estimate where there is no delta Z
    
    tic
    columns = length(w0_left);
    kernelWidth = 6;            % for front camera, TODO: put this in params

    kernel=ones(1,kernelWidth)/kernelWidth;
    
    % Add padding of one smoothing kernel width
    w0_left = [ones(1,kernelWidth)*w0_left(1),w0_left,ones(1,kernelWidth)*w0_left(end)];
    w0_right = [ones(1,kernelWidth)*w0_right(1),w0_right,ones(1,kernelWidth)*w0_right(end)];
    
    %smooth
    w0_left=conv(w0_left,kernel,'same'); 
    w0_right=conv(w0_right,kernel,'same'); 
    
    %remove padding
    w0_left = w0_left(kernelWidth+1:end-kernelWidth);
    w0_right = w0_right(kernelWidth+1:end-kernelWidth);
    
    [numRows, numCols] = size(w0_left);
    centerCol = round(numCols/2);
    centerRow = round(numRows/2);
    if w0_left(centerCol) < 0
        w0_left = -w0_left;
        w0_right = -w0_right;
    end
    
    % find initial offset
%     plot(w0_left)
%     hold all
%     plot(w0_right)
    
    for Z = params.minZ:10:params.maxZ
        m = (params.f_r/params.f_l)*(Z+params.d)/(Z);
        rightPixel_1 = centerCol;
        x_r = (rightPixel_1-numCols/2)*.006;
        x_l = (-x_r*Z*params.f_l+params.b*params.f_r*params.f_l)/(params.f_r*Z+params.f_r*params.d);
        leftPixel_1 = numCols/2-x_l/.006;        % this is a fractional pixel

        flowLeft = m*w0_right(centerRow, rightPixel_1);
%         plot(round(leftPixel_1),flowLeft,'*')
        
        crossingSign = sign(flowLeft - w0_left(centerRow,round(leftPixel_1)));
        if Z ~= params.minZ
            if lastCrossingSign ~= crossingSign;
                break;
            end  
        end
        lastCrossingSign = crossingSign;
            
%        title(num2str(crossingSign))
  
    end

    deltaXest = Z*flowLeft/params.f_l;
    
    zEst0 = (deltaXest*params.f_l)./w0_left;


end

