function [zEst0, zEst1, rms_flow] = imrStereo_withGraphics(w0_left, w0_right, zEst0, params )
% version 15 October 2015
% this only deals with the first estimate where there is no delta Z
    
    columns = length(w0_left);
    kernelWidth = 6;            % for front camera, TODO: put this in params

    kernel=ones(1,kernelWidth)/kernelWidth;
    u = (0:1:columns-1);
    u = (u - round(columns/2))*params.pixelDim;
    b = params.b;               %base line
    
    [p, h] = computeP(zEst0, params);      %comput P and h
    pixel_l = 1:1:columns;
    %pixel_l = pixel_r + h;                      % compute pixel cooredinates for right image using h
    vl_adj = w0_right;                       % adjust by p
    %[ adjX, vl_adj2 ] = fixXaxisReversal_Rref(pixel_l, vl_adj ); % fix reversals

%     figure
%     plot(p.*w0_left)
%     hold all
%     plot(pixel_l + h, w0_right)
    
    pixel_l = 1:640;



    [ right_adj ] = resample2Grid( vl_adj, pixel_l + h);
    
%     plot(right_adj)
    
    for i = 1:length(right_adj);
        if isnan(right_adj(i))
            right_adj(i) = right_adj(1);
        end
    end
    
    
    % h0 = h0*0-(params.fl*params.b)/zEstCenter;
     
    % w0_rightShifted = w0_right;

    %clf
    I0z = right_adj - p.*w0_left;           % this is the error that I'm trying to minimize
    %plot(I0z)
    
    

%     I0z = w0_left-w0_rightShifted;                   % this is the error that I'm trying to minimize
%     w0_rightShifted = [w0_rightShifted(3),w0_rightShifted(2),w0_rightShifted,w0_rightShifted(end-2),w0_rightShifted(end-3)];  % pad
%     I0zx = conv(w0_rightShifted,[-0.5 0 0.5],'same');
%     I0zx = I0zx(3:end-2);
%     w0_rightShifted = w0_rightShifted(3:end-2);

    % find laplacian of Z
%     h0 = [h0(3),h0(2),h0,h0(end-2),h0(end-3)];  % pad
%     dxx_h0 = conv(h0,0.5*[1 -2 1],'same');      % 2nd derivative
%     dxx_h0 = dxx_h0(3:end-2);        
%     h0 = h0(3:end-2);                           % remove padding

    zEst0 = [zEst0(3), zEst0(2), zEst0, zEst0(end-2), zEst0(end-3)];
    dxx_Z = conv(zEst0,0.5*[1 -2 1],'same');
    dxx_Z = dxx_Z(3:end-2);
    zEst0 = zEst0(3:end-2);
    dp = -(params.f_r/params.f_l)*(params.d)./(zEst0.*zEst0);
    
    h1 = (params.f_r/params.f_l).*u.*zEst0 + params.f_r*params.b + u.*zEst0 - u*params.d;
    h2 = (params.f_r/params.f_l)*u+u;
    
    dhdz = h2./zEst0 - h1/(zEst0.*zEst0);
    
    I0zx = dp.* w0_left*params.pixelDim + p .* (-w0_left*params.pixelDim./zEst0) - (-right_adj./(zEst0 - params.d)).*dhdz;
    
    zEst0 = zEst0 + params.lambda*I0z.*I0zx + params.alpha^2*dxx_Z;

    % calcuate statistics
    rms_flow = sqrt(mean((I0z(200:600)).*conj(I0z(200:600))));
    %rms_Z = sqrt(mean((zEst0(10:end-10)-z0_left(10:end-10)).*conj(zEst0(10:end-10)-z0_left(10:end-10))))/mean(z0_left(10:end-10));


     subplot(2,3,1)
     plot(w0_left)
     hold all
     plot(w0_right)
     plot(right_adj./p)
     hold off
     %title({'Flow matching',strcat('Gradient Decent 1 Iterations:',num2str(j),' of:',num2str(params.iterations1))})
     title('Flow matching')
     legend('wL0', 'wR0','wL(x+h)') 
 
     subplot(2,3,2)
     plot(I0z)
     title({'I0z',strcat('RMS Flow Error:',num2str(rms_flow),' pixels')})
    %axis([0,400,-1e-3,1e-3])

    subplot(2,3,3)
    plot(h/params.pixelDim)
    title({'h',strcat('Alpha:',num2str(params.alpha),' Lamda:',num2str(params.lambda))}) 
    ylabel('Pixels')

    subplot(2,3,4)
    plot(params.lambda*I0z.*I0zx)
    hold all
    plot(params.alpha^2*dxx_Z)
    title('Components of EL')
    legend('I0xz','Smooth')
    hold off
    %axis([0,400,-1e-4,1e-4])

    subplot(2,3,5)
    plot(- params.lambda*I0z.*I0zx + params.alpha^2*dxx_Z)
    title('Gradient Descent Update')
    %axis([0,400,-1e-4,1e-4])

    subplot(2,3,6)
    plot(zEst0)
    %hold all
    %plot(z0_left)
    %hold off
    %title({'Z0 estimate vs. Z0 actual';strcat('RMS Z Error: ',num2str(rms_Z*100),'%')})


    %pause

    
    % shift the z estimate by the optical flow and resample.  This works because there is no delta Z
    % between times 0 and 1
    
    u = 1:1:length(zEst0);
    zEst1_u = u + w0_left;      % zEst1_u are the non-grid pixel coordinates of the depth estimate shifted from zEst0
    
    % now resample onto pixel grid
    zEst1 = resample2Grid(zEst0,zEst1_u);
    
    
end







