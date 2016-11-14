function [ rmsErr ] = computErrorGC( leftFlow, depthMap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


left = 100;
right = 540;
top = 50;
bottom = 430;

fl = 1293.74844*.0048; 

depthMapXY = -depthMap;

imshow((-depthMapXY-500)/2500);


hold on
plot(left:right,top:top,'r')
plot(right:right,top:bottom,'r')
plot(left:right,bottom:bottom,'r')
plot(left:left,top:bottom,'r')
hold off;

stepSize = 1;

depth = [];
flow = [];
Vel = [];
if stepSize == 1
    index = 1;
    for i = left:right
        for j = top:bottom
            depth(index) = depthMapXY(j,i);
            flow(index) = leftFlow(j,i);
            Vel(index) = flow(index)*.006.*depth(index)/fl ;
            index = index + 1;
        end
    end
end
    

vMean = 20;
rmsErr = sqrt(mean((Vel-vMean).*conj(Vel-vMean)));
disp(strcat('RMS Error:',num2str((rmsErr/20)*100),'%'))

vMean = mean(Vel);
rmsErr = sqrt(mean((Vel-vMean).*conj(Vel-vMean)));
disp(strcat('RMS Error:',num2str((rmsErr/20)*100),'%'))

end

