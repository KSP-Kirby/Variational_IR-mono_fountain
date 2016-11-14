function [ w_bOut ] = resample2Grid( w_b, p_b )
%Takes the back flow (w_b) and the subpizel (nongrid) coordinate of each
% pixel and returns a grid sampling 
%w_b and p_b must be the same length
%p_b(1) must = 0
%extra pixel coordinates are filled with last value of w_b

    w_bOut = zeros(size(w_b));
    w_bOut(1) = w_b(1);
    for i=2:length(w_b)
        for j = 2:length(p_b)
            if p_b(j) > i
                i1 = p_b(j-1);
                i2 = p_b(j);
                y1 = w_b(j-1);
                y2 = w_b(j);
                break
            end
        end
        w_bOut(i) = (i-i1)/(i2-i1)*(y2-y1)+y1;     
    end
end

