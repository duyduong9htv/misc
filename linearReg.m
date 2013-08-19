function y = linearReg(x, x1, x2, y1, y2)
% function y = linearReg(x, x1, x2, y1, y2)
% calculates the linear regression (y - x1)/(x - x1) = (y2 - y1)/(x2 - x1)
    y = y1 + (y2 - y1)/(x2 - x1)*(x - x1); 
end
