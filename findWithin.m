function [xWithin, inds] =  findWithin(xStart, xEnd, xData)
%function [xWithin, inds]  = findWithin(xStart, xEnd, xData)
%returns all the values of xData that lies within the interval [xStart, xEnd]
%xData is an input vector
inds =find((xData(:)-xStart).*(xData(:)-xEnd)<=0); 
xWithin = xData(inds); 

end
