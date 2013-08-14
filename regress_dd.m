function [b1,b0, SE] = regress_dd(x,y)
% function [b1,b0, SE] = regress_dd(x,y)
% performs the linear regression to get y = b1*x + b0
% SE is the standard error of the slope estimate 
% x, y are two input vectors of same dimension 

x_bar = mean(x); 
y_bar = mean(y); 
b1 = (sum((x-x_bar).*(y-y_bar)))/sum((x-x_bar).^2); 
b0 = y_bar - b1*x_bar; 

y_hat = b1*x + b0; 
numer = (1/length(x))*sum((y - y_hat).^2); 
denom = sum((x - x_bar).^2); 
SE = sqrt(numer/denom); 

end 