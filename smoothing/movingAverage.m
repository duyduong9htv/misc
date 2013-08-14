function out = movingAverage(t, x, t_window)
% function out = movingAverage(t, x, t_window)
% finds the moving average within a time window (in seconds)
% INPUTS: t (seconds)
%         x(data)
%         t_window (seconds)
% OUTPUT: the moving averaged value
% written by DD Tran. Last updated Aug 9, 2013. 

t_end = t(end); 

inds1 = find(t_end - t(:) > t_window); 

ind_end = inds1(end); 
out = x; 

for k = 1:ind_end
    inds = find((t(:) - t(k)- t_window).*(t(:) - t(k)) <= 0);
    out(k) = mean(x(inds));     
end



end
