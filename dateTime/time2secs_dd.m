function sec = time2secs_dd(time_stamp)
%function sec = time2secs_dd(time_stamp)
%returns the time in seconds

time_stamp = mod(time_stamp, 1000000); 

hh = floor(time_stamp/10000); 
mm = floor((time_stamp - hh*10000)/100); 
ss = mod(time_stamp, 100); 

sec = hh*3600 + mm*60 + ss; 