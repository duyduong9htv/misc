% time2secs.m converts the time written out as jd hh mm ss and converts it to seconds 

function seconds = time2secs(time)

jd = floor(time/1000000);
tm = mod(time,1000000);
hh = floor(tm/10000);
ms = mod(tm,10000);
mm = floor(ms/100);
ss = mod(ms,100);
seconds = jd*24*3600 + hh*3600 + mm*60 + ss;
