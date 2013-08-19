function secs = dateTimeToSeconds(dateTime)
%function secs = dateTimeToSeconds(dateTime)
%converts the string dateTime in a whale
%property to seconds from 00:00 of the same day. 
%answer is in EDT. 

hh = dateTime(12:13); 
if str2num(hh) < 10
   hh_f = str2num(hh)+24; 
else
    hh_f = str2num(hh);
end

mm = dateTime(15:16);
ss = dateTime(18:end); 
secs = hh_f*60*60 + str2num(mm)*60 + str2num(ss); 

if secs < 3600*4
    secs = 24*3600 + secs - 3600*4;
else 
    secs = secs - 3600*4; 
end

end 