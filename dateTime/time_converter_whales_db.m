function [yyyy, month, dd, hh, mm, ss] = time_converter_whales_db(dateTime, TimeFlag)
% function [yyyy, mm, dd, hh, mm, ss] = time_converter_whales_db(dateTime)
% converts the attribute dateTime (string)
% obtained from a whale call ID
% into year, month, date and time in double type numbers. 
% TimeFlag by default is set to zero (EDT)
% if TimeFlag ~= 0 then the values are converted into GMT
% for example: 
% whales(10).dateTime
% 
% ans =
% 
% 2006-10-01 00:21:28.0
% [yyyy, mm, dd, hh, mm, ss] = time_converter_whales_db(dateTime, 1)
% yyyy = 2006
% month = 10
% dd = 1
% hh = 0
% mm = 21
% ss = 28.0

yyyy = str2double(dateTime(1:4));
month = str2double(dateTime(6:7));
dd = str2double(dateTime(9:10));
hh = str2double(dateTime(12:13));
mm = str2double(dateTime(15:16));
ss = str2double(dateTime(18:end));

if (nargin == 1) 
   TimeFlag = 0 ;  
end

if TimeFlag ==0
    if hh < 4
        hh = hh+24-4;
        if dd > 1
            dd = dd - 1;
        else 
            dd = 30;
            month = 9;
        end
    else
        hh = hh - 4; 
    end
   
end


            
        


