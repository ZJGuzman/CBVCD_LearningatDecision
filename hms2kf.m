function [ Snum ] = hms2kf( time_string1,time_string2, Tsec )
%HMS2KF calculates how many segments of Tsec length are in the time interval [time_string1:time_string2]
%   Inputs:
%       time_string1:hh:mm:ss
%		time_string2:hh:mm:ss
%
%   Output:
%       Snum: number of segments
%       
    seconds1=3600*str2double(time_string1(1:2))+60*str2double(time_string1(4:5))+ str2double(time_string1(7:min(8,length(time_string1))));
	seconds2=3600*str2double(time_string2(1:2))+60*str2double(time_string2(4:5))+ str2double(time_string2(7:min(8,length(time_string2))));
    if seconds1<=180 && seconds2<=180
        Snum=floor((seconds2-seconds1)/Tsec);
    else
        Snum=0;
    end
    %Number of initial keyframe at time_string time, maximum error of +fps frames

end

