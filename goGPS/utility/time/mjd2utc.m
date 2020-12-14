%**************************************************************************
%DESCRIPTION:                                                              *
%           mjd2utc converts Modified Julian Date(MJD) to UTC time [Year,..*
%           Month,Day,Hour Minute,Seconds].Thus it returns the ...         *
%           Gregorian calendar date (year, month,day, hour, minute,second) *
%           corresponding to the Julian day number JDAY.                   *
%NOTE:
%      Start of the JD (Julian day) count is from 0 at 12 noon 1 JAN -4712
%      (4713 BC), Julian proleptic calendar.Note that this day count conforms
%      with the astronomical convention starting the day at noon, in contrast
%      with the civil practice where the day starts with midnight.
%
%      Astronomers have used the Julian period to assign a unique number to
%      every day since 1 January 4713 BC.  This is the so-called Julian Day
%      (JD). JD 0 designates the 24 hours from noon UTC on 1 January 4713 BC
%      (Julian calendar) to noon UTC on 2 January 4713 BC.

%   Sources:  - http://tycho.usno.navy.mil/mjd.html
%             - The Calendar FAQ (http://www.faqs.org)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
% USAGE:                                                                   *
%       [Yr,Mn,D,Hr,Min,Sec]= mjd2utc(mjd);                                *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%                                                                          *
% INPUT:                                                                   *
%       mjd : Modified Julian Date                                         *

%***OUTPUTs:                                                               *
%1.         Yr : Year                                                      *
%2.         Mn : Month                                                     *
%3.         D  : Day                                                       *
%4.         Hr : Hour                                                      *
%5.        Min : Minute                                                    *
%6.        Sec : Seconds                                                   * 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================
% REFERENCE:                                                               *
%1.         Jean Meeus -"ASTRONOMICAL ALGORITHMS" second edition           *
%2.         Oliver Montenbruck-"Practical Ephemeris Calculations"          *
%3.         Peter Duffett-Smith-"Practical Astronomy with Your Calculator" *
%                            (4th Edition),Cambridge Univeristy Press,1988 * 
%4.         B.Hofmann-Wellenhof, H.Lichtenegger and J.Collins: GPS Theory..*
%                              and practice.2001. Fifth revised edition.   * 
%                                        Springer, Wien, New York.pp.37-38.         
%5.         Guochang Xu, Dr.-Ing:GPS Theory,Algorithms and Applications,...*
%                              June 2007. Second Edition.Springer,pp.37-38.*
%6.%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 12:50:30 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam
%==========================================================================+
%WRITTEN by: Samuel Osah,Msc Geomatic Engineering ,2016                    *    
%             Email: osahsamuel@yahoo.ca                                   *
%             Tel:+233 (0)246137410/+233 (0)509438484                      * 
%==========================================================================+
%***************************************************************************
function [Yr,Mn,D,Hr,Min,Secs] = mjd2utc(mjd)

%1. CALCULATE TIME [Hr Mun Secs] FROM MJD 

try
   Hr = floor((mjd-floor(mjd))*24);   % get hours
   Min = floor((((mjd-floor(mjd))*24)-Hr)*60);   % get minutes
   Secs = (((((mjd-floor(mjd))*24)-Hr)*60)-Min)*60;   % get seconds

catch %ALTERNATIVE
      %Call the "mjd2hms" fxn
      [Hr, Min, Secs] = mjd2hms(mjd);
     
end
 
%*****CHANGE SECs, MIN, HOUR WHOSE SEC == 60
Min(Secs==60) = Min(Secs==60)+1;
Secs(Secs==60) = 0;
Hr(Min==60) = Hr(Min==60)+1;
Min(Min==60)=0;
mjd(Hr==24)=mjd(Hr==24)+1;

%2.*******CALCULATION OF DAY,MONTH YEAR OF CALENDAR FROM JULIAN DATE

try
   %CALCULATE JULIAN DATE(JD)
   jd_all = mjd+2400000.5;

   % integer Julian date
   jd_all_int = floor(jd_all+0.5);

   aa = jd_all_int+32044;
   bb = floor((4*aa+3)/146097);
   cc = aa-floor((bb*146097)/4);
   dd = floor((4*cc+3)/1461);
   ee = cc-floor((1461*dd)/4);
   mm = floor((5*ee+2)/153);

   D  = ee-floor((153*mm+2)/5)+1; %Day of the Month
   Mn = mm+3-12*floor(mm/10);     %Month of the Year
   Yr = bb*100+dd-4800+floor(mm/10);%Year


catch %ALTERNATIVE
     %a,b, c, d, and e are auxiliary numbers. 
     a = floor(JD+0.5);

     try
        b = a+1537;
        c = floor((b-122.1)./365.25);
        d = floor(365.25.*c);
        e = floor((b-d)./30.6001);

        %**CALCULATE DAY OF THE MONTH
        try
           D = b-d-floor(30.6001.*e)+rem(JD+0.5,1);%DAY OF MONTH(DAY)
        catch 
             D = b - d - floor(30.6001.*e) + (JD + 0.5) - floor(JD + 0.5);
        end  %try

        %**MONTH OF THE YEAR
        Mn = e - 1 - 12.*floor(e./14);% MONTH OF YEAR

        %**YEAR
        Yr = c - 4715 - floor((7+Mn)./10);%YEAR

        %*****AN ALTERNATIVE CALCULATION
     catch 
          %Algorithm from "Practical Ephemeris Calculations"page34 by Oliver Montenbruck
          if a < 2299161
             c = a + 1524;
          else  
              b = floor( (a-1867216.25) / 36524.25 );
              c = a + b - floor(b/4) + 1525;
          end 
             d = floor( (c-122.1)/365.25 );
             e = floor(365.25.*d);
             f = floor( (c-e) ./ 30.6001 );
             D = c - e - floor(30.6001.*f) + rem((JD+0.5),a);
            Mn = f - 1 - 12.*floor(f./14);
            Yr = d - 4715 - floor( (7+Mn)./10 );
            
     end  %try
     
end

%SUBROUTINE TO CONVERT MJD TO [HOUR, MINUTE, SECONDS]
 function [hour, minute, second] = mjd2hms(mjd)
%**************************************************************************
%DESCRIPTION:
%            mjd2hms Convert days/mjd into hours, minutes, and seconds.
%
%           [HOUR, MINUTE, SECOND] = DAYS2HMS(DAYS) converts the number of 
%           days to hours, minutes, and seconds.
%The following holds (to within rounding precision):
%DAYS = HOUR / 24 + MINUTE / (24 * 60) + SECOND / (24 * 60 * 60)
%     = (HOUR + (MINUTE + SECOND / 60) / 60) / 24 
%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 12:52:02 +0100
 %   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam
 
fmjd = mjd - floor(mjd);%GET FRACTION PART
second = 86400 * fmjd;
hour   = fix(second/3600);           % get number of hours
second = second - 3600*hour;         % remove the hours
minute = fix(second/60);             % get number of minutes
second = second - 60*minute;         % remove the minutes


%SUBROUTINE TO CONVERT JD TO DATE[YEAR,MONTH,DAY,HOUR, MINUTE, SECONDS]
function [year, month, day, hour, minute, second] = jd2date(jd)
    
%[YEAR, MONTH, DAY, HOUR, MINUTE, SECOND] = JD2DATE(JD) returns the
% Gregorian calendar date (year, month, day, hour, minute, and second)
% corresponding to the Julian day number JD.
    
ijd = floor(jd + 0.5); % integer part
fjd = jd - ijd + 0.5;    % fraction part
[hour, minute, second] = days2hms(fjd);
 
%The following algorithm is from the Calendar FAQ.

a = ijd + 32044;
b = floor((4 * a + 3) / 146097);
c = a - floor((b * 146097) / 4);

d = floor((4 * c + 3) / 1461);
e = c - floor((1461 * d) / 4);
m = floor((5 * e + 2) / 153);

day   = e - floor((153 * m + 2) / 5) + 1;
month = m + 3 - 12 * floor(m / 10);
year  = b * 100 + d - 4800 + floor(m / 10);
