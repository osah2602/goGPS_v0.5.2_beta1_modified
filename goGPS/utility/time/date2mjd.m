%**************************************************************************
%DESCRIPTION:                                                              *
%           "date2mjd" converts UTC Date/Time[Year,Month,Day,Hour Minute,...
%           Seconds] to Julian Date(JD)and Modified Julian Date (MJD)      *                                 *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
% USAGE:                                                                   *
%       [MJD,JD]= date2mjd(date);                                          *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%                                                                          *
% INPUT:                                                                   *
%       date = date [year, month, day, hour, min, sec]                     *

%***OUTPUTs:                                                               *
%           JD  = julian day                                               *
%           MJD = modified julian day                                      * 
%==========================================================================+
%***************************************************************************
function [MJD,JD]=date2mjd(date)

%(2)****GET VARIOUS TIME INPUTs 

%(2.1)******DATE[YEAR MONTH DAY]
Yr  = date(:,1);%get Hour
Mn  = date(:,2);%get Month
Day = date(:,3);%get Day

%(2.2)******TIME[HOUR MINUTE SECONDs]
if size(date,2) == 6
   H    = date(:,4);%get Hour
   MIN  = date(:,5);%get Minute
   SECs = date(:,6);%get Seconds
   
   %CHANGE SECs,MIN,HOUR WHOSE SEC==60
   MIN(SECs==60) = MIN(SECs==60)+1;
   H(MIN==60) = H(MIN==60)+1;
     
elseif size(date,2) == 3
       H    = zeros(size(date,1),1);%create Hour of zeros
       MIN  = zeros(size(date,1),1);%get Minute of zeros
       SECs = zeros(size(date,1),1);%get Seconds of zeros    
end
    
%(5)********COMPUTE MODIFIED JULIAN DATE(MJD)
%Call the "utc2JulianDay_DoY.m" function
 [MJD ,JD]=utc2JulianDay(Yr,Mn,Day,H,MIN,SECs);

%SUBROUTINE TO COMPUTE JULIAN & MODIFIED DAY(JD,MJD)
function [MJD,JD]=utc2JulianDay(Year,Month,Day,Hour,Minute,Seconds)

%**************************************************************************
%***DESCRIPTION:
%               The function "JulianDay" Computes Julian &  Modified 
%               Julian Day as well as Day of Year (DoY) from Calender date 
 
%USAGE: [MJD,JD]=JulianDay(Year,Month,Day,Hour,Minute,Seconds)
  
%**INPUTs
%DATE
%1.  Year: Two digit calendar year representing the range 1980-2079...
%                                           (e.g.,99 = 1999 and 01 = 2001).
%2. Month: Month number of the year, must be in range 1-12 ....
%                                   (January = 1, Feb = 2, March = 3, etc.)
%3. Day  : day number of the month(calendar day), must be in range 1-31

%TIME
%4.  Hour    - hour (UTC), must be in range 0-24.
%5.  Minutes - minutes (UTC), must be in range 0-59.
%6.  Seconds - seconds (UTC), must be in range 0-59.
  
%OUTOUT :
  %1.     JD -------->Julian Date
  %2.    MJD -------->Modified Julian Date
 
%REFERENCE:
%1.          Jean Meeus -"ASTRONOMICAL ALGORITHMS" second edition
%2.          Oliver Montenbruck-"Practical Ephemeris Calculations"
%3.          Peter Duffett-Smith-"Practical Astronomy with Your Calculator"
%                            (4th Edition),Cambridge Univeristy Press, 1988
%**************************************************************************+
%==========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%**************************************************************************+
%**************************************************************************+    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%********************ASSIGNMENT 
%*******DATE
Yr=Year;%Assigning Year to Yr
Mn=Month;%Assigning Month to Mn
D=Day;%Assigning Day to D
          
%***********CHECK FOR 2 DIGIT YEARS & CONVERT THE YEAR TO 4 DIGITS
%NOTE:                               +            ================
%     Two digit year represents a year in the range 1980-2079.

%Call the FourDigitYear.m function
[Yr] = FourDigitYear(Yr);

%********TIME
H=Hour;%Assigning Hour to H    
MIN=Minute;%Assigning Minute to MIN 
SECs=Seconds;%Assigning Seconds to SECs 

%*********CONVERT TIME IN [Hour Min Sec] TO HOURs               
UT = H + (MIN / 60) + (SECs / 3600);% Time in Hours
             
%**************COMPUTE JULIAN DATE ......
          
%*****If The Month is 1 or 2 then we consider january and february .... 
%                              to be the 13th or 14th of the preceding year          
if (Mn<=2)  
  Y=Yr-1;  
  M=Mn+12;  
else     
    M=Mn;
    Y=Yr;
end  %(Mn<=2)
         
if any(UT(:) ~=0) 
    
  try
     JD = floor( (365.25 * Y) ) + floor( (30.6001 * (M+1)) ) + D + ( UT / 24 ) + 1720981.5;
  catch      
       JD = floor( (365.25 * Y) ) + floor( (30.6001 * (M+1)) )- 15 + 1720996.5 + D + ( UT / 24 );    
  end   %try
                   
%****CALCULATION OF MODIFIED JULIAN DATE(MJD)
MJD = JD-2400000.5;%Modified Julian date

else
    JD  = juliandate(Yr, Mn, D);
    MJD = mjuliandate(Yr, Mn, D);
end

%*********SUBROUTINE TO CONVERT TWO DIGITs TO YEAR 4 DIGITS
function [yearOUT] = FourDigitYear(yearIN)
%***************************************************************************
%***************************************************************************
%DESCRIPTION: 
%            "FourDigitYear"Checks for two digit year values & convert them
%             to four-digit year values
%NOTE:                   
%     Two digit year represents a year in the range 1980-2079.
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% USAGE:
%       [yearOUT] = four_digit_year(yearIN);

%INPUT:
%      yearIN = two-digit year

%OUTPUT:
%   yearOUT = four-digit year
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================

yearOUT = zeros(size(yearIN));

try
   yearOUT(yearIN >= 80 & yearIN <= 99) = yearIN(yearIN >= 80 & yearIN <= 99)  + 1900;
yearOUT(yearIN >= 0 & yearIN <= 79)  = yearIN(yearIN >= 0 & yearIN <= 79)  + 2000; 

catch
     yearOUT(yearIN > 79)  = yearIN(yearIN > 79)  + 1900;
     yearOUT(yearIN <= 79) = yearIN(yearIN <= 79) + 2000;
end
%***CHECK FOR UNCONVERTED YEARs
%This happens when a particular year value is in four digits already

%Get index for an already four digit year value.It is zero in the yearOUT
%matrix
Index=(yearOUT==0);
Index_v=yearIN(Index);%Get the year value of that index from yearIN
yearOUT(Index)=Index_v;%Assign the year value to corresponding row in yearOUT

%******************************************FourDigitYear.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-