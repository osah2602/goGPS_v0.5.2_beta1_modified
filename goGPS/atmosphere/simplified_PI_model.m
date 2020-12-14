%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *
%            "simplified_PI_model" is a subroutine that Computes the       *
%            Water Conversion Factor(PI) needed to convert Zenith Wet      * 
%            Delays(ZWD) to Precipitable Water Vapor(PWV) using User       * 
%            location/position and Time of observation for the computation * 
%            day-of-year(DoY).Generally, the simplified PI model uses      * 
%            user latitude and day-of-year(DoY) to compute PI(pie)- a      *
%            dimensionless conversion factor needed for the retrieval of   *
%            precipitable water vapor (PWV) from Global Navigation         *
%            Satellite System (GNSS)signal.
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%USAGE:                                                                    *                                                            *
%      [PI] = simplified_PI_model(Time,Lat,h)                              *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%1.     Time  = Time of observation                                        *
% Time format :[year month day hour minute seconds] OR [year month day]    *                  
%2.       lat = Latitude of Station in degrees  OR [D M S]                 *
%3.         h = Height of Station                                          *                                                

%OUTPUT:                                                                   *
%       PI = Water Conversion factor                                       *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================
%REFERENCE:                                                                *
%          Manandhar, S., Lee, Y. H., Meng, Y. S., & Ong, J. T. (2017). A  * 
%          simplified model for the retrieval of precipitable water vapor  *
%          from GPS signal. IEEE Transactions on Geoscience and Remote ... *
%          Sensing, 55(11), 6245-6253.
% =========================================================================
%  Code Written by:                                                        *
%                  OSAH SAMUEL, Msc Geomatic Engineering ,2016             *    
%                  Email: osahsamuel@yahoo.ca                              *
%                  Tel:+233 (0)246137410/+233 (0)509438484                 * 
%==========================================================================
%***************************************************************************
%***************************************************************************

function [PI] = simplified_PI_model(Time,Lat,h)

%**********CHECK & REFORMAT INPUTs DATA
switch nargin
    
   case {3,2} %Various inputs format 
       
        if (any(nargin==[3,2])) 
            
           %*****CHECK TIME FORMAT
           ncol=size(Time,2);% finding Number of columns 
           
           switch ncol              
               case  {6,5,4,3}
                  
                   if (any(ncol==[6,5,4,3]))
                    
                     Date=Time(:,1:3);%Assigning column 1 to 3 of time to date                     
                     %***Date
                     Yr=Date(:,1);%Assigning 1st column to years
                     Mn=Date(:,2);%Assigning 2ND column to months
                     Day=Date(:,3);%Assigning 3RD column to days
               
                     %***********CHECK FOR 2 DIGIT YEARS & CONVERT THE YEAR TO 4 DIGITS
                     %                              +            ================
                     %NOTE:                   
                     %     Two digit year represents a year in the range 1980-2079.

                     if (Yr >= 80 & Yr <= 99)
                       Yr = 1900 + Yr;
                     end  
         
                     if (Yr >= 0 & Yr <= 79)
                       Yr = 2000 + Yr;
                     end   
                     
                     if ncol==6                                                 
                       time = Time(:,4:6);%Assigning column 4 to 6 of tyme to time                                            
                       %***Time
                       Hr   = time(:,1);%Assigning 1st column to hours    
                       Min  = time(:,2);%Assigning 2nd column to min  
                       Secs = time(:,3);%Assigning 3rd column to sec
                  
                     elseif   ncol==5                         
                           time = Time(:,4:5);%Assigning column 4 to 5 of to tyme                          
                           %Time
                           Hr   = time(:,1);%Assigning 1st column to hours    
                           Min  = time(:,2);%Assigning 2nd column to min               
                           Secs = zeros(size(Hr,1),1);%Assigning zeros to Seconds 
                       
                     elseif   ncol==4 
                           time = Time(:,4);%Assigning column 4 of time to tyme
                           Hr   = time(:,1);%Assigning 1st column to hours 
                           Min  = zeros(size(Hr,1),1);%Assigning zeros to Minute
                           Secs = zeros(size(Hr,1),1);%Assigning zeros to Seconds
                       
                     elseif   ncol==3
                         
                             Hr   = zeros(size(Yr,1),1);%Assigning zeros to Hour
                             Min  = zeros(size(Yr,1),1);%Assigning zeros to Minute
                             Secs = zeros(size(Yr,1),1);%Assigning zeros to Seconds
                       
                     end  %if ncol==6 
                     
                     %*****CHANGE SECs, MIN, HOUR WHOSE SEC == 60
                     Min(Secs==60) = Min(Secs==60)+1;
                     Secs(Secs==60) = 0;
                     Hr(Min==60) = Hr(Min==60)+1;
                     Min(Min==60)=0;
                     
                    %CONCATENATE TIME ELEMENTS
                    UTCtime = [Yr Mn Day Hr Min Secs];%UTC time
                    
                   end   %if (any(ncol==[6,5,4,3]))
                   
           end %switch ncol
           
           if nargin == 3
             
              lon = [];%#ok<*NASGU> %Assigning empty([]) matrix to longitude coord
            
              %****CHECK LATITUDE INPUT 
              ncol_lat= size(Lat,2); %Get number of latitude entry  
   
              switch ncol_lat
                  case  3     %if input is in DMS,Convert to degrees
                      lat = dms2degrees(Lat);%Lat in degrees 
                      
                  case  2   %if input is in DM,Convert to degrees    
                      lat = dm2degrees(Lat);%Lat in degrees
                      
                  otherwise  
                           lat = Lat;
              end 
            
           elseif nargin == 2
             
                  %CHECK COORD ENTRY FORMAT
                  nCOL_lat=size(Lat,2); %Get number of Columns entry
                   
                      switch nCOL_lat
                          case 7
                               lat = dms2degrees(Lat(:,1:3));%Latitude in DMS (columns 1-3) converted to degrees 
                               lon = dms2degrees(Lat(:,4:6));%Longitude in DMS (columns 4-6) converted to degrees   
                                 h = Lat(:,end);%Assigning end (7th) column to heights 
                   
                          case 6   
                               lat = dms2degrees(Lat(:,1:3));%Latitude in DMS (columns 1-3) converted to degrees 
                               lon = dms2degrees(Lat(:,4:6));%Longitude in DMS (columns 4-6) converted to degrees   
                                 h = zeros(size(Lat,1),1);%Assigning zeros to  heights
                             
                          case 5
                               lat = dm2degrees(Lat(:,1:2));%Latitude in DM (columns 1-2) converted to degrees 
                               lon = dm2degrees(Lat(:,3:4));%Longitude in DM (columns 3-4) converted to degrees   
                                 h = Lat(:,end);%Assigning end (5th) column to heights
                             
                          case 4
                               lat = dm2degrees(Lat(:,1:2));%Latitude in DM (columns 1-2) converted to degrees 
                               lon = dm2degrees(Lat(:,3:4));%Longitude in DM (columns 3-4) converted to degrees   
                                 h = zeros(size(Lat,1),1);%Assigning zeros to  heights  
                             
                          case 3
                               lat = Lat(:,1);%Latitude in degrees(columns 1) 
                               lon = Lat(:,2);%Longitude in degrees(columns 2)
                                 h = Lat(:,end);%Assigning end (3rd) column to heights       
    
                          case 2
                               lat = Lat(:,1);%Latitude in degrees(columns 1) 
                               lon = Lat(:,2);%Longitude in degrees(columns 2)
                                 h = zeros(size(Lat,1),1);%Assigning zeros to  heights 
            
                      end  %switch nCOL_lat
                      
                      
           end   %//if nargin == 3
             
        end %//if (any(nargin==[3,2]))
        
    otherwise  
             beep
             fprintf('\n\nInsuficient Inputs for computation of PI. Check inputs & Try again.\n\n')     
             PI = [];
             return
        
end %switch nargin
%===============================END OF INPUT REFORMATTING
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

%********COMPUTE PI USING THE SIMPLIFIED PI MODEL
%Call the "getPI.m" Function
[PI] = getPI(UTCtime,lat,h);


%============================VARIOUS SUB-ROUTINES==========================
%                            --------------------
function [PI] = getPI(UTCtime,lat,h)

%GET SIZE OF USER INPUT TIME & POSITION
nrow_time = size(UTCtime,1);
nrow_pos = size(lat,1);

%1.******UTCdate
Yr  = UTCtime(:,1);%get Year
Mn  = UTCtime(:,2);%get Month
Day = UTCtime(:,3);%get Day

if isequal(nrow_time,nrow_pos)
      
   if nrow_time > 1 %IF SIZE OF TIME IS GREATER 1
       
      %TESTING FOR IDENTICAL NUMBERS IN(YR,Mn,DAY) MATRIX ROWs
      %--------------------------------------------------------------------
      %NOTE:
      %*Y= diff(X) calculates differences between adjacent elements of X along 
      %the first array dimension whose size does not equal 1: 
      %The elements of Y are the differences between adjacent elements of X.
      %Y = [X(2)-X(1) X(3)-X(2) ... X(m)-X(m-1)].
      %The not (~) converts the differences to logicals(0 difference->true). 
      %Then all sees if all the differences in a column are 0
      %*all(bsxfun(@eq,Yr,Yr(1,:))) ====> (This actually has one advantage, 
      %in that it would work even if A has only one row!)
      %--------------------------------------------------------------------
      Identical = any([all([all(~diff(Yr)),all(~diff(Mn)),all(~diff(Day))]),...
                       all([all(bsxfun(@eq,Yr,Yr(1,:))),all(bsxfun(@eq,Mn,Mn(1,:))),...
                       all(bsxfun(@eq,Day,Day(1,:)))])]);               
   else 
       Identical = all([all(bsxfun(@eq,Yr,Yr(1,:))),all(bsxfun(@eq,Mn,Mn(1,:))),...
                        all(bsxfun(@eq,Day,Day(1,:)))]);
   end 
       
  if isequal(Identical,1) %if ROWS ARE IDENTICAL
       
      %*****INITIALIZE OUTPUTs 
      PI = zeros(nrow_time,1);
   
      for i = 1:nrow_time
          %Call the "PI_pwv.m" Function
          PI(i,1) = PI_pwv(UTCtime(i,:),lat(i),h(i)); 
          
      end
      
  else  
       %*****INITIALIZE OUTPUTs 
       PI = zeros(nrow_pos,nrow_time);
      
     for i = 1:nrow_time %LOOP OVER TIME
         
        for j = 1:nrow_pos %LOOP OVER POSITIONS
            
            %Call the "PI_pwv.m" Function  
            PI(j,i) = PI_pwv(UTCtime(i,:),lat(j),h(j));
            
        end 
        
     end 
     
  end %//if isequal(Identical,1)
  
else 
    %*****INITIALIZE OUTPUTs 
    PI = zeros(nrow_pos,nrow_time);
    
    for i = 1:nrow_time %LOOP OVER TIME
        
        for j = 1:nrow_pos %LOOP OVER POSITIONS
              
            %Call the "PI_pwv.m" Function  
            PI(j,i) = PI_pwv(UTCtime(i,:),lat(j),h(j));
            
        end
        
    end
    
end %//if isequal(nrow_time,nrow_pos)


function [PI] = PI_pwv(UTCtime,lat,h)

%******GET DATE[YEAR MONTH DAY]
Yr  = UTCtime(:,1);%get Hour
Mn  = UTCtime(:,2);%get Month
Day = UTCtime(:,3);%get Day

%**********COMPUTE DAY OF YEAR(DoY)
%Call the "DayofYear function.m"
DoY=DayofYear(Yr,Mn,Day);

%*****************DETERMINE hfac
%NB:
%  hfac is introduced to account for the difference between Northern & 
%  Southern Hemispheres.It is selected depending on whether latitude(lat) 
%  values are +ve(from North) or -ve(from south). The term sign(lat)
%  denotes whether the station is from the Northern hemisphere or the
%  Southern hemisphere.sign(lat) is 1 for stations from the Northern
%  hemisphere and -1 for stations from the Southern hemisphere.
    if sign(lat) == 1 %IF Northern Hemisphere
        hfac = 1.48;
    elseif sign(lat) == -1 %IF Southern Hemisphere
           hfac = 1.25;       
    end
        
%DETERMINE f
%NB:
%   f is neglgible or 0 for stations the low-altitude(h < 1000m) region and
%   f = -2.38e-6*h for stations from the high-altitude(h > 1000m.
if h < 1000
   f = 0;
elseif h > 1000
    
       f = -2.38e-6*h;       
end

c = -1.*sign(lat).*1.7e-5.*abs(lat).^(hfac-0.0001);
d = 0.165-(1.7e-5).*abs(lat).^(1.65) + f;

%COMPUTE PI USING THE SIMPLIFIED PI MODEL
PI = c.*cos(((DoY-28)./365.25).*2.*pi) + d;



%A.*********SUBROUTINE TO COMPUTE DAY OF YEAR(DoY)
function DoY=DayofYear(Year,Month,Day)

switch nargin    
    case 1
          if size(Year,2)>=3
             Yr=Year(:,1);%Assigning 1st column to years
             Mn=Year(:,2);%Assigning 2ND column to months
             D=Year(:,3);%Assigning 3RD column to days             
          end          
    otherwise        
             Yr=Year;%Assigning Year to Yr
             Mn=Month;%Assigning Month to Mn
             D=Day;%Assigning Day to D
             
end %switch nargin
   
%***********CHECK FOR 2 DIGIT YEARS & CONVERT THE YEAR TO 4 DIGITS
%NOTE:                               +            ================
%     Two digit year represents a year in the range 1980-2079.
%Call the FourDigitYear.m function
[Yr] = FourDigitYear(Yr);


%******CALCULATION OF DAY OF YEAR(DoY)
try
   I = ones(length(Yr),1);
 DoY = datenum([Yr Mn D]) - datenum([Yr I I]) + 1;
catch    
     %Compute Julian Date
         jDay = juliandate(Yr,Mn,D,0,0,0);%Midnight Morning of this Day
     jDayJan1 = juliandate(Yr,1,1,0,0,0);%Midnight Morning of January 1st
    dayNumber = jDay - jDayJan1 + 1;
     DoY=dayNumber;%Day of the Year     
end %try
%******************************************END OF DayofYear.m 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%A.1*********SUBROUTINE TO CONVERT TWO DIGITs TO YEAR 4 DIGITS
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
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-