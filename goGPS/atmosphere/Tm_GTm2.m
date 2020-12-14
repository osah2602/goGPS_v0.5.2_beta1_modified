%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *
%            "Tm_GTm2" is a subroutine estimates Weighted-Mean Temperature *
%             based on Spherical Harmonics up to degree and order  9       *                                           
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%USAGE:                                                                    *
%      Tm = Tm_GTm2(Time,lat, lon, hgt)                                    *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%1.     Time = Time of observation                                         *
% Time format:[year month day hour minute seconds] OR [year month day]     *                  
%2.      lat = Latitude of Station in degrees OR [D M S]                   *
%3.      lon = Longitude of Station in degrees OR [D M S]                  *
%4.        h = ellipsoidal height of Station  in m                         * 

%OUTPUT:                                                                   *
%       Tm = Weighted Mean Temperature                                     *
%       Tm : The key variable for calculating the exact conversion factor  * 
%            to map zenith wet delaysonto precipitable water               *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================
%REFERENCE:                                                                *
% YAO YiBin, ZHANG Bao, YUE ShunQiang
% Global empirical model for mapping zenith wet delays onto precipitable water
% =========================================================================
%  Code Editted by:                                                        *
%                  OSAH SAMUEL, Msc Geomatic Engineering ,2016             *    
%                  Email: osahsamuel@yahoo.ca                              *
%                  Tel:+233 (0)246137410/+233 (0)509438484                 * 
%==========================================================================
%***************************************************************************
%***************************************************************************
function Tm = Tm_GTm2(Time, Lat, Lon, h)

%**********CHECK & REFORMAT INPUTs DATA
switch nargin
    
   case {4,2} %Various inputs format 
       
        if (any(nargin==[4,2])) 
            
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
           
            if nargin==4 %IF all inputs are provided
                
               hgt=h;%Assigning h to hgt
               
               %****CHECK LATITUDE INPUT 
               ncol_lat= size(Lat,2); %Get number of latitude entry  
   
              switch ncol_lat
                  case  3     %if input is in DMS,Convert to degrees
                      lat=dms2degrees(Lat);%Lat in degrees      
                  case  2   %if input is in DM,Convert to degrees    
                      lat=dm2degrees(Lat);%Lat in degrees
                      
                  otherwise
                      lat=Lat;
              end       
          
            %*****CHECK LONGITUDE INPUT 
            ncol_lon= size(Lon,2); %Get number of longitude entry 
   
           switch ncol_lon
               case  3     %if input is in DMS,Convert to degrees
                    lon=dms2degrees(Lon);%Lat in degrees      
               case  2   %if input is in DM,Convert to degrees    
                   lon=dm2degrees(Lon);%Lat in degrees
                   
               otherwise 
                      lon=Lon;
           end       
               
            elseif  nargin == 2 %If input are two
                
                  %CHECK COORD ENTRY FORMAT
                   nCOL_lat=size(Lat,2); %Get number of Columns entry
                   
                   switch nCOL_lat
                       case 7
                           lat=dms2degrees(Lat(:,1:3));%Latitude in DMS (columns 1-3) converted to degrees 
                           lon=dms2degrees(Lat(:,4:6));%Longitude in DMS (columns 4-6) converted to degrees   
                           hgt=Lat(:,end);%Assigning end (7th) column to heights 
                   
                       case 6   
                           lat=dms2degrees(Lat(:,1:3));%Latitude in DMS (columns 1-3) converted to degrees 
                           lon=dms2degrees(Lat(:,4:6));%Longitude in DMS (columns 4-6) converted to degrees   
                           hgt=zeros(size(Lat,1),1);%Assigning zeros to  heights
                             
                       case 5
                           lat=dm2degrees(Lat(:,1:2));%Latitude in DM (columns 1-2) converted to degrees 
                           lon=dm2degrees(Lat(:,3:4));%Longitude in DM (columns 3-4) converted to degrees   
                           hgt=Lat(:,end);%Assigning end (5th) column to heights
                             
                      case 4
                           lat=dm2degrees(Lat(:,1:2));%Latitude in DM (columns 1-2) converted to degrees 
                           lon=dm2degrees(Lat(:,3:4));%Longitude in DM (columns 3-4) converted to degrees   
                           hgt=zeros(size(Lat,1),1);%Assigning zeros to  heights  
                             
                      case 3
                           lat=Lat(:,1);%Latitude in degrees(columns 1) 
                           lon=Lat(:,2);%Longitude in degrees(columns 2)
                           hgt=Lat(:,end);%Assigning end (3rd) column to heights       
    
                       case  2
                           lat=Lat(:,1);%Latitude in degrees(columns 1) 
                           lon=Lat(:,2);%Longitude in degrees(columns 2)
                           hgt=zeros(size(Lat,1),1);%Assigning zeros to  heights 
            
                   end %switch nCOL_lat
                   
            end %if nargin==4
            
        end %if (any(nargin==[4,2])) 
        
end %switch nargin
%===============================END OF INPUT REFORMATTING
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

%********COMPUTE Tm USING Global Weighted Mean Temperature II(GTmII)MODEL
%Call the "get_GTm2_Tm.m" Function
[Tm] = get_GTm2_Tm(UTCtime,lat,lon,hgt);


%=======================VARIOUS SUB-ROUTINES TO COMPUTE Tm=================
%                        ---------------------------------
function [ Tm ] = get_GTm2_Tm(UTCtime,lat,lon,h)

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
      Tm = zeros(nrow_time,1);
   
      for i = 1:nrow_time
          %Call the "PI_pwv.m" Function
          Tm(i,1) = GTm2_Tm(UTCtime(i,:),lat(i),lon(i),h(i)); 
          
      end
      
  else  
       %*****INITIALIZE OUTPUTs 
       Tm = zeros(nrow_pos,nrow_time);
      
     for i = 1:nrow_time %LOOP OVER TIME
         
        for j = 1:nrow_pos %LOOP OVER POSITIONS
            
            %Call the "PI_pwv.m" Function  
            Tm(j,i) = GTm2_Tm(UTCtime(i,:),lat(j),lon(j),h(j));
            
        end 
        
     end 
     
  end %//if isequal(Identical,1)
  
else 
    %*****INITIALIZE OUTPUTs 
    Tm = zeros(nrow_pos,nrow_time);
    
    for i = 1:nrow_time %LOOP OVER TIME
        
        for j = 1:nrow_pos %LOOP OVER POSITIONS
              
            %Call the "PI_pwv.m" Function  
            Tm(j,i) = GTm2_Tm(UTCtime(i,:),lat(j),lon(j),h(j));
            
        end
        
    end
    
end %//if isequal(nrow_time,nrow_pos)



function [Tm] = GTm2_Tm(UTCtime,lat,lon,hgt)
%**********MODEL COEFFICIENTS
alpha2 = -0.0060;
atm_mean = [...
   278.1606841258,0.7178150444,0.6347292577,-20.6515652902,0.9217568983,...
     0.1104730401,0.0460198181,0.1999965527,-0.0388377202,-0.0434545323,...
    -3.1375419701,0.2558079858,0.0385776907,-0.0142717844,-0.0029932841,...
     0.8331898140,-0.0052323500,0.0753666000,0.0006510563,0.0037060096,...
     0.0002524233,1.6935254615,0.0436064501,0.0404884560,0.0008855617,...
    -0.0001132782,0.0002322566,0.0000006505,0.4630789426,0.1497767500,...
     0.0194410417,0.0015387650,-0.0001215838,-0.0001072817,0.0000061194,...
    -0.0000012592,-1.2664469277,0.0120317046,0.0054847939,0.0019691557,...
     0.0001221583,0.0000064825,0.0000023866,0.0000011281,-0.0000001534,...
    -0.2249162117,0.1157422629,0.0055214419,-0.0013378342,-0.0000201478,...
    -0.0000154242,0.0000015295,0.0000004154,-0.0000000417,-0.0000000119];

atm_amp = [...
    -1.6278159856,-7.4866076088,-0.0814355070,-3.1321617531,-0.0396907415,...
     0.2103815364,-0.6682335585,0.3061679249,0.1958189200,-0.0099172634,...
     1.6730359926,0.3753883740,0.1923523875,-0.0094417306,0.0003882706,...
     0.8460217355,0.4287598411,0.0999903136,-0.0102842700,0.0003784148,...
     0.0001508969,1.4571398265,0.1053248869,0.0487426812,0.0000717202,...
    -0.0003073714,-0.0000277809,-0.0000103370,-2.0842918080,-0.0459665790,...
    -0.0043725881,-0.0015750142,0.0000552977,0.0000004218,-0.0000091303,...
    -0.0000002335,0.1284351248,0.1230643114,-0.0119477584,0.0000275078,...
     0.0002528090,-0.0000350604,-0.0000008112,-0.0000001873,0.0000000230,...
     1.1217545132,-0.0272857023,-0.0039719071,-0.0001943104,0.0000455921,...
    -0.0000019361,-0.0000015467,-0.0000001882,0.0000000052,0.0000000054];

btm_mean = [...
     0.0000000000,0.0000000000,0.1907979422,0.0000000000,0.0517578164,...
     0.1191773590,0.0000000000,-0.6855046479,0.1447426651,0.0283218455,...
     0.0000000000,0.2209564018,0.1146746488,0.0011285756,0.0047977552,...
     0.0000000000,-0.2490039047,0.0296616295,-0.0007635774,-0.0031564269,...
     0.0002694002,0.0000000000,0.3355813929, 0.0213109457,-0.0065350703,...
    -0.0007793893,-0.0002023246,0.0000033951,0.0000000000,0.0207759629,...
    -0.0126697133,-0.0023273625,-0.0003685592,-0.0000267648,-0.0000013347,...
    -0.0000019366,0.0000000000,0.0727089964,0.0158125363,0.0005017080,...
    -0.0001337302,0.0000324743,0.0000063718,0.0000000882,-0.0000001896,...
     0.0000000000,0.0487462426,-0.0141586861,0.0008277932,-0.0001081133,...
     0.0000105427,0.0000008221,0.0000001935,0.0000000090,0.0000000024];
 
btm_amp = [...
     0.0000000000,0.0000000000,-0.6588223039,0.0000000000,-0.8367533213,...
    -0.0481135800,0.0000000000,-0.4353679203,0.0632088817,-0.0148026445,...
     0.0000000000,-0.0111156025,0.0257839475,-0.0283506523,-0.0047639005,...
     0.0000000000,0.1422551799,0.0419245011,-0.0154521873,-0.0013611148,...
    -0.0000104406,0.0000000000,0.0539713509,0.0217023885,-0.0049635583,...
    -0.0005499290,0.0000434872,0.0000129990,0.0000000000,0.1149473270,...
    -0.0058603648,-0.0023545749,-0.0004258941,0.0000554082,0.0000034280,...
     0.0000020058,0.0000000000,0.0361976479,-0.0061496007,0.0007999286,...
    -0.0002155948,0.0000001467,-0.0000013105,0.0000007230,-0.0000000063,...
     0.0000000000,0.0889756187,0.0068989439,-0.0006929923,-0.0000073735,...
     0.0000085972,-0.0000016778,0.0000001579,0.0000000038,-0.0000000073];

a_geoid = [ ... 
-5.6195d-001,-6.0794d-002,-2.0125d-001,-6.4180d-002,-3.6997d-002, ...
+1.0098d+001,+1.6436d+001,+1.4065d+001,+1.9881d+000,+6.4414d-001, ...
-4.7482d+000,-3.2290d+000,+5.0652d-001,+3.8279d-001,-2.6646d-002, ...
+1.7224d+000,-2.7970d-001,+6.8177d-001,-9.6658d-002,-1.5113d-002, ...
+2.9206d-003,-3.4621d+000,-3.8198d-001,+3.2306d-002,+6.9915d-003, ...
-2.3068d-003,-1.3548d-003,+4.7324d-006,+2.3527d+000,+1.2985d+000, ...
+2.1232d-001,+2.2571d-002,-3.7855d-003,+2.9449d-005,-1.6265d-004, ...
+1.1711d-007,+1.6732d+000,+1.9858d-001,+2.3975d-002,-9.0013d-004, ...
-2.2475d-003,-3.3095d-005,-1.2040d-005,+2.2010d-006,-1.0083d-006, ...
+8.6297d-001,+5.8231d-001,+2.0545d-002,-7.8110d-003,-1.4085d-004, ...
-8.8459d-006,+5.7256d-006,-1.5068d-006,+4.0095d-007,-2.4185d-008]; 

b_geoid = [ ... 
+0.0000d+000,+0.0000d+000,-6.5993d-002,+0.0000d+000,+6.5364d-002, ...
-5.8320d+000,+0.0000d+000,+1.6961d+000,-1.3557d+000,+1.2694d+000, ...
+0.0000d+000,-2.9310d+000,+9.4805d-001,-7.6243d-002,+4.1076d-002, ...
+0.0000d+000,-5.1808d-001,-3.4583d-001,-4.3632d-002,+2.2101d-003, ...
-1.0663d-002,+0.0000d+000,+1.0927d-001,-2.9463d-001,+1.4371d-003, ...
-1.1452d-002,-2.8156d-003,-3.5330d-004,+0.0000d+000,+4.4049d-001, ...
+5.5653d-002,-2.0396d-002,-1.7312d-003,+3.5805d-005,+7.2682d-005, ...
+2.2535d-006,+0.0000d+000,+1.9502d-002,+2.7919d-002,-8.1812d-003, ...
+4.4540d-004,+8.8663d-005,+5.5596d-005,+2.4826d-006,+1.0279d-006, ...
+0.0000d+000,+6.0529d-002,-3.5824d-002,-5.1367d-003,+3.0119d-005, ...
-2.9911d-005,+1.9844d-005,-1.2349d-006,-7.6756d-009,+5.0100d-008];

%****CONVERT LATITUDE & LONGITUDE COORDS IN DEGREES TO RADIAN
lat=deg2rad(lat); %Latitude coord in radian
lon=deg2rad(lon); %Longitude coord in radian

%******COMPUTE DAY OF YEAR(doy) 
%***CALL THE DayofYear.m FUNCTION
 doy=DayofYear(UTCtime);

%REFERENCE DAY IS 28 JANUARY.THIS IS TAKEN FROM NIELL (1996) TO BE CONSISTENT
%***NEW DAY OF YEAR
doy=doy+1-28;

%***PARAMETERS t
t = sin(lat);

%****DEGREE n and ORDER m FOR LEGENDRE POLYNOMIAL
n = 9;
m = 9;

%***DETERMINE n! (faktorielle)  moved by 1
dfac(1) = 1;

for i = 1:(2*n + 1)
    dfac(i+1) = dfac(i)*i;
end 

%***DETERMINE LEGENDRE FUNCTIONS 
%REFERENCE:Heiskanen and Moritz,Physical Geodesy, 1967, eq. 1-62
for i = 0:n
    for j = 0:min(i,m)
       ir = floor((i - j)/2);
      sum = 0;
         for k = 0:ir
             sum = sum + (-1)^k*dfac(2*i - 2*k + 1)/dfac(k + 1)/......
                dfac(i - k + 1)/dfac(i - j - 2*k + 1)*t^(i - j - 2*k);
         end
         
       %LEGENDRE FUNCTIONs moved by 1
       P(i + 1,j + 1) = 1.d0/2^i*sqrt((1 - t^2)^(j))*sum;
    end 
end 

%*****SPHERICAL HARMONICS
i = 0;
for n = 0:9
   for m = 0:n
       i = i + 1;
       aP(i) = P(n+1,m+1)*cos(m*lon);
       bP(i) = P(n+1,m+1)*sin(m*lon);
   end 
end 
    
%****GEOIDAL HEIGHT
undu = 0.d0;
for i = 1:55
    undu = undu + (a_geoid(i)*aP(i) + b_geoid(i)*bP(i));
end 

%***ORTHOMETRIC HEIGHT
hort = hgt - undu;

%***SURFACE TEMPERATURE ON THE GEOID
atm = 0.d0;
ata = 0.d0;
for i = 1:55
    atm = atm + (atm_mean(i)*aP(i) + btm_mean(i)*bP(i));
    ata = ata + (atm_amp(i) *aP(i) + btm_amp(i) *bP(i));
end 

Tm0 =  atm + ata*cos(doy/365.25d0*2*pi);

%****HEIGHT CORRECTION
Tm = Tm0 +alpha2*hort;


%*********SUBROUTINE TO COMPUTE DAY OF YEAR(DoY)
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