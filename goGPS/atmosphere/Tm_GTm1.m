%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *
%            "Tm_GTm1" is a subroutine estimates Weighted-Mean Temperature*
%            (Tm) based on Spherical Harmonics up to degree and order  9   *                               *            
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%USAGE:                                                                    *
%      Tm = Tm_GTm1(Time,lat, lon, hgt)                                    *
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
% YAO YiBin, ZHU Shuang, YUE ShunQiang                                     *
% A globally-applicable, season-specific model for estimating the weighted *
% mean temperature of the atmosphere,to be submitted to Journal of Geodesy,*
% 2012.4                                                                   *
% =========================================================================
%  Code Editted by:                                                        *
%                  OSAH SAMUEL, Msc Geomatic Engineering ,2016             *    
%                  Email: osahsamuel@yahoo.ca                              *
%                  Tel:+233 (0)246137410/+233 (0)509438484                 * 
%==========================================================================
%***************************************************************************
%***************************************************************************
function Tm = Tm_GTm1(Time,Lat,Lon,h)

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

%********COMPUTE Tm USING Global Weighted Mean Temperature I(GTmI)MODEL
%Call the "get_GTm1_Tm.m" Function
[Tm] = get_GTm1_Tm(UTCtime,lat,lon,hgt);


%=======================VARIOUS SUB-ROUTINES TO COMPUTE Tm=================
%                        ---------------------------------
function [ Tm ] = get_GTm1_Tm(UTCtime,lat,lon,h)

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
          Tm(i,1) = GTm1_Tm(UTCtime(i,:),lat(i),lon(i),h(i)); 
          
      end
      
  else  
       %*****INITIALIZE OUTPUTs 
       Tm = zeros(nrow_pos,nrow_time);
      
     for i = 1:nrow_time %LOOP OVER TIME
         
        for j = 1:nrow_pos %LOOP OVER POSITIONS
            
            %Call the "PI_pwv.m" Function  
            Tm(j,i) = GTm1_Tm(UTCtime(i,:),lat(j),lon(j),h(j));
            
        end 
        
     end 
     
  end %//if isequal(Identical,1)
  
else 
    %*****INITIALIZE OUTPUTs 
    Tm = zeros(nrow_pos,nrow_time);
    
    for i = 1:nrow_time %LOOP OVER TIME
        
        for j = 1:nrow_pos %LOOP OVER POSITIONS
              
            %Call the "PI_pwv.m" Function  
            Tm(j,i) = GTm1_Tm(UTCtime(i,:),lat(j),lon(j),h(j));
            
        end
        
    end
    
end %//if isequal(nrow_time,nrow_pos)



function [Tm] = GTm1_Tm(UTCtime,lat,lon,hgt)
%*********MODEL COEFFICIENS
alpha2=-0.0041;

a_geoid = [ ...
    -5.6195e-001,-6.0794e-002,-2.0125e-001,-6.4180e-002,-3.6997e-002, ...
    +1.0098e+001,+1.6436e+001,+1.4065e+001,+1.9881e+000,+6.4414e-001, ...
    -4.7482e+000,-3.2290e+000,+5.0652e-001,+3.8279e-001,-2.6646e-002, ...
    +1.7224e+000,-2.7970e-001,+6.8177e-001,-9.6658e-002,-1.5113e-002, ...
    +2.9206e-003,-3.4621e+000,-3.8198e-001,+3.2306e-002,+6.9915e-003, ...
    -2.3068e-003,-1.3548e-003,+4.7324e-006,+2.3527e+000,+1.2985e+000, ...
    +2.1232e-001,+2.2571e-002,-3.7855e-003,+2.9449e-005,-1.6265e-004, ...
    +1.1711e-007,+1.6732e+000,+1.9858e-001,+2.3975e-002,-9.0013e-004, ...
    -2.2475e-003,-3.3095e-005,-1.2040e-005,+2.2010e-006,-1.0083e-006, ...
    +8.6297e-001,+5.8231e-001,+2.0545e-002,-7.8110e-003,-1.4085e-004, ...
    -8.8459e-006,+5.7256e-006,-1.5068e-006,+4.0095e-007,-2.4185e-008];

b_geoid = [ ...
    +0.0000e+000,+0.0000e+000,-6.5993e-002,+0.0000e+000,+6.5364e-002, ...
    -5.8320e+000,+0.0000e+000,+1.6961e+000,-1.3557e+000,+1.2694e+000, ...
    +0.0000e+000,-2.9310e+000,+9.4805e-001,-7.6243e-002,+4.1076e-002, ...
    +0.0000e+000,-5.1808e-001,-3.4583e-001,-4.3632e-002,+2.2101e-003, ...
    -1.0663e-002,+0.0000e+000,+1.0927e-001,-2.9463e-001,+1.4371e-003, ...
    -1.1452e-002,-2.8156e-003,-3.5330e-004,+0.0000e+000,+4.4049e-001, ...
    +5.5653e-002,-2.0396e-002,-1.7312e-003,+3.5805e-005,+7.2682e-005, ...
    +2.2535e-006,+0.0000e+000,+1.9502e-002,+2.7919e-002,-8.1812e-003, ...
    +4.4540e-004,+8.8663e-005,+5.5596e-005,+2.4826e-006,+1.0279e-006, ...
    +0.0000e+000,+6.0529e-002,-3.5824e-002,-5.1367e-003,+3.0119e-005, ...
    -2.9911e-005,+1.9844e-005,-1.2349e-006,-7.6756e-009,+5.0100e-008];

atm_amp= [ ...
    -2.7728e+000,-7.5283e+000,+2.1999e-002,-2.5937e+000,-1.6267e-001, ...
    -7.4794e-002,-1.3088e-001,+9.8054e-001,+1.3614e-001,+6.2514e-002, ...
    +4.1121e-001,+2.6168e-001,+3.3232e-001,-8.2088e-003,+2.2860e-002, ...
    -1.7374e+000,+4.3731e-001,-4.3279e-002,-4.0835e-002,+2.1399e-004, ...
    +1.6214e-003,+4.6134e+000,+5.6725e-001,+1.8581e-002,+2.7329e-002, ...
    -4.7682e-003,-1.4869e-005,+5.8632e-005,+6.2653e-003,-3.1505e-001, ...
    +4.2535e-003,+1.8852e-003,-5.5224e-004,+1.7282e-004,-8.3607e-005, ...
    -5.7291e-006,-2.1930e+000,-1.2753e-001,-7.2903e-003,-2.6333e-003, ...
    +8.4671e-004,-1.0379e-005,+4.5349e-007,+1.4445e-006,+1.3709e-008, ...
    +6.2457e-001,+8.6543e-002,+3.2264e-003,-9.6226e-005,-5.8426e-005, ...
    +5.1657e-007,-7.8229e-007,+4.3494e-007,+4.8682e-008,+8.8366e-009]; 

atm_mean= [ ...
    +2.8016e+002,-5.6872e+000,+3.8909e-002,-1.7315e+001,+3.2029e+000, ...
    +2.0950e-001,-4.2781e-001,-1.9802e+000,+3.2656e-001,-6.5104e-002, ...
    -9.4437e+000,+2.7176e+000,-7.2632e-001,-1.0385e-003,-3.6031e-002, ...
    +1.2180e+001,-3.8980e-001,+3.3896e-001,+1.0117e-001,-2.0203e-003, ...
    +4.0536e-003,-3.2068e+000,-1.2381e+000,+1.3275e-001,-7.7234e-002, ...
    +1.9781e-003,-1.8740e-005,+2.2050e-005,-1.8464e+000,+2.6792e-001, ...
    -9.0697e-002,+2.4401e-002,+8.1835e-004,-6.1739e-004,+7.5870e-005, ...
    -2.9726e-005,-1.0808e+000,+1.5702e-001,+3.0995e-002,+9.8326e-003, ...
    -1.9646e-003,+3.9563e-004,-1.7089e-005,-6.0835e-007,-2.4168e-007, ...
    -9.9527e-001,+1.6400e-001,+1.7958e-002,-6.3732e-003,+5.7614e-004, ...
    -5.2672e-005,+1.8227e-006,+1.1834e-006,-2.5700e-008,-4.3692e-008];

btm_amp= [ ...
    +0.0000e+000,+0.0000e+000,-1.4327e+000,+0.0000e+000,-1.5105e+000, ...
    -8.9129e-001,+0.0000e+000,+4.9271e-001,-4.0104e-002,-3.4000e-002, ...
    +0.0000e+000,-3.2911e-001,+2.6720e-001,-3.7570e-002,-1.5440e-002, ...
    +0.0000e+000,-3.4506e-001,+2.4431e-001,-2.8199e-002,+2.1516e-003, ...
    +1.5038e-003,+0.0000e+000,+6.0108e-001,-8.8919e-002,-1.1045e-002, ...
    +2.2734e-003,-3.5372e-004,+2.0760e-004,+0.0000e+000,+3.1097e-001, ...
    -8.3211e-002,-2.6425e-003,-5.8182e-005,-5.7868e-005,-1.3714e-006, ...
    +2.0388e-005,+0.0000e+000,-3.8079e-001,+3.3248e-002,+4.1894e-003, ...
    -4.9002e-004,+7.4654e-005,-1.6207e-005,+2.3754e-006,+6.0765e-008, ...
    +0.0000e+000,+1.8702e-001,+1.0714e-003,-8.9522e-004,+2.6448e-005, ...
    -4.5731e-005,+1.0969e-006,+9.6145e-008,-2.5005e-007,-7.5742e-009];

btm_mean= [ ...
    +0.0000e+000,+0.0000e+000,-1.0560e+000,+0.0000e+000,+4.0234e+000, ...
    +1.1051e+000,+0.0000e+000,-7.2625e+000,-2.1040e-001,-3.4823e-001, ...
    +0.0000e+000,+4.3705e+000,+7.7373e-001,+1.4833e-002,+5.2877e-003, ...
    +0.0000e+000,+3.9907e-001,-5.1635e-001,+3.6154e-002,+4.8566e-003, ...
    +2.0707e-003,+0.0000e+000,-1.9797e+000,+2.2000e-001,-9.3264e-003, ...
    -7.0120e-003,+4.9429e-004,-2.5810e-004,+0.0000e+000,+8.5808e-001, ...
    -7.3579e-003,-2.0269e-002,+4.3399e-003,-2.3579e-004,-1.7385e-005, ...
    -8.5941e-007,+0.0000e+000,+5.6363e-001,-3.2139e-003,+1.2155e-002, ...
    -1.4205e-003,-2.4200e-004,+4.7961e-005,-4.0606e-006,+4.5887e-007, ...
    +0.0000e+000,-4.7424e-001,-1.9732e-002,-1.9259e-003,+1.7334e-004, ...
    +8.1834e-005,-1.8045e-006,+5.7849e-007,+2.0996e-007,-3.1102e-008];

%difference between  modified julian date and Serial Date Number in Matlab
Tconst=678942;

%****CONVERT LATITUDE & LONGITUDE COORDS IN DEGREES TO RADIAN
lat=deg2rad(lat); %Latitude coord in radian
lon=deg2rad(lon); %Longitude coord in radian

%******COMPUTE DAY OF YEAR(doy) 
try
   %***CALL THE DayofYear.m FUNCTION
   doy=DayofYear(UTCtime);
 
catch
     %****COMPUTE MODIFIED JULIAN DATE
     mjd=mjuliandate(UTCtime);
     doy=mjd-mjuliandate(year(mjd+Tconst),1,1)+1;
end

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
undu = 0.0;
for i = 1:55
    undu = undu + (a_geoid(i)*aP(i) + b_geoid(i)*bP(i));
end

%***ORTHOMETRIC HEIGHT
hort = hgt - undu;

%***SURFACE TEMPERATURE ON THE GEOID
atm = 0.0;
ata = 0.0;
for i = 1:55
    atm = atm + (atm_mean(i)*aP(i) + btm_mean(i)*bP(i));
    ata = ata + (atm_amp(i) *aP(i) + btm_amp(i) *bP(i));
end
Tm0 =  atm + ata*cos(doy/365.25*2*pi);

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