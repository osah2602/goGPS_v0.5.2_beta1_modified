%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *
%            "Tm_TVGG_ERAI" is a subroutine that estimates Time-Varying    *
%             Global Gridded Ts-Tm(TCGG)[Weighted-Mean Temperature(Tm)] at *
%             a given User Position & time of observation.                 *                                     
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%USAGE:                                                                    *                                                           *
%      [Tm] = Tm_TVGG_ERAI(Time,Lat,Lon,Temp,gridV_TVGG)                   *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%1.     Time = Time of observation                                         *
% Time format= [year month day hour minute seconds]                        *                  
%2.      Lat = Latitude of Station in degrees OR [D M S]                   *
%3.      Lon = Longitude of Station in degrees OR [D M S]                  *
%4.     Temp = Surface Temperature in degree Celcius(°C)                   * 
%5.gridV_TVGG= TVGG-Tm model grid values extracted from Coeff_TVGG_ERAI.mat*
%              file                                                        *

%OUTPUT:                                                                   *
%       Tm = Weighted Mean Temperature                                     *
%       Tm : The key variable for calculating the exact conversion factor  * 
%            to map zenith wet delays onto precipitable water              *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================
%REFERENCE:                                                                *
%          Jiang, P., Ye, S., Lu, Y., Liu, Y., Chen, D., & Wu, Y. (2019).  *
%          Development of time-varying global gridded T s–T m model for    * 
%          precise GPS–PWV retrieval. Atmospheric Measurement Techniques,..
%          12(2), 1233-1249.                                               *
% =========================================================================
%  Code Editted by:                                                        * 
%                  OSAH SAMUEL, Msc Geomatic Engineering ,2016             *    
%                  Email: osahsamuel@yahoo.ca                              *
%                  Tel:+233 (0)246137410/+233 (0)509438484                 * 
%==========================================================================
%***************************************************************************
%***************************************************************************

function [Tm] = Tm_TVGG_ERAI(Time,Lat,Lon,Temp,gridV_TVGG)

%**********CHECK & REFORMAT INPUTs DATA
switch nargin
    
   case {5,4,3} %IF ALL INPUTS ARE PROVIDED 
       
       if (any(nargin==[5,4,3]))
           %*****CHECK TIME FORMAT
           ncol = size(Time,2);% finding Number of columns 
           
           switch ncol
               
               case  {6,5,4,3}
                  
                   if (any(ncol==[6,5,4,3]))
                    
                     Date = Time(:,1:3);%Assigning column 1 to 3 of time to date                     
                     %***Date
                     Yr  = Date(:,1);%Assigning 1st column to years
                     Mn  = Date(:,2);%Assigning 2ND column to months
                     Day = Date(:,3);%Assigning 3RD column to days
               
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
                   
           end  %switch ncol
           
           if nargin == 5 %IF ALL INPUTS ARE PROVIDED
                
              %CONVERT TEMPERATURE(T) IN Celcius(°C) to Kelvin(K) 
               T = Temp + 273.15; %Temperature in kelvin(K)

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
          
               %*****CHECK LONGITUDE INPUT 
               ncol_lon = size(Lon,2); %Get number of longitude entry 
   
               switch ncol_lon
                   case  3     %if input is in DMS,Convert to degrees
                       lon = dms2degrees(Lon);%Lat in degrees      
                   case  2   %if input is in DM,Convert to degrees    
                       lon = dm2degrees(Lon);%Lat in degrees
                   
                   otherwise   
                            lon = Lon;
               end 
                  
           elseif nargin == 4 %IF 4 INPUTS ARE PROVIDED
                   
                  %CHECK IF TEMP POSITION IS NOT TVGG GRID VALUES(gridV_TVGG)
                  %NOTE:
                  %     THIS HAPPENS WHEN USER DECIDE TO PROVIDE LAT & LON
                  %     COORD AS A SINGLE INPUT,I.E.[Lat Lon].IN THAT CASE,Temp
                  %     WILL TAKE THE POSITION OF LON & gridV_TVGG THAT OF Temp
                  try
                     if all([size(Temp,1) ~= 131072,size(Temp,1) ~= 8])
               
                     %CONVERT TEMPERATURE(T) IN Celcius(°C) to Kelvin(K) 
                     T          = Temp + 273.15; %Temperature in kelvin(K) 
                     gridV_TVGG = [];%Assigning empty matrix ([]) to gridV_TVGG
                     
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
          
                     %*****CHECK LONGITUDE INPUT 
                     ncol_lon = size(Lon,2); %Get number of longitude entry 
   
                     switch ncol_lon
                         case  3     %if input is in DMS,Convert to degrees
                             lon = dms2degrees(Lon);%Lat in degrees      
                         case  2   %if input is in DM,Convert to degrees    
                             lon = dm2degrees(Lon);%Lat in degrees
                   
                         otherwise    
                                  lon = Lon;
                     end  
                       
                  else %IF gridV_TVGG ASSUMES Temp POSITION. Temp ALSO ASSUMES Lon POSITION 
                      
                      %ASSIGNMENT
                      gridV_TVGG = Temp; %ASSUMING Temp POSITION TO BE gridV_TVGG
                      Temp       = Lon;%ASSUMING Lon POSITION TO BE TEMP
                      
                      %CONVERT TEMPERATURE(T) IN Celcius(°C) to Kelvin(K) 
                      T = Temp + 273.15; %Temperature in kelvin(K)
                  
                      %CHECK COORD ENTRY FORMAT
                      nCOL_lat=size(Lat,2); %Get number of Columns entry
                   
                      switch nCOL_lat
                       
                          case {7,6}
                              lat = dms2degrees(Lat(:,1:3));%Latitude in DMS (columns 1-3) converted to degrees 
                              lon = dms2degrees(Lat(:,4:6));%Longitude in DMS (columns 4-6) converted to degrees 
                               
                          case {5,4}
                              lat = dm2degrees(Lat(:,1:2));%Latitude in DM (columns 1-2) converted to degrees 
                              lon = dm2degrees(Lat(:,3:4));%Longitude in DM (columns 3-4) converted to degrees   
                                 
                          case {3,2}
                              lat = Lat(:,1);%Latitude in degrees(columns 1) 
                              lon = Lat(:,2);%Longitude in degrees(columns 2)
 
                      end   %switch nCOL_lat
                      
                     end %//if all([size(Temp,1) ~= 131072,size(Temp,1) ~= 8])
                  
                  catch %INCASE Temp ISN'T A 'struct' AND ERRORS ARE DISPLAYED.
                        %IN ORDER TO AVOID THOSE ERROR MSGs, WE USE THE TRY
                        %& CATCH FUNCTION
                      
                       %CONVERT TEMPERATURE(T) IN Celcius(°C) to Kelvin(K) 
                       T          = Temp + 273.15; %Temperature in kelvin(K) 
                       gridV_TVGG = [];%Assigning empty matrix ([]) to gridV_TVGG
                     
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
          
                       %*****CHECK LONGITUDE INPUT 
                       ncol_lon = size(Lon,2); %Get number of longitude entry 
   
                       switch ncol_lon
                           case  3     %if input is in DMS,Convert to degrees
                               lon = dms2degrees(Lon);%Lat in degrees      
                           case  2   %if input is in DM,Convert to degrees    
                               lon = dm2degrees(Lon);%Lat in degrees
                   
                           otherwise     
                                    lon = Lon;
                       end 
                       
                  end %//try
                  
              
           elseif nargin == 3 %IF INPUTS ARE THREE(3) 
                  
                  Temp = Lon;%ASSUMING Lon POSITION TO BE TEMP
                 
                  %CONVERT TEMPERATURE(T) IN Celcius(°C) to Kelvin(K) 
                  T = Temp + 273.15; %Temperature in kelvin(K)
                
                  gridV_TVGG = [];%Assigning empty matrix ([]) to gridV_TVGG
                 
                  %CHECK COORD ENTRY FORMAT
                   nCOL_lat=size(Lat,2); %Get number of Columns entry
                   
                   switch nCOL_lat
                       
                       case {7,6}
                           lat = dms2degrees(Lat(:,1:3));%Latitude in DMS (columns 1-3) converted to degrees 
                           lon = dms2degrees(Lat(:,4:6));%Longitude in DMS (columns 4-6) converted to degrees 
                               
                       case {5,4}
                           lat = dm2degrees(Lat(:,1:2));%Latitude in DM (columns 1-2) converted to degrees 
                           lon = dm2degrees(Lat(:,3:4));%Longitude in DM (columns 3-4) converted to degrees   
                                 
                      case {3,2}
                           lat = Lat(:,1);%Latitude in degrees(columns 1) 
                           lon = Lat(:,2);%Longitude in degrees(columns 2)
 
                   end  %switch nCOL_lat
                   
            end %//if nargin == 5
            
        end %if (any(nargin==[5,4,3])) 
        
    otherwise 
             beep
             fprintf('\n\nInsuficient Inputs for computation of Tm. Check inputs & Try again.\n\n')     
             Tm = [];
             return
        
end  %switch nargin

%===============================END OF INPUT REFORMATTING
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

%********COMPUTE Tm USING Time-Varying Global Gridded(TVGG)MODEL
%Call the "get_TVGG_ERAI_Tm.m" Function
[Tm] = get_TVGG_ERAI_Tm(UTCtime,lat,lon,T,gridV_TVGG);


%=======================VARIOUS SUB-ROUTINES TO COMPUTE Tm=================
%                        ---------------------------------
function [ Tm ] = get_TVGG_ERAI_Tm(UTCtime,lat,lon,T,gridV_TVGG)

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
          Tm(i,1) = TVGG_ERAI_Tm(UTCtime(i,:),lat(i),lon(i),T(i),gridV_TVGG); 
          
      end
      
  else  
       %*****INITIALIZE OUTPUTs 
       Tm = zeros(nrow_pos,nrow_time);
      
     for i = 1:nrow_time %LOOP OVER TIME
         
        for j = 1:nrow_pos %LOOP OVER POSITIONS
            
            %Call the "PI_pwv.m" Function  
            Tm(j,i) = TVGG_ERAI_Tm(UTCtime(i,:),lat(j),lon(j),T(j),gridV_TVGG);
            
        end 
        
     end 
     
  end %//if isequal(Identical,1)
  
else 
    %*****INITIALIZE OUTPUTs 
    Tm = zeros(nrow_pos,nrow_time);
    
    for i = 1:nrow_time %LOOP OVER TIME
        
        for j = 1:nrow_pos %LOOP OVER POSITIONS
              
            %Call the "PI_pwv.m" Function  
            Tm(j,i) = TVGG_ERAI_Tm(UTCtime(i,:),lat(j),lon(j),T(j));
            
        end
        
    end
    
end %//if isequal(nrow_time,nrow_pos)



function [Tm] = TVGG_ERAI_Tm(UTCtime,lat,lon,T,gridV_TVGG)

%========STEPS TO COMPUTE MEAN TEMPERATURE(Tm) USING TVGG MODEL============
%        ------------------------------------------------------
%NOTE: TVGG - Time-Varying Global Gridded Ts-Tm model

%(1)**************************EXTRACTION OF Tm COEFFICIENTS
%CONCATENATE LAT LON COORDS
POS = [lon lat];

%GET THE # ROWS[i.e. #SITEs)
sitecount = size(POS,1);

%****DEFINE CONSTANTs 
R_Earth = 6378.137e3; %Radius Earth [m]; WGS-84 

ERAI_NodeTotel = 512*256;% Totol grid number in ERA-Interim's resolution

ERAI_LonGridSpan = 360/512;%Longtitude interval in ERA-Interim's resolution

ERAI_XNum = 512;%Grid number along latitude line in ERA-Interim's resolution

%******************LATITUDE SERIES IN N128 GAUSSIAN GRID 
N128_Gaussian_Grid_Lat = [89.4628220000000,88.7669510000000,88.0669720000000,87.3660630000000,86.6648030000000,85.9633720000000,85.2618460000000,84.5602610000000,83.8586380000000,83.1569880000000,82.4553190000000,81.7536350000000,...
                          81.0519400000000,80.3502370000000,79.6485270000000,78.9468110000000,78.2450910000000,77.5433670000000,76.8416400000000,76.1399100000000,75.4381770000000,74.7364430000000,74.0347070000000,73.3329700000000,...
                          72.6312310000000,71.9294910000000,71.2277500000000,70.5260080000000,69.8242650000000,69.1225220000000,68.4207770000000,67.7190330000000,67.0172870000000,66.3155410000000,65.6137950000000,64.9120480000000,...
                          64.2103010000000,63.5085540000000,62.8068060000000,62.1050580000000,61.4033090000000,60.7015600000000,59.9998110000000,59.2980620000000,58.5963130000000,57.8945630000000,57.1928140000000,56.4910640000000,...
                          55.7893140000000,55.0875630000000,54.3858130000000,53.6840620000000,52.9823120000000,52.2805610000000,51.5788100000000,50.8770590000000,50.1753080000000,49.4735570000000,48.7718060000000,48.0700540000000,...
                          47.3683030000000,46.6665510000000,45.9648000000000,45.2630480000000,44.5612960000000,43.8595450000000,43.1577930000000,42.4560410000000,41.7542890000000,41.0525370000000,40.3507850000000,39.6490330000000,...
                          38.9472800000000,38.2455280000000,37.5437760000000,36.8420240000000,36.1402710000000,35.4385190000000,34.7367670000000,34.0350140000000,33.3332620000000,32.6315090000000,31.9297570000000,31.2280040000000,...
                          30.5262520000000,29.8244990000000,29.1227460000000,28.4209940000000,27.7192410000000,27.0174880000000,26.3157360000000,25.6139830000000,24.9122300000000,24.2104770000000,23.5087250000000,22.8069720000000,...
                          22.1052190000000,21.4034660000000,20.7017130000000,19.9999600000000,19.2982070000000,18.5964550000000,17.8947020000000,17.1929490000000,16.4911960000000,15.7894430000000,15.0876900000000,14.3859370000000,...
                          13.6841840000000,12.9824310000000,12.2806780000000,11.5789250000000,10.8771720000000,10.1754190000000,9.47366600000000,8.77191300000000,8.07016000000000,7.36840700000000,6.66665400000000,5.96490100000000,...
                          5.26314800000000,4.56139500000000,3.85964200000000,3.15788900000000,2.45613600000000,1.75438300000000,1.05263000000000,0.350877000000000,-0.350877000000000,-1.05263000000000,-1.75438300000000,-2.45613600000000,...
                          -3.15788900000000,-3.85964200000000,-4.56139500000000,-5.26314800000000,-5.96490100000000,-6.66665400000000,-7.36840700000000,-8.07016000000000,-8.77191300000000,-9.47366600000000,-10.1754190000000,-10.8771720000000,...
                          -11.5789250000000,-12.2806780000000,-12.9824310000000,-13.6841840000000,-14.3859370000000,-15.0876900000000,-15.7894430000000,-16.4911960000000,-17.1929490000000,-17.8947020000000,-18.5964550000000,-19.2982070000000,...
                          -19.9999600000000,-20.7017130000000,-21.4034660000000,-22.1052190000000,-22.8069720000000,-23.5087250000000,-24.2104770000000,-24.9122300000000,-25.6139830000000,-26.3157360000000,-27.0174880000000,-27.7192410000000,...
                          -28.4209940000000,-29.1227460000000,-29.8244990000000,-30.5262520000000,-31.2280040000000,-31.9297570000000,-32.6315090000000,-33.3332620000000,-34.0350140000000,-34.7367670000000,-35.4385190000000,-36.1402710000000,...
	                      -36.8420240000000,-37.5437760000000,-38.2455280000000,-38.9472800000000,-39.6490330000000,-40.3507850000000,-41.0525370000000,-41.7542890000000,-42.4560410000000,-43.1577930000000,-43.8595450000000,-44.5612960000000,...
	                      -45.2630480000000,-45.9648000000000,-46.6665510000000,-47.3683030000000,-48.0700540000000,-48.7718060000000,-49.4735570000000,-50.1753080000000,-50.8770590000000,-51.5788100000000,-52.2805610000000,-52.9823120000000,...
	                      -53.6840620000000,-54.3858130000000,-55.0875630000000,-55.7893140000000,-56.4910640000000,-57.1928140000000,-57.8945630000000,-58.5963130000000,-59.2980620000000,-59.9998110000000,-60.7015600000000,-61.4033090000000,...
	                      -62.1050580000000,-62.8068060000000,-63.5085540000000,-64.2103010000000,-64.9120480000000,-65.6137950000000,-66.3155410000000,-67.0172870000000,-67.7190330000000,-68.4207770000000,-69.1225220000000,-69.8242650000000,...
	                      -70.5260080000000,-71.2277500000000,-71.9294910000000,-72.6312310000000,-73.3329700000000,-74.0347070000000,-74.7364430000000,-75.4381770000000,-76.1399100000000,-76.8416400000000,-77.5433670000000,-78.2450910000000,...
	                      -78.9468110000000,-79.6485270000000,-80.3502370000000,-81.0519400000000,-81.7536350000000,-82.4553190000000,-83.1569880000000,-83.8586380000000,-84.5602610000000,-85.2618460000000,-85.9633720000000,-86.6648030000000,...
	                      -87.3660630000000,-88.0669720000000,-88.7669510000000,-89.4628220000000];

%INTERPOLATION
NodeID(sitecount,4) = 0;%four neighboring grids' serial number

NodeInterpC(sitecount,4) = 0;%four neighboring grids' interpolation coefficients

%********LOOP OVER EACH POSITION GIVEN
%At each given site, determine the interpolation coefficients of each 
%neighboring grid to site's position
for i = 1 : sitecount
    lon = POS(i,1);
    lat = POS(i,2);
    
    %ONLY +VE LONGITUDE IN [degrees]
    if lon < 0
        lon = lon + 360;
    end
    
    %Latitude > 89.462822
    if lat > N128_Gaussian_Grid_Lat(1) 
        NodeID(i,:) = 1;
        NodeInterpC(i,1) = 1;
        
    %Latitude < -89.462822  
    elseif lat < N128_Gaussian_Grid_Lat(256) 
        NodeID(i,:) = ERAI_NodeTotel;
        NodeInterpC(i,1) = 1;  
        
    else  % -89.462822 < Latitude < 89.462822
        iLon_grid = floor(lon/ERAI_LonGridSpan);
        if iLon_grid > 511
            iLon_grid = 511;
        end
        
        iLat_grid = 1;
        while N128_Gaussian_Grid_Lat(iLat_grid) > lat
            iLat_grid = iLat_grid + 1;
        end
        
        cosB = cosd(lat);
        sinB = sind(lat);
        
        %  p3   p4
        %  p1   p2       
        cosPsi = sind(N128_Gaussian_Grid_Lat(iLat_grid-1))*sinB...
                 + cosd(N128_Gaussian_Grid_Lat(iLat_grid-1))*cosB*cosd(ERAI_LonGridSpan * (iLon_grid+1) - lon);
        
          psi = acos(cosPsi);
          
           w4 = 1 / (R_Earth * psi);
        
        cosPsi = sind(N128_Gaussian_Grid_Lat(iLat_grid-1))*sinB...
                 + cosd(N128_Gaussian_Grid_Lat(iLat_grid-1))*cosB*cosd(ERAI_LonGridSpan * iLon_grid - lon);
        
          psi = acos(cosPsi);
          
           w3 = 1 / (R_Earth * psi);
        
        cosPsi = sind(N128_Gaussian_Grid_Lat(iLat_grid))*sinB...
                 + cosd(N128_Gaussian_Grid_Lat(iLat_grid))*cosB*cosd(ERAI_LonGridSpan * (iLon_grid+1) - lon);
           
           psi = acos(cosPsi);
           
            w2 = 1 / (R_Earth * psi);
        
        cosPsi = sind(N128_Gaussian_Grid_Lat(iLat_grid))*sinB...
                 + cosd(N128_Gaussian_Grid_Lat(iLat_grid))*cosB*cosd(ERAI_LonGridSpan * iLon_grid - lon);
        
          psi = acos(cosPsi);
          
           w1 = 1 / (R_Earth * psi);
        
            w = w1 + w2 + w3 + w4;
        
        NodeID(i,1) = (iLat_grid-1) * ERAI_XNum + iLon_grid+1;
        NodeInterpC(i,1) = w1 / w;
        
        NodeID(i,3) = (iLat_grid-2) * ERAI_XNum + iLon_grid+1; 
        NodeInterpC(i,3) = w3 / w;
        
        if 511 == iLon_grid 
            NodeID(i,2) = (iLat_grid-1)*ERAI_XNum+1;
            NodeInterpC(i,2) = w2/w;
            
            NodeID(i,4) = (iLat_grid-2)*ERAI_XNum+1;
            NodeInterpC(i,4) = w4/w;
        else
            NodeID(i,2) = (iLat_grid-1)*ERAI_XNum+iLon_grid+2;
            NodeInterpC(i,2) = w2/w;
            
            NodeID(i,4) = (iLat_grid-2)*ERAI_XNum+iLon_grid+2;
            NodeInterpC(i,4) = w4/w;
        end
        
        clear w1 w2 w3 w4 w cosPsi psi cosB sinB iLat_grid iLon_grid;        
    end    
    clear lon lat;
end
clear  N128_Gaussian_Grid_Lat;

%********GET FUNCTION COEFFICIENTS OF ALL 512*256 ERA-INTERIM GRIDS
%CHECK IF gridV_TVGG IS NOT EMPTY([]) 
if ~isempty(gridV_TVGG)
    
   Coeff_TVGG = gridV_TVGG;%ASSIGN gridV_TVGG TO Coeff_TVGG_ERAI
  
else %If gridV_TVGG IS EMPTY([]),read  Coeff_TVGG_ERAI.mat FILE 
    
     %SEARCH FOR TVGG-Tm model GRID FILE(Coeff_TVGG_ERAI.mat) & DIRECTORY & LOAD VALUES
     %Call the "SearchTVGGgrid.m" fxn
     [~,Coeff_TVGG] = SearchTVGGgrid();
     
end %//if ~isempty(gridV_TVGG)   

%INTERPOLATE FUNCTION COEFFICIENTS TO SITES' POSITIONS 
for j = 1 : sitecount
    C_temp = Coeff_TVGG(NodeID(j,:),:);    
    Coeffs(:,j) = C_temp'*NodeInterpC(j,:)';
    clear C_temp;
end

Coeffs = Coeffs';

%(2)*********COMPUTE WEIGHTED MEAN TEMPERATURE(Tm)
%            -------------------------------------
%(2.1)****GET VARIOUS TIME INPUTs 

%(2.1.1)******DATE[YEAR MONTH DAY]
Yr  = UTCtime(:,1);%get Hour
Mn  = UTCtime(:,2);%get Month
Day = UTCtime(:,3);%get Day

%(2.1.2)******TIME[HOUR MINUTE SECONDs]
if size(UTCtime,2) == 6
   H    = UTCtime(:,4);%get Hour
   MIN  = UTCtime(:,5);%get Minute
   SECs = UTCtime(:,6);%get Seconds
   
   %CHANGE SECs,MIN,HOUR WHOSE SEC==60
   MIN(SECs==60) = MIN(SECs==60)+1;
   H(MIN==60) = H(MIN==60)+1;
      
elseif size(UTCtime,2) == 3
       H    = zeros(size(UTCtime,1),1);%create Hour of zeros
       MIN  = zeros(size(UTCtime,1),1);%get Minute of zeros
       SECs = zeros(size(UTCtime,1),1);%get Seconds of zeros    
end

%**********COMPUTE DAY OF YEAR(DoY)
%Call the "DayofYear function.m"
 DoY=DayofYear(Yr,Mn,Day);
 
%*********CONVERT TIME IN [Hour Min Sec] TO HOURs (UT)              
UT = H + (MIN / 60) + (SECs / 3600);% Time in Hours
              
%CONCATENATE TEMPERATURE(T),TIME(UT) & DOY
INPUT = [T DoY UT];

%GENERATE COEFFICIENTS IN TIME-VARYING GLOBAL GRIDDED (TVGG) MODEL AT GIVEN 
%POSITION BY HORIZONTALLY INTERPOLATING COEFFICIENTS OF NEIGHBORING ERA-INTERIM GRIDS 
modelfun=@(f,X)f(1)*X(:,1)+f(2)+f(3)*cos(X(:,2)*2*pi/365.25)+f(4)*sin(X(:,2)*2*pi/365.25)+...
           f(5)*cos(X(:,2)*4*pi/365.25)+f(6)*sin(X(:,2)*4*pi/365.25)+f(7)*cos(X(:,3)*2*pi/24)+f(8)*sin(X(:,3)*2*pi/24);
       
%   where f(1) slope constant,
%         f(2) intercpet constant,
%         f(3) and f(4) for annual variations terms,
%         f(5) and f(6) for semiannual variations terms,
%         f(7) and f(8) for diurnal variations terms,
%   and   X is the function input array,
%         X(:,1), first column in array X, surface temperature input series
%         X(:,2), UTC day of year(doy),
%         X(:,3), UTC hour of day .

%CALCULATE Tm AT THE SELECTED SITE 
Tm = modelfun(Coeffs,INPUT);


%A*********SUBROUTINE TO COMPUTE DAY OF YEAR(DoY)
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


%              --------------------------------------------------------
%B**************SUB-ROUTINE FOR CHECKING THE EXISTENCE OF TVGG GRID FILES
%              --------------------------------------------------------
function [TVGG_grid,gridV_TVGG,folder] = SearchTVGGgrid()
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *                                              
%            The "SearchTVGGgrid" subroutine searches Time-Varying Global..+
%            Gridded Mean Temperature(TVGG-Tm)model Grid file in a given...+
%            Directory/Folder as indicated in the sub-routine. Situation...+ 
%            where TVGG-Tm file is not found in the default directory/...  +   
%            folder,it searches recursively through all Sub-Directories ...+  
%            /Folders of the given Directory.TVGG-Tm file is extracted ... + 
%            by looping through all the listed files in the provided folder+
%            or sub-folders.Finally,if TVGG-Tm file is still not found in  + 
%            the default directory/folder and its sub-folders,the search...+
%            for TVGG-Tm file is extended to the current directory.        +           
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%%USAGE:                                                                   +
%       [TVGG_grid,gridV_TVGG,folder] = SearchTVGGgrid()                   +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%NOTE:                                                                     *
%     THE FF FOLDERS/DIRECTORIES ARE SEARCHED FOR TVGG GRID FILE           *
%1.   'TVGG-Tm file' folder; main folder for TVGG GRID                     *
%2.   'TropoGRIDS' folder; main folder for Tropospheric & Tm GRIDs         *
%3.   'data' folder; folder for goGPS data                                 *
%4.   'pwd'; Current folder/Directory. In goGPS, it is the 'goGPS' folder  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%--------------------------------------------------------------------------
%OUTPUTs                                                                   *
%1.   TVGG_grid : Flag to indicate presence of TVGG grid file.1 to indicate*
%                 there is file in the folder & 0 to mean absence of file  *
%2.  gridV_TVGG :  Extracted TVGG grid Values                              *
%3.     folder  : The found directory or folder(File Path) for TVGG grid...* 
%                 file.e.g.:                                               *
%                 'C:...data\TropoGRIDS\TVGG-Tm file\Coeff_TVGG_ERAI.mat   *       
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
% WRITTEN BY: 
%            OSAH SAMUEL, MSC GEOMATIC ENGINEERING (PhD STUDENT)           +
%            Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com                +
%            Phone:+233(0)246137410 / +233(0)509438484                     +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%**************************************************************************+

%1.SEARCH TVGG-Tm GRIDs from 'TVGG-Tm' folder
[TVGG_grid,TVGGgDir] = SearchTVGGGrid('../data/TropoGRIDS/TVGG-Tm file');

if isequal(TVGG_grid,0)%if TVGG grids are not in 'TVGG-Tm file' folder, try the 'TropoGRIDS'
    
   %2.SEARCH TVGG-Tm GRIDs from 'TropoGRIDS' folder 
   [TVGG_grid,TVGGgDir] = SearchTVGGGrid('../data/TropoGRIDS');

   if isequal(TVGG_grid,0)%if TVGG grids are not in 'TropoGRIDS' folder, try the 'data' folder
    
      [TVGG_grid,TVGGgDir] = SearchTVGGGrid('../data');%searching the 'data' folder
   
      if isequal(TVGG_grid,0) %if TVGG grids are not in the 'data' folder, try the 'current' directory
       
         [TVGG_grid,TVGGgDir] = SearchTVGGGrid(pwd);
      
         if isequal(TVGG_grid,0)%if TVGG grids are not in the 'current' directory, then issue a message
          
            beep %Give a beep sound
            errmsg2{1}=sprintf('No TVGG grid file found in goGPS directory.\n');
            errmsg2{2}=sprintf('Please Provide TVGG grid file & Try Again.\n');
            warndlg(errmsg2,'TVGG grid file Error','modal') 
            
            %ASSIGN EMPTY([]) MATRIX
            gridV_TVGG = [];
            folder     = [];
            
            return
            
         else  
             folder = TVGGgDir;%TVGG grid directory/full path [EX.C:\Users\...data\TropoGRIDS\TVGG-Tm file\Coeff_TVGG_ERAI.mat]
      
         end  
      
      else  
          folder = TVGGgDir; %TVGG grid directory/full path
       
      end  
   
   else   
       folder = TVGGgDir; %TVGG grid directory/full path
       
   end 
       
else 
     folder = TVGGgDir; %TVGG grid directory/full path  
    
end  

%            -------------------------------------------------
%************IF A DIRECTORY/FOLDER WITH TVGG GRIDs IS FOUND....
%            -------------------------------------------------
 if ~isempty(folder) %If Folder is not Empty
    
    %****LOAD TVGG-Tm MODEL COEFFICIENTS
    GridV_TVGG = load(folder);%512*256 ERA-INTERIM GRIDS
    gridV_TVGG = GridV_TVGG.Coeff_TVGG_ERAI;
         
 end
 
 %SAVE DIRECTORY & GRID VALUES 
 setappdata(0,'TVGGgDir',folder) 
 setappdata(0,'gridV_TVGG',gridV_TVGG) 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 
 
%            -----------------------------------------------------------
%***********MAIN SUB-ROUTINE FOR CHECKING THE EXISTENCE OF TVGG GRID FILE
%            ------------------------------------------------------------

function [TVGG_grid,Folder] = SearchTVGGGrid(directory)

%ASSIGNMENT
folder = directory;  %folder/Directory suppose to be containing TVGG FILE in goGPS

%****FIRST CHECK IF FOLDER IS EMPTY OR NOT
SizeListing=size(dir(folder)); %Size of folder content

%********CHECK IF FOLDER IS EMPTY OR NOT
if any([SizeListing(1,1)==2,SizeListing(1,1)==0])% size of 2 means folder is empty & 0 means folder doesn't exist
                  
   TVGG_grid = 0; %Flag to indicate the absence of TVGG grid file in directory
   Folder = []; %File directory is empty([]) if Coeff_TVGG_ERAI.mat is not found
   
   return             
      
else  
     %***GET THE LIST OF FOLDER CONTENTS:
     %===(listFiles is a struct array)=== 
     listFiles = dir(folder);%Get List of Folder/current directory Contents
     
end   %//if any([SizeListing(1,1)==2,SizeListing(1,1)==0])
                          
if exist('listFiles','var') %If Folder is not Empty
                  
   %**CHECK IF FOLDER HAS FILES & SUBFOLDERS
   dirIndex = [listFiles.isdir];%#Find the index for directories(zeros(0's) indicates files,& ones(1s) for folders or pointers)
          
    %GET LIST OF FILES
    fileList = {listFiles(~dirIndex).name}'; %Get list of files & remove pointers('.','..')

    %GET LIST OF SUB-FOLDERS/DIRECTORIES
    subDirs={listFiles(dirIndex).name};%Get list of Sub-folders/directories within parent folder 
                                       %NOTE:filenames={listFiles.name} gives the the list files & subfolders 
                                       %subdirs={listFiles(dirIndex).name} gives list of folders & poniter('.','..') OR                                             
                                       %subdirs1 = filenames([listFiles.isdir])
    %GET INDEX OF SUB-DIRECTORIES
    validIndex = ~ismember(subDirs,{'.','..'});%Find index of subdirectories that are not '.' or '..'
          
    if ~isempty(fileList)% If fileList is not empty
              
       %SORT FILES & REMOVE DUPLICATE FILES
       fileList=unique(fileList);
          
       %*****GET INDEXEs FOR THE TVGG GRID FILE(i.e.Integer Val /[])
       TVGGIndex = regexp(fileList,'\w*Coeff_TVGG_ERAI.mat'); %matches any words with Coeff_TVGG_ERAI.mat.  
          
       %************GET THE GRID FILES
       TVGGFiles  = fileList(~cellfun(@isempty,TVGGIndex));%TVGG files(eg:'Coeff_TVGG_ERAI.mat')
       
    else
        TVGGFiles = [];
           
    end  %if ~isempty(fileList)
    
    %IF NO TVGG FILES EXIST IN THE CURRENT FOLDER, CHECK SUB-FOLDERS/              
    %SUB-DIRECTORIES FOR TVGG GRID FILEs
    %                              +                        ==========
     if isempty(fileList)
              
        checkSubfolders = 1; %Flag to search Sub-Directories/folders
           
     elseif isempty(TVGGFiles)
                   
            checkSubfolders = 1;%Flag to search Sub-Directories/folders 
                         
     else
          if ~isempty(TVGGFiles)
          
             checkSubfolders=0;%Flag not to search Sub-Directories/folders
             
             TVGG_grid = 1; %Flag to indicate the presence of TVGG grid file in directory
             
             %***CONSTRUCT THE FULL PATH [File directory [if Coeff_TVGG_ERAI.mat file is found)]       
             Folder = cellfun(@(x) fullfile(folder,x),TVGGFiles,'UniformOutput',false);%Prepend path to files [File directory [if GTVR-Tm.txt file is found)] 
             
             %SORT FILES & REMOVE DUPLICATE FILES
             [TVGGFiles,i_TVGG] = unique(TVGGFiles);%GET TVGG FILES with no repetitions & INDEX(i_TVGG)
             Folder =(Folder(i_TVGG));%FILES WITH PATH
             
             %CONVERT CELL ARRAY TO CHARACTER
             if iscell(Folder)
                Folder=char(Folder);
             end    
             
             if iscell(TVGGFiles)
                TVGGFiles = char(TVGGFiles);
             end    
             
          end
           
     end  %\\if isempty(fileList)
                
     
     if checkSubfolders==1 %if this condition is true,then either files in the current
                           %Directory were not TVGG grid file or there were no files at all.        
                           %System has to look through Sub-Folders or Sub-directories for TVGG Grid File                                                              
                                                                                                                                                                                                                                                                                  
        %***OPEN ALL SUB-FOLDERS & GET TVGG GRID FILE
                
        p=1;%Loop index
        q=1;
                
        for iDir = find(validIndex) % Loop over valid sub-directories
                   
            %***GET THE LIST OF SUB-FOLDER CONTENTS:
            %   ===(sublist is a struct array)===    
            subList = dir(strcat(folder,'\',subDirs{iDir}));
                                                                                            
            %SIZE OF SUB-FOLDER CONTENT
            Sizesublist=size(subList);

            %**IF  SUB-FOLDER NOT EMPTY 
            if Sizesublist(1,1)~=2 % size of 2 means folder is empty 
                        
               %************CHECK FOR ONLY FILES IN THE SUB-FOLDER
               subList([subList.isdir]) = [];%remove non-directories /current and parent directory.
               subIndex = [subList.isdir];%#Find the index for directories  
                      
               %GET LIST OF FILES
               fileLists = {subList(~subIndex).name}';%Get list of files & remove pointers('.','..')                                            
                
               %CONSTRUCT THE FULLFILE PATH USING THE SUB-DIRECTORY & THE FILENAME
               try
                  Filepath=fullfile(strcat(folder,'\',subDirs{iDir},'\',fileLists));%Prepend path to files
           
               catch   
                     Filepath= cellfun(@(x) fullfile(strcat(folder,'\',subDirs{iDir}),x),...  
                                                     fileLists,'UniformOutput',false); %Prepend path to files  
                   
               end   %try
               
               Filelists{p,1}= fileLists; %#ok<*AGROW>
               Filepaths{p,1}=Filepath;
               
               p=p+1;%Update index
                       
                 %COMBINE ALL GRID FILES FROM VARIOUS SUB-FOLDERS  
                 filelists = vertcat(Filelists{:});     
                 filepaths = vertcat(Filepaths{:}) ;
                 
                  if ~isempty(fileLists)   % If fileLists is not empty
            
                     j=1;%Loop index for the TVGG grid file type
                     
                      for i = 1:length(filelists)%****LOOP OVER THE FILES & CHECK FOR TVGG GRID FILE
                                                 %Identify the needed files and extract them separately
                                                 %into different file name(TVGGfiles)
                                                 %Note that dir also lists the directories,so you have to check for them. 
                                                         
                           %***GET NAME OF FILES FROM THE LIST OF FILES FOR EACH LOOP               
                           
                           FileName = filelists{i,1};
                           filepath = filepaths{i,1};

                           %********MATCH ANY WORD WITH Coeff_TVGG_ERAI.mat
                           if regexpi(FileName,'\w*Coeff_TVGG_ERAI.mat') 

                             %When Tested & the file is a TVGG file
                             TVGGfiles{j,1}=FileName; %#ok<*NASGU>
                             TVGGpath{j,1}=filepath;%TVGG file Path in cells
                             
                             j=j+1;%Update index
                                                                           
                           end  %//if regexpi(FileName,'\w*Coeff_TVGG_ERAI.mat') 
                           
                      end %//for i = 1:length(fileLists) 
                                                                                                               
                  end   %//if ~isempty(fileLists)
                      
            else                          
                 %WRITE ERROR MESSAGE
                 erMsg{q,1}=sprintf('The Folder %s and  It''s Sub-Folder, %s , are empty / Contain No TVGG Grid file .\n',folder,subDirs{iDir});
                         
                 q=q+1; %update index
                        
            end   %if Sizesublist(1,1)~=2
                   
        end    %for iDir = find(validIndex)
        
      if ~exist('TVGGfiles','var')
         
         TVGG_grid = 0;%Flag to indicate the absence of TVGG grid file in directory
         
         Folder = []; %File directory is empty([]) if Coeff_TVGG_ERAI.mat file is not found
         
         return
     
      else 
          TVGG_grid = 1; %Flag to indicate the presence of TVGG grid file in directory
          
          Folder = TVGGpath; %if Coeff_TVGG_ERAI.mat file is found
          
          %SORT FILES & REMOVE DUPLICATE FILES
          [TVGGFiles,i_TVGG] = unique(TVGGfiles);%GET TVGG FILES with no repetitions & INDEX(i_TVGG)
           Folder =(Folder(i_TVGG));%FILES WITH PATH
         
          %CONVERT CELL ARRAY TO CHARACTER
          if iscell(Folder)
             Folder=char(Folder);
          end 
          
          if iscell(TVGGFiles)
             TVGGFiles = char(TVGGFiles);
          end     
                    
      end %//~exist('TVGGfiles','var')                                        
           
     end  %//if checkSubfolders==1
     
end %//if exist('listFiles','var')

%=========================================END OF SearchTVGGGrid.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%                        END OF SearchTVGGgrid
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
