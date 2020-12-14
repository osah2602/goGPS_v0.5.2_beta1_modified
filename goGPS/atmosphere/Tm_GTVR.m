%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *
%            "Tm_GTVR" is a subroutine that Computes the weighted-mean...  *
%             temperature(Tm) based on global-vertical-correction-model    * 
%             The model coefficients are stored in the grid points on a    * 
%             global scale (87.5°S?90°N and 180°W?177.5°E) in the file     *
%             GTVR-Tm.txt                                                  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%USAGE:                                                                    *
%      Tm = Tm_GTVR(Time,lat,lon,h,gridV_GTVR,GPT2w_grid)                  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%1.     Time  = Time of observation                                        *
% Time format :[year month day hour minute seconds] OR [year month day]    *                  
%2.       lat = Latitude of Station in degrees  OR [D M S]                 *
%3.       lon = Longitude of Station in degrees OR [D M S]                 *
%4.         h = Height of Station                                          *
%5.gridV_GTVR = GTVR-Tm model grid values extracted from GTVR-Tm.txt file  *
%6.GPT2w_grid = grid values in cells extracted from GPT2w 5°x5° external   *  
%               grid file                                                  *
%NOTE:                                                                     *                                                            
%     Surface Tm0 and/or h0 is/are provided by default using GPT2w model...* 
%     (Böhm et al., 2014).Where  GPT2w_grid is not provided(i.e.empty([])  *
%     either the GTmIII model or Bevis_GPT model is used to compute Tm0    *                                                        *
%     and/or h0.
%WHERE :                                                                   *
%       Tm0 = Surface weighted mean temperature                            *
%        h0 = Surface Height(i.e. Orthometric height for GTmIII & Bevis_GPT*
%             Models but orthometric grid height for GPT2w model)          *

%OUTPUT:
%       Tm = Weighted Mean Temperature                                     *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================
%REFERENCE:                                                                *
%Yao, Y., Sun, Z., Xu, C., Xu, X., & Kong, J. (2018).Extending a model for *
%water vapor sounding by ground-based GNSS in the vertical direction.      * 
%Journal of Atmospheric and Solar-Terrestrial Physics, 179, 358-366.       *
% =========================================================================
%  Code Editted by:                                                        * 
%                  OSAH SAMUEL, Msc Geomatic Engineering ,2016             *    
%                  Email: osahsamuel@yahoo.ca                              *
%                  Tel:+233 (0)246137410/+233 (0)509438484                 * 
%==========================================================================
%***************************************************************************
%***************************************************************************
function Tm = Tm_GTVR(Time,Lat,Lon,h,gridV_GTVR,GPT2w_grid)

%**********CHECK & REFORMAT INPUTs DATA
switch nargin
    
   case {6,5,4,3,2} %Various inputs format 
       
        if (any(nargin==[6,5,4,3,2])) 
            
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
           
            if (any(nargin==[6,5])) %IF 6 or 5 INPUTS ARE PROVIDED
                
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
           

               if nargin == 5 %IF 5 INPUTS ARE PROVIDED
          
                  GPT2w_grid = [];%Assigning empty matrix ([]) to GPT2w_grid
                  
               end
                  
            elseif nargin == 4 %IF 4 INPUTS ARE PROVIDED
                   
                  %CHECK IF h & Lon POSITIONS ARE GTVR-Tm & GPT2w GRID
                  %VALUES RESPECTIVELY(gridV_GTVR,GPT2w_grid)
                  %--------------------------------------------------------
                  %NOTE:
                  %     THIS HAPPENS WHEN USER DECIDE TO PROVIDE LAT,LON & h
                  %     COORD AS A SINGLE INPUT,I.E.[Lat Lon h].IN THAT CASE,
                  %     gridV_GTVR WILL TAKE THE POSITION OF h & GPT2w_grid (Time,Lat,Lon,h
                  %     THAT OF Lon.
                  %--------------------------------------------------------
                  if all([~isempty(Lon),~isempty(h)])
                      
                     if all([all([size(Lon,1) == 10368,size(Lon,2) == 11]),all([length(h) == 10,iscell(h)])])
                      
                        gridV_GTVR = Lon;%ASSUMING Lon POSITION TO BE gridV_GTVR
                        GPT2w_grid = h;%ASSUMING h POSITION TO BE GPT2w_grid
                      
                     elseif all([all([size(h,1) == 10368,size(h,2) == 11]),all([length(Lon) == 10,iscell(Lon)])])
                      
                            gridV_GTVR = h;%ASSUMING h POSITION TO BE gridV_GTVR
                            GPT2w_grid = Lon;%ASSUMING Lon POSITION TO BE GPT2w_grid
                      
                     else 
                         gridV_GTVR = [];%Assigning empty matrix ([]) to gridV_TVGG
                         GPT2w_grid = [];%Assigning empty matrix ([]) to GPT2w_grid
                         
                     end %//if all([all([size(Lon,1) == 10368,size(Lon,2) == 11]),all([length(h) == 10,iscell(h)])])    
                     
                  elseif all([~isempty(Lon),isempty(h)])
                          
                         %CHECK IF Lon POSITIONS IS GTVR-Tm GRID VALUES 
                         if all([size(Lon,1) == 10368,size(Lon,2) == 11])%CHECK DIMs for GTVR-Tm GRID VALUES
                         
                            gridV_GTVR = Lon;%ASSUMING Lon POSITION TO BE gridV_TVGG
                            GPT2w_grid = [];%Assigning empty matrix ([]) to GPT2w_grid
                            
                         %CHECK IF Lon POSITIONS IS GPT2w GRID VALUES    
                         elseif all([length(Lon) == 10,iscell(Lon)])
                             
                                GPT2w_grid = Lon;%ASSUMING Lon POSITION TO BE GPT2w_grid
                                gridV_GTVR = [];%Assigning empty matrix ([]) to gridV_TVGG
                                
                         end %//if all([size(Lon,1) == 10368,size(Lon,2) == 11])
                         
                  elseif all([isempty(Lon),~isempty(h)])
                      
                         %CHECK IF h POSITIONS IS GTVR-Tm GRID VALUES 
                         if all([size(h,1) == 10368,size(h,2) == 11])%CHECK DIMs for GTVR-Tm GRID VALUES
                         
                            gridV_GTVR = h;%ASSUMING h POSITION TO BE gridV_TVGG
                            GPT2w_grid = [];%Assigning empty matrix ([]) to GPT2w_grid
                            
                         %CHECK IF Lon POSITIONS IS GPT2w GRID VALUES    
                         elseif all([length(h) == 10,iscell(h)])
                             
                                GPT2w_grid = h;%ASSUMING h POSITION TO BE GPT2w_grid
                                gridV_GTVR = [];%Assigning empty matrix ([]) to gridV_TVGG
                                
                         end %//if all([size(h,1) == 10368,size(h,2) == 11])
                         
                  else
                      gridV_GTVR = [];%Assigning empty matrix ([]) to gridV_TVGG
                      GPT2w_grid = [];%Assigning empty matrix ([]) to GPT2w_grid
                      
                  end %//if all([~isempty(Lon),~isempty(h)])
                  
                  %-------------------------------------------------------
                  %CHECK IF Lon & h POSITIONS WERE NOT GTVR-Tm & GPT2w GRID
                  %VALUES BUT ACTUAL COORD VALUES
                  %--------------------------------------------------------
                  if all([isempty(gridV_GTVR),isempty(GPT2w_grid)])%IF gridV_GTVR & GPT2w_grid ARE EMPTY([])
                      
                     if all([all([~isempty(Lon),~isempty(h)]),any([all([~all([size(Lon,1) == 10368,size(Lon,2) == 11]),...
                             ~all([length(h) == 10,iscell(h)])]),all([~all([size(h,1) == 10368,size(h,2) == 11]),...
                             ~all([length(Lon) == 10,iscell(Lon)])])])])%CHECK THE SOURCE OF EMPTINESS.IF SOURCE IS FROM...
                                                                        %WHEN Lon & h POSITIONS ARE ACTUAL COORD VALUES
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
                       
                     end %//if all([all([~isempty(Lon),~isempty(h)]),any([all([~all([size(Lon,1) == 10368,size(Lon,2) == 11]),...
                         %          ~all([length(h) == 10,iscell(h)])]),all([~all([size(h,1) == 10368,size(h,2) == 11]),...
                         %          ~all([length(Lon) == 10,iscell(Lon)])])])])
                     
                  else
                       if any([all([~isempty(gridV_GTVR),~isempty(GPT2w_grid)]),any([~isempty(gridV_GTVR),~isempty(GPT2w_grid)])])
                                
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
                          
                       end %//if any([all([~isempty(gridV_GTVR),~isempty(GPT2w_grid)]),any([~isempty(gridV_GTVR),~isempty(GPT2w_grid)])])
                           
                  end  %//if all([isempty(gridV_GTVR),isempty(GPT2w_grid)])
                      
            elseif  (any(nargin==[3,2])) %IF 3 or 2 INPUTS ARE PROVIDED ASSUMING POSITION IS IN [LAT LON h]
                   
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
            
                      end %switch nCOL_lat
                      
                      if nargin == 3 %IF 3 INPUTS ARE PROVIDE
                         
                         %CHECK IF Lon POSITIONS IS GTVR-Tm GRID VALUES 
                         if all([size(Lon,1) == 10368,size(Lon,2) == 11])%CHECK DIMs for GTVR-Tm GRID VALUES
                         
                            gridV_GTVR = Lon;%ASSUMING Lon POSITION TO BE gridV_TVGG
                            GPT2w_grid = [];%Assigning empty matrix ([]) to GPT2w_grid
                            
                         %CHECK IF Lon POSITIONS IS GPT2w GRID VALUES    
                         elseif all([length(Lon) == 10,iscell(Lon)])
                             
                                GPT2w_grid = Lon;%ASSUMING Lon POSITION TO BE GPT2w_grid
                                gridV_GTVR = [];%Assigning empty matrix ([]) to gridV_TVGG
                                
                         else
                             gridV_GTVR = [];%Assigning empty matrix ([]) to gridV_TVGG
                             GPT2w_grid = [];%Assigning empty matrix ([]) to GPT2w_grid
                             
                         end 
                         
                      elseif nargin == 2 %IF 2 INPUTS ARE PROVIDE
                          
                             gridV_GTVR = [];%Assigning empty matrix ([]) to gridV_TVGG
                             GPT2w_grid = [];%Assigning empty matrix ([]) to GPT2w_grid
                             
                      end %//if nargin == 3
                      
            end %if (any(nargin==[6,5]))
            
        end %//if (any(nargin==[6,5,4,3,2]))
        
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
[Tm] = get_GTVR_Tm(UTCtime,lat,lon,h,gridV_GTVR,GPT2w_grid);


%=======================VARIOUS SUB-ROUTINES TO COMPUTE Tm=================
%                        ---------------------------------
function [ Tm ] = get_GTVR_Tm(UTCtime,lat,lon,h,gridV_GTVR,GPT2w_grid)

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
          Tm(i,1) = GTVR_Tm(UTCtime(i,:),lat(i),lon(i),h(i),gridV_GTVR,GPT2w_grid); 
          
      end
      
  else  
       %*****INITIALIZE OUTPUTs 
       Tm = zeros(nrow_pos,nrow_time);
      
     for i = 1:nrow_time %LOOP OVER TIME
         
        for j = 1:nrow_pos %LOOP OVER POSITIONS
            
            %Call the "PI_pwv.m" Function  
            Tm(j,i) = GTVR_Tm(UTCtime(i,:),lat(j),lon(j),h(j),gridV_GTVR,GPT2w_grid);
            
        end 
        
     end 
     
  end %//if isequal(Identical,1)
  
else 
    %*****INITIALIZE OUTPUTs 
    Tm = zeros(nrow_pos,nrow_time);
    
    for i = 1:nrow_time %LOOP OVER TIME
        
        for j = 1:nrow_pos %LOOP OVER POSITIONS
              
            %Call the "PI_pwv.m" Function  
            Tm(j,i) = GTVR_Tm(UTCtime(i,:),lat(j),lon(j),h(j),gridV_GTVR,GPT2w_grid);
            
        end
        
    end
    
end %//if isequal(nrow_time,nrow_pos)

function [ Tm ] = GTVR_Tm(UTCtime,lat,lon,h,gridV_GTVR,GPT2w_grid)

%*****CHECK IF GPT2w_grid IS EMPTY([])
if isempty(GPT2w_grid)
   %RETRIEVE STORED GPT2w GRID VALUES FROM goGPS GUI FIG
   GPT2w_grid = getappdata(0,'GPT2w_gridV') ; 
end

%CHECK IF gridV_GTVR IS NOT EMPTY([]) 
if ~isempty(gridV_GTVR)
    
  c = gridV_GTVR;%ASSIGN gridV_GTVR TO c
  
else %IF gridV_GTVR IS EMPTY([]),read  GTVR-Tm.txt FILE 
     %*****GET GTVR-Tm.txt FILE & DIRECTORY                   
     %Call the "searchGTVRgrid.m" function
     [~,GTVRfile] = SearchGTVRgrid();

     if ~isempty(GTVRfile) 
    
        %**IMPORT FILE & ASSIGN TO c
        c = importdata(GTVRfile);% FILE IS 10368x11
        
     end
     
end %//if ~isempty(gridV_GTVR)

%******COMPUTE DAY OF YEAR(doy) 
%***CALL THE DayofYear.m FUNCTION
 doy=DayofYear(UTCtime);
 
%GET THE LATITUDE AND LONGITUDE OF THE FOUR CORNER POINTS OF THE GRID 
%POINT WHERE THE LATITUDE AND LONGITUDE OF THE TARGET POINT IS LOCATED
B1 = floor(lat/2.5)*2.5;
B2 = ceil(lat/2.5)*2.5;
L1 = floor(lon/2.5)*2.5;
L2 = ceil(lon/2.5)*2.5;

if B1==B2&&L1==L2
    
    a = c(c(:,1)==lat&c(:,2)==lon,:);
    
elseif B1==B2
    
    a1 = c(c(:,1)==B1&c(:,2)==L1,:);
    a2 = c(c(:,1)==B1&c(:,2)==L2,:);
    a = ((L2-lon)*a1+(lon-L1)*a2)/(L2-L1);    
    
elseif L1==L2
    
    a1 = c(c(:,1)==B1&c(:,2)==L1,:);
    a2 = c(c(:,1)==B2&c(:,2)==L1,:);
    a = ((B2-lat)*a1+(lat-B1)*a2)/(B2-B1);   
    
else
    %****GET THE COEFFICIENTS OF THE FOUR GRID POINTS
    a1 = c(c(:,1)==B1&c(:,2)==L1,:);
    a2 = c(c(:,1)==B1&c(:,2)==L2,:);
    a3 = c(c(:,1)==B2&c(:,2)==L1,:);
    a4 = c(c(:,1)==B2&c(:,2)==L2,:);

    %***BILINEAR INTERPOLATION TO OBTAIN TARGET POINT COEFFICIENTS
    a = ((B2-lat)*(L2-lon)*a1+(lat-B1)*(lon-L1)*a4+(lat-B1)*(L2-lon)*a3+(B2-lat)*(lon-L1)*a2)/((B2-B1)*(L2-L1));    
end

X = a(:,3:9)';

%COMPUTE SURFACE Tm0 & SURFACE HEIGHT(h0)
if ~isempty(GPT2w_grid)

   Timevar = 0; %APPLY time variation (annual and semiannual terms) 
   [Tm0,h0] = GPT2w_5_x_5(UTCtime,deg2rad(lat),deg2rad(lon),h,GPT2w_grid,Timevar);   

else  
    try
        [Tm0,h0] = GTmIII(UTCtime, lat, lon, h);
        
    catch
         [Tm0,h0] = Bevis_GPT(UTCtime,lat,lon,h);
    
    end
    
end


%CONVERT HEIGHTS IN METERS TO KM
h  = h./1000; %[Km]
h0 = h0./1000;%[Km]

%***COMPUTING DESIGN MATRIX
dh = h-h0;%height difference(km)

len = length(h(:,1));
i = 1;
I = repmat(i,len,1); %#ok<*REPMAT>

C1 = cos(2*pi*doy/365.25);
S1 = sin(2*pi*doy/365.25);
C2 = cos(4*pi*doy/365.25);
S2 = sin(4*pi*doy/365.25);
C3 = cos(2*pi*h/19000)-cos(2*pi*h0/19000);
S3 = sin(2*pi*h/19000)-sin(2*pi*h0/19000);

H = [I.*dh,C1.*dh,S1.*dh,C2.*dh,S2.*dh,C3,S3];

Tm = H*X+Tm0;


%A.SUBROUTINE TO COMPUTE METEOROLOGICAL PARAMETERS FROM 5x5 GRID GPT2w MODEL
function [Tm,Hg,horth] = GPT2w_5_x_5(UTCtime,lat,lon,hell,grid,Timevar)

%**************************************************************************
%DESCRIPTION:
%            This subroutine determines pressure, temperature, temperature *
%            lapse rate, mean temperature of the water vapor, water vapor  * 
%            pressure, hydrostatic and wet mapping function coefficients   * 
%            ah and aw, water vapour decrease factor and geoid undulation  * 
%            for specific sites near the Earth surface. It is based on a   *
%            5 x 5 degree external grid file ('gpt2_5w.grd') with mean     *
%            values as well as sine and cosine amplitudes for the annual   * 
%            and semiannual variation of the coefficients.
%USAGE:
%      [Tm,Hg,horth] = GPT2w_5(UTCtime,lat,lon,hell,Timevar)               *
%INPUTs:                                                                   *
%1.   UTCtime:  UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *                                                                
%2.       lat:  Station geodetic latitude in radians [-pi/2:+pi/2] (n x 1) *
%3.       lon:  Station geodetic longitude in radians [-pi:pi] or [0:2pi]  *
%4.      hell:  ellipsoidal height in [meters] (n x 1)                     *
%5    Timevar:  case 1: no time variation but static quantities            *
%               case 0: with time variation (annual and semiannual terms)  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% OUTPUTs:
%1.      Tm: mean temperature of the water vapor in degrees Kelvin         * 
%2.      Hg: Grid orthometric Height(i.e ho FOR THE GTVR MODEL )           *
%2.   horth: Orthometric Height                                            *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%%REFERENCE:                                                               *
%           J. Böhm, G. Möller, M. Schindelegger, G. Pain, R. Weber,...    *
%           Development of an improved blind model for slant delays in the *
%           troposphere (GPT2w),GPS Solutions, 2014,...                    *
%           doi:10.1007/s10291-014-0403-7                                  *

%Original codes by Böhm et al 2014                                         *
%Modified by: Samuel Osah,Msc Geomatic Engineering ,2016                   *    
%             Email: osahsamuel@yahoo.ca                                   *
%             Tel:+233 (0)246137410/+233 (0)509438484         
%==========================================================================+
%**************************************************************************+
%**************************************************************************+
%********COMPUTE MODIFIED JULIAN DATE(MJD)
%GET UTC TIME COMPONENTs
%1.******UTCdate
  Yr = UTCtime(:,1);%get Year
  Mn = UTCtime(:,2);%get Month
 Day = UTCtime(:,3);%get Day
%2.UTCtime
   H = UTCtime(:,4);%get Hour
 MIN = UTCtime(:,5);%get Minute
SECs = UTCtime(:,6);%get Seconds

%****MODIFIED JULIAN DATE
[~, dmjd ,~]=utc2JulianDay_DoY(Yr,Mn,Day,H,MIN,SECs);

%change the reference epoch to January 1 2000
dmjd1 = dmjd-51544.5;

%*******DEFINE CONSTANTs
gm = 9.80665;%Mean Gravity in m/s^2
dMtr = 28.965*10^-3;%molar mass of dry air in kg/mol
Rg = 8.3143;% universal gas constant in J/K/mol

%Factors for amplitudes
if (Timevar==1) % then  constant parameters
    cosfy = 0;
    coshy = 0;
    sinfy = 0;
    sinhy = 0;
else 
    cosfy = cos(dmjd1/365.25*2*pi);
    coshy = cos(dmjd1/365.25*4*pi);
    sinfy = sin(dmjd1/365.25*2*pi);
    sinhy = sin(dmjd1/365.25*4*pi);
end

%READ RESPECTIVE DATA FROM grid 
p_grid    = grid{1};    % pressure in Pascal
T_grid    = grid{2};    % temperature in Kelvin
Q_grid    = grid{3};    % specific humidity in kg/kg
dT_grid   = grid{4};    % temperature lapse rate in Kelvin/m
u_grid    = grid{5};    % geoid undulation in m
Hs_grid   = grid{6};    % orthometric grid height in m
ah_grid   = grid{7};    % hydrostatic mapping function coefficient, dimensionless
aw_grid   = grid{8};   % wet mapping function coefficient, dimensionless
la_grid   = grid{9};   % water vapor decrease factor, dimensionless
Tm_grid   = grid{10};   % mean temperature in Kelvin

%DETERMINE NUMBER OF STATIONS
nstat = length(lat);%Number of stations

%INITIALIZATION OF NEW VECTORS
p     =  zeros([nstat, 1]);
T     =  zeros([nstat, 1]);
dT    = zeros([nstat, 1]);
Tm    = zeros([nstat, 1]);
e     =  zeros([nstat, 1]);
ah    = zeros([nstat, 1]);
aw    = zeros([nstat, 1]);
la    = zeros([nstat, 1]);
undu  = zeros([nstat, 1]);
Hg    = zeros([nstat, 1]);
horth = zeros([nstat, 1]);

%LOOP OVER STATIONs
for k = 1:nstat
    
    %ONLY +VE LONGITUDE IN [degrees]
    if lon(k) < 0
        plon = (lon(k) + 2*pi)*180/pi;
    else
        plon = lon(k)*180/pi;
    end
    %TRANSFORM TO POLAR DISTANCE IN [degrees]
    ppod = (-lat(k) + pi/2)*180/pi; 

    %FIND THE INDEX (LINE IN THE GRID FILE) OF THE NEAREST POINT
    ipod = floor((ppod+5)/5); 
    ilon = floor((plon+5)/5);
    
    %NORMALIZED (TO ONE) DIFFERENCES, CAN BE +VE or -VE
    diffpod = (ppod - (ipod*5 - 2.5))/5;
    difflon = (plon - (ilon*5 - 2.5))/5;
    
    %added by HCY
    if ipod == 37
        ipod = 36;
    end
    
	%added by GP
    if ilon == 73
		ilon = 1;
    end
    if ilon == 0
		ilon = 72;
    end 

    %GET THE NUMBER OF THE CORRESPONDING LINE
    indx(1) = (ipod - 1)*72 + ilon;
    
    %Near the poles: nearest neighbour interpolation, otherwise: bilinear
    bilinear = 0;
    if ppod > 2.5 && ppod < 177.5 
           bilinear = 1;          
    end          
    
    %Case of nearest neighbourhood
    if bilinear == 0

        ix = indx(1);
        
        %Transforming ellipsoidial height to orthometric height
        undu(k) = u_grid(ix);%Geoid Undulation
       horth(k) = hell(k)-undu(k);%Orthomethric Height
          Hg(k) = Hs_grid(ix);%Grid Orthomethric Height
        
        %Pressure, temperature at the heigtht of the grid
        T0 = T_grid(ix,1) + ...
             T_grid(ix,2)*cosfy + T_grid(ix,3)*sinfy + ...
             T_grid(ix,4)*coshy + T_grid(ix,5)*sinhy;
        p0 = p_grid(ix,1) + ...
             p_grid(ix,2)*cosfy + p_grid(ix,3)*sinfy+ ...
             p_grid(ix,4)*coshy + p_grid(ix,5)*sinhy;
         
        %Specific humidity
        Q = Q_grid(ix,1) + ...
            Q_grid(ix,2)*cosfy + Q_grid(ix,3)*sinfy+ ...
            Q_grid(ix,4)*coshy + Q_grid(ix,5)*sinhy;
            
        %lapse rate of the temperature
        dT(k) = dT_grid(ix,1) + ...
                dT_grid(ix,2)*cosfy + dT_grid(ix,3)*sinfy+ ...
                dT_grid(ix,4)*coshy + dT_grid(ix,5)*sinhy; 

               
        %Station height - grid height
        redh = horth(k) - Hg(k);

        %Temperature at station height in Celsius
        T(k) = T0 + dT(k)*redh - 273.15;
        
        %Temperature lapse rate in degrees / km
        dT(k) = dT(k)*1000;

        %Virtual temperature in Kelvin
        Tv = T0*(1+0.6077*Q);       
        c = gm*dMtr/(Rg*Tv);
        
        %Pressure in hPa
        p(k) = (p0*exp(-c*redh))/100;
            
        %Hydrostatic coefficient ah 
        ah(k) = ah_grid(ix,1) + ...
                ah_grid(ix,2)*cosfy + ah_grid(ix,3)*sinfy+ ...
                ah_grid(ix,4)*coshy + ah_grid(ix,5)*sinhy;
            
        %Wet coefficient aw
        aw(k) = aw_grid(ix,1) + ...
                aw_grid(ix,2)*cosfy + aw_grid(ix,3)*sinfy + ...
                aw_grid(ix,4)*coshy + aw_grid(ix,5)*sinhy;
		
		%Water vapour decrease factor la - added by GP
        la(k) = la_grid(ix,1) + ...
                la_grid(ix,2)*cosfy + la_grid(ix,3)*sinfy + ...
                la_grid(ix,4)*coshy + la_grid(ix,5)*sinhy; 
		
		%Mean temperature Tm - added by GP
        Tm(k) = Tm_grid(ix,1) + ...
                Tm_grid(ix,2)*cosfy + Tm_grid(ix,3)*sinfy + ...
                Tm_grid(ix,4)*coshy + Tm_grid(ix,5)*sinhy;
		
		%Water vapor pressure in hPa - changed by GP
		e0 = Q*p0/(0.622+0.378*Q)/100; % on the grid
		e(k) = e0*(100*p(k)/p0)^(la(k)+1);   % on the station height - (14) Askne and Nordius, 1987
		
     else % bilinear interpolation
        
        ipod1 = ipod + sign(diffpod);
        ilon1 = ilon + sign(difflon);
        if ilon1 == 73
            ilon1 = 1;
        end
        if ilon1 == 0
            ilon1 = 72;
        end
        
        %Get the number of the line
        indx(2) = (ipod1 - 1)*72 + ilon;  % along same longitude
        indx(3) = (ipod  - 1)*72 + ilon1; % along same polar distance
        indx(4) = (ipod1 - 1)*72 + ilon1; % diagonal
        
        for l = 1:4
               
            %Transforming ellipsoidial height to orthometric height:
            % Hortho = -N + Hell
            undul(l) = u_grid(indx(l));
            horth1(l) = hell(k)-undul(l);%Orthomethric Height
            Hs1(l) = Hs_grid(indx(l));%Grid Orthomethric Height
            
            %Pressure, Temperature at the heigtht of the grid
            T0 = T_grid(indx(l),1) + ...
                 T_grid(indx(l),2)*cosfy + T_grid(indx(l),3)*sinfy + ...
                 T_grid(indx(l),4)*coshy + T_grid(indx(l),5)*sinhy;
            p0 = p_grid(indx(l),1) + ...
                 p_grid(indx(l),2)*cosfy + p_grid(indx(l),3)*sinfy + ...
                 p_grid(indx(l),4)*coshy + p_grid(indx(l),5)*sinhy;

            %Humidity 
            Ql(l) = Q_grid(indx(l),1) + ...
                    Q_grid(indx(l),2)*cosfy + Q_grid(indx(l),3)*sinfy + ...
                    Q_grid(indx(l),4)*coshy + Q_grid(indx(l),5)*sinhy;
 
            %Reduction = stationheight - gridheight
            redh = horth1(l) - Hs1(l);

            %lapse rate of the temperature in degree / m
            dTl(l) = dT_grid(indx(l),1) + ...
                     dT_grid(indx(l),2)*cosfy + dT_grid(indx(l),3)*sinfy + ...
                     dT_grid(indx(l),4)*coshy + dT_grid(indx(l),5)*sinhy; 

            %Temperature reduction to station height
            Tl(l) = T0 + dTl(l)*redh - 273.15;

            %Virtual temperature
            Tv = T0*(1+0.6077*Ql(l));  
            c = gm*dMtr/(Rg*Tv);
            
            %Pressure in hPa
            pl(l) = (p0*exp(-c*redh))/100;
            
            %Hydrostatic coefficient ah
            ahl(l) = ah_grid(indx(l),1) + ...
                     ah_grid(indx(l),2)*cosfy + ah_grid(indx(l),3)*sinfy + ...
                     ah_grid(indx(l),4)*coshy + ah_grid(indx(l),5)*sinhy;
            
            %Wet coefficient aw
            awl(l) = aw_grid(indx(l),1) + ...
                     aw_grid(indx(l),2)*cosfy + aw_grid(indx(l),3)*sinfy + ...
                     aw_grid(indx(l),4)*coshy + aw_grid(indx(l),5)*sinhy;
					 
			%Water vapour decrease factor la - added by GP
			lal(l) = la_grid(indx(l),1) + ...
					 la_grid(indx(l),2)*cosfy + la_grid(indx(l),3)*sinfy + ...
					 la_grid(indx(l),4)*coshy + la_grid(indx(l),5)*sinhy;
			
			%Mean temperature Tm - added by GP
			Tml(l) = Tm_grid(indx(l),1) + ...
					 Tm_grid(indx(l),2)*cosfy + Tm_grid(indx(l),3)*sinfy + ...
					 Tm_grid(indx(l),4)*coshy + Tm_grid(indx(l),5)*sinhy;
            
			%Water vapor pressure in hPa - changed by GP
			e0 = Ql(l)*p0/(0.622+0.378*Ql(l))/100; % on the grid
			el(l) = e0*(100*pl(l)/p0)^(lal(l)+1);  % on the station height - (14) Askne and Nordius, 1987			
        end            
        dnpod1 = abs(diffpod); % distance nearer point
        dnpod2 = 1 - dnpod1;   % distance to distant point
        dnlon1 = abs(difflon);
        dnlon2 = 1 - dnlon1;
        
        %Pressure
        R1 = dnpod2*pl(1)+dnpod1*pl(2);
        R2 = dnpod2*pl(3)+dnpod1*pl(4);
        p(k) = dnlon2*R1+dnlon1*R2;
            
        %Temperature
        R1 = dnpod2*Tl(1)+dnpod1*Tl(2);
        R2 = dnpod2*Tl(3)+dnpod1*Tl(4);
        T(k) = dnlon2*R1+dnlon1*R2;
        
        %Temperature in degree per km
        R1 = dnpod2*dTl(1)+dnpod1*dTl(2);
        R2 = dnpod2*dTl(3)+dnpod1*dTl(4);
        dT(k) = (dnlon2*R1+dnlon1*R2)*1000;
            
        %Water vapor pressure in hPa - changed by GP
		R1 = dnpod2*el(1)+dnpod1*el(2);
        R2 = dnpod2*el(3)+dnpod1*el(4);
        e(k) = dnlon2*R1+dnlon1*R2;
            
        %Hydrostatic
        R1 = dnpod2*ahl(1)+dnpod1*ahl(2);
        R2 = dnpod2*ahl(3)+dnpod1*ahl(4);
        ah(k) = dnlon2*R1+dnlon1*R2;
           
        %Wet
        R1 = dnpod2*awl(1)+dnpod1*awl(2);
        R2 = dnpod2*awl(3)+dnpod1*awl(4);
        aw(k) = dnlon2*R1+dnlon1*R2;
        
        %undulation
        R1 = dnpod2*undul(1)+dnpod1*undul(2);
        R2 = dnpod2*undul(3)+dnpod1*undul(4);
        undu(k) = dnlon2*R1+dnlon1*R2;
		
        %Grid Orthomethric Height
        R1 = dnpod2*Hs1(1)+dnpod1*Hs1(2);
        R2 = dnpod2*Hs1(3)+dnpod1*Hs1(4);
        Hg(k) = dnlon2*R1+dnlon1*R2;
        
        %Orthomethric Height
        R1 = dnpod2*horth1(1)+dnpod1*horth1(2);
        R2 = dnpod2*horth1(3)+dnpod1*horth1(4);
        horth(k) = dnlon2*R1+dnlon1*R2;
        
		%Water vapour decrease factor - added by GP
        R1 = dnpod2*lal(1)+dnpod1*lal(2);
        R2 = dnpod2*lal(3)+dnpod1*lal(4);
        la(k) = dnlon2*R1+dnlon1*R2;
		
		%Mean temperature Tm - added by GP
        R1 = dnpod2*Tml(1)+dnpod1*Tml(2);
        R2 = dnpod2*Tml(3)+dnpod1*Tml(4);
        Tm(k) = dnlon2*R1+dnlon1*R2;
                    
    end %//if bilinear == 0
    
end%//for k = 1:nstat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GPT2w_5 x 5.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%B.******SUBROUTINE TO COMPUTE SURFACE MEAN TEMPERATURE(Tm0)BY GTm-III MODEL
function [Tm,horth] = GTmIII(UTCtime,lat,lon,h)

%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *
%            "GTmIII" is a subroutine estimates Weighted-Surface Mean      *
%             Temperature(Tm0) based on Spherical Harmonics up to degree   * 
%             and order  9                                                 *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================
%REFERENCE:                                                                *
%YAO YiBin, XUChaoQian,ZHANGBao, CAONa                                     *
%GTm-III: A new Global Empirical Model for Mapping Zenith Wet Delays onto  *
%Precipitable Water // Vapor,to be submitted to Journal of Geodesy,2013.8  *
%==========================================================================
%***************************************************************************
       
%**********MODEL COEFFICIENTS
alpha2 = -0.00506964011908224;

a_geoid = [-5.6195e-001, -6.0794e-002, -2.0125e-001, -6.4180e-002, -3.6997e-002,...
           +1.0098e+001, +1.6436e+001, +1.4065e+001, +1.9881e+000, +6.4414e-001,...
           -4.7482e+000, -3.2290e+000, +5.0652e-001, +3.8279e-001, -2.6646e-002,...
           +1.7224e+000, -2.7970e-001, +6.8177e-001, -9.6658e-002, -1.5113e-002,...
           +2.9206e-003, -3.4621e+000, -3.8198e-001, +3.2306e-002, +6.9915e-003,...
           -2.3068e-003, -1.3548e-003, +4.7324e-006, +2.3527e+000, +1.2985e+000,...
           +2.1232e-001, +2.2571e-002, -3.7855e-003, +2.9449e-005, -1.6265e-004,...
           +1.1711e-007, +1.6732e+000, +1.9858e-001, +2.3975e-002, -9.0013e-004,...
           -2.2475e-003, -3.3095e-005, -1.2040e-005, +2.2010e-006, -1.0083e-006,...
           +8.6297e-001, +5.8231e-001, +2.0545e-002, -7.8110e-003, -1.4085e-004,...
           -8.8459e-006, +5.7256e-006, -1.5068e-006, +4.0095e-007, -2.4185e-008 ];

b_geoid =[ +0.0000e+000, +0.0000e+000, -6.5993e-002, +0.0000e+000, +6.5364e-002,...
           -5.8320e+000, +0.0000e+000, +1.6961e+000, -1.3557e+000, +1.2694e+000,...
           +0.0000e+000, -2.9310e+000, +9.4805e-001, -7.6243e-002, +4.1076e-002,...
           +0.0000e+000, -5.1808e-001, -3.4583e-001, -4.3632e-002, +2.2101e-003,...
           -1.0663e-002, +0.0000e+000, +1.0927e-001, -2.9463e-001, +1.4371e-003,...
           -1.1452e-002, -2.8156e-003, -3.5330e-004, +0.0000e+000, +4.4049e-001,...
           +5.5653e-002, -2.0396e-002, -1.7312e-003, +3.5805e-005, +7.2682e-005,...
           +2.2535e-006, +0.0000e+000, +1.9502e-002, +2.7919e-002, -8.1812e-003,...
           +4.4540e-004, +8.8663e-005, +5.5596e-005, +2.4826e-006, +1.0279e-006,...
           +0.0000e+000, +6.0529e-002, -3.5824e-002, -5.1367e-003, +3.0119e-005,...
           -2.9911e-005, +1.9844e-005, -1.2349e-006, -7.6756e-009, +5.0100e-008 ];
       
atm_mean =[278.913188097120,2.14950743203835,0.590789720042615,-23.7213627729675,0.865639902548481,...
           0.0782363927308371,1.34011243531550,0.197325826484017,-0.00374338738355388,-0.0308245520466987,...
          -3.97903875730563,0.378818877339062,0.0155068205133418,-0.000286763890860013,-0.00343631725823704,...
           1.11343491015239,-0.100201690860026,0.0443293022873576,-0.000620208658102620,0.00226756123984507,...
          -0.000231089998735478,1.58439212360353,-0.0244004753637887,0.0371853690903250,0.00653788758591591,...
           3.08096590704719e-05,3.07283541378535e-05,4.88417273193241e-06,0.658551743733080,0.129963128694739,...
           0.0249732954770270,0.00243291791905656,-0.000131687055671787,-1.42030714026868e-05,-2.40143675934510e-06,...
           4.27102336631301e-07,-2.47057610179647,-0.0626461425848206,0.00463942369757520,0.000714564693414340,...
           5.37125465908080e-06,-1.91359548073436e-05,4.64974129002234e-07,6.69943590224833e-07,-6.42730788836686e-08,...
          -0.663055710321518,0.0507518707649335,0.00304755447037716,0.000199964085915060,-7.97784721071275e-05,...
          -5.58180766637502e-06,1.09176210336729e-07,1.43122290211307e-07,3.21522865290533e-08,-2.01673232702633e-09];
       

atm_amp1 =[-1.50359389351680,-6.88539220090626,-0.137083414956074,-2.74581891273023,0.0497961507652041,...
            0.156436781744063,-1.52599020754738,0.337718889670002,0.172296955381897,-0.0178892880305058,...
            1.60590922033741,0.356611756250802,0.157434614486986,-0.0101076395762882,-0.000758982585000300,...
            1.05035822477085,0.367845073753737,0.0877647923232471,-0.00570724647955695,6.35382969348539e-05,...
            8.53672459556992e-05,1.49167770671077,0.124902997768231,0.0411918943587793,0.00192567180964808,...
           -0.000313334280728821,-3.10417840388115e-05,-2.04298751130846e-05,-2.39105162790819,-0.0374356604909260,...
           -0.00790309554079800,-0.000165806836224461,-8.01664860203305e-05,2.23873891703414e-05,-7.68595400032295e-06,...
           -1.32235241583523e-06,-0.394056610212945,0.0467100172747354,-0.00741051421600200,-0.000326420986704923,...
            0.000107521617123011,3.70022565236300e-07,-9.50788780617002e-07,-8.90319629661535e-08,-1.24917830067203e-08,...
            1.30803171224721,0.0111643754597963,-0.00686501200384918,-0.000435553920383797,5.36134667051492e-05,...
           -7.47404605496918e-07,-6.66100663451935e-07,-1.07455659116843e-07,-1.67909060562311e-09,1.01948145736122e-09];


atm_amp2 = [-0.490443092915572,-3.60816285126190,0.208940276909841,-1.29240122329504,0.0651764725821021,...
            -0.0671919714309861,0.256550989854838,0.124736654638083,-0.00744972427308259,0.00372968850089232,...
             0.647702817482626,0.0703705713545525,0.0348942999939424,-0.00190164431118696,0.00464005548254637,...
             1.27942449645463,0.0685230334953456,0.0208336318596124,-4.69816591257502e-05,0.000160270849444238,...
             0.000139111029983704,0.0171054672800839,-0.0449765925560391,0.0190704387042213,0.00238656051791216,...
             0.000121471736930913,-2.18113832173575e-05,-7.33574049451827e-06,-1.58998366067295,-0.0328191320235705,...
             0.000342351495488140,0.00118323369629953,3.33351788572796e-05,-2.15632106740390e-07,-2.61230770397944e-06,...
             1.22257094816678e-07,-0.343001802806075,0.0260803703841844,0.00364735408491998,0.000156107413452177,...
            -8.27248245073960e-05,-4.57147728228947e-07,4.17608193766559e-07,-1.22781207150326e-08,3.19793739216335e-08,...
             1.18600525888166,0.0292564420727342,-3.20220346835973e-05,0.000177648413691716,-5.65628641643874e-05,...
            -2.44248095974132e-06,2.52083535728056e-07,-1.33324406972961e-08,1.26859473536693e-09,9.54852803555643e-10];
        
atm_amp3 = [0.122999739955732,0.0587934970674077,-0.100744802072331,1.03168809037487,-0.0302876909315986,...
            0.0308833057556208,0.177564244396706,0.000691437489596723,0.0243946612367816,-1.99785872107000e-05,...
            0.905317141266527,0.0267526503120486,0.00272292153791285,-0.00196071038227297,0.000235831934108012,...
           -0.124810830757577,-0.0667527261503340,-0.00101642828686318,0.00109323585135110,-0.000524070936690564,...
            5.19126158635204e-06,0.282085532581196,-0.0315209538663076,-0.000834567768081796,-0.000296180960171001,...
           -0.000141947235079848,-1.91683337705855e-05,-5.13684067083098e-06,-0.250474634671335,-0.00998574216662525,...
            0.000291909246496190,0.000269918772274343,-3.52754756345339e-05,-8.93832029569276e-06,-1.45499055653068e-06,...
           -4.42924152090358e-08,0.149476024972113,0.0577416328026315,-0.000977327652658303,-0.000182035188699039,...
           -3.32340073421311e-05,-1.95592884361174e-06,-3.64948792784629e-08,-8.77693964452551e-08,-4.75928449957812e-09,...
           -0.0108373951815197,0.0523617648755689,0.000341742984358881,-0.000231292838925289,-9.33023315426893e-06,...
           -1.00757684418208e-06,1.44435338392427e-08,5.85966613524236e-09,-3.23345930510065e-09,-4.80742762259579e-10];  
        
atm_amp4 = [0.183761741467557,0.206274750668958,-0.0332189959428305,0.466508922534718,-0.0424519990774646,...
            0.0223976076607926,0.460472613070245,-0.102139663462668,0.0247806840790088,-0.00107632175935673,...
           -0.334713301354530,-0.00816915435026069,0.0115475754204437,0.00192297914032037,-0.000121798609127049,...
           -0.188546030616976,-0.0519763982725866,0.0110043850188291,0.00131297750780608,-0.000315734285798994,...
           -9.78124190718262e-05,0.158970118976141,0.0210951562883347,-0.00614933096375057,0.00171447046536296,...
           -0.000117657184395670,-1.83342357484406e-05,-1.02433050574879e-06,0.0734322016299409,0.0429313125044200,...
           -0.00141837011199688,-0.000986108700356642,3.12772657246579e-05,-2.21936613163248e-05,3.95716056136771e-07,...
           -1.10366905030294e-07,0.102202103204556,0.0453626398101101,-0.00178124793072085,3.17982157222843e-05,...
            1.47855216308489e-05,-7.08133745124385e-06,1.90093142636056e-07,2.45961230137684e-09,1.49196608554077e-08,...
            0.211274879890438,0.0363550493327364,0.000403700680597285,0.000231556379543750,-1.85190424472124e-06,...
           -2.45921096312441e-06,6.02652083005793e-08,1.77996395490908e-10,4.79836084813624e-09,-8.35316426777732e-11]; 
       
atm_amp5 = [-0.0205617229457509,-0.0231485419903774,-0.0524563982463164,-0.0174700606215218,-0.0139312886359497,...
            -0.00428446255857715,0.0289417616270992,0.0121393964234784,0.00534112004061830,-0.00225727799185354,...
             0.0269351330714480,0.0102973483325833,0.00168671414654829,0.00185694203356950,0.000324620114105981,...
             0.00474575549207467,-0.00351645870331169,-0.00133623697070204,0.000128365429388483,-3.64131943949197e-06,...
             6.02792147916078e-05,-0.0302546827731224,0.000344790269017128,-0.000400360446164854,-0.000257087039585194,...
             1.52198201318880e-05,-2.65899286603112e-06,4.12534640525621e-06,0.0273684511832503,-0.00206738176270874,...
             0.000347053214784983,2.16383055317555e-05,1.35612995980041e-05,-3.25583149093547e-06,-5.99829699064995e-07,...
             2.55285432867403e-07,0.0466895193483676,0.00107937714295722,0.000616209143039752,-2.26202866915416e-05,...
             9.20906913956784e-07,-7.11862904656690e-07,-4.91266804400284e-08,-7.56185887704408e-08,7.43795898693108e-10,...
             -0.0443928754454000,4.64455852785867e-05,1.20723616571776e-05,5.48431384656341e-06,2.28399258374416e-06,...
             5.98654236903030e-07,-5.47915083398008e-08,1.79899033304983e-08,-1.35035335322408e-09,-4.15545717329992e-10];
         
atm_amp6 = [-0.0336801519300105,-0.0132874927409527,-0.292432001129928,0.0444614171708192,-0.0211611485476957,...
            -0.0180129456949471,0.0227856583016144,0.0150354131576962,-0.00443115949115174,0.00143931501411051,...
            -0.0262696869951009,-0.00896664793049747,0.000925214619412443,-0.000629362342698296,1.25070804350188e-05,...
            -0.0184576429069697,-0.00603109638341317,-0.000136599431827293,0.000363722259117332,3.83843512437173e-05,...
             1.75277400470200e-05,0.0166813885578653,0.00790880799731665,0.000254400382414843,0.000333083773654716,...
            -3.39186994970671e-05,-1.02786463168231e-05,1.76209529278774e-06,0.00766282407253595,-0.00775106332791737,...
             0.00101793054497267,-3.41950494345595e-06,-6.19461610731495e-06,-1.96938257323253e-06,-4.11768566599016e-08,...
             1.12048379370978e-07,-0.00733847005172744,0.00399006673693712,-0.000414658491881626,-1.87725519454771e-05,...
             4.07109886538015e-06,3.95370149191396e-07,-9.49701889355786e-08,-6.13528930812264e-09,1.21410179053204e-08,...
            -0.00721087626443461,0.00227735849201357,9.22822228237458e-05,8.14297418448394e-06,-3.67025714801191e-06,...
            -5.35815612992430e-07,8.65023392931230e-08,-3.41053136045033e-09,-5.01692410982460e-09,1.00507807741180e-10];
        
btm_mean = [0,0,0.0778582718705853,0,0.305714476415429,0.181800239555099,0,-0.392798583933330,0.117498505504886,0.0165902608450345,...
            0,0.274711939658912,0.0440567344113229,-0.00119802230097914,0.00102872187886491,0,-0.271268923821498,0.0381865215878222,...
            0.000623455354470002,-0.00171954070502172,0.000142018254231164,0,0.316276320171668,0.0220633097853836,-0.00578445253303862,...
           -0.000447347872722017,-6.79718468611822e-05,8.56761405949400e-06,0,0.0606486211597418,0.00249090787714564,-0.00234443729512925,...
           -0.000200622667186758,2.91859909256101e-05,-6.87997556260304e-07,-1.15688947960933e-06,0,0.104690495378494,0.00268958989540007,...
            0.000218861158288698,-0.000111452042297002,1.33604955920894e-05,2.14389240746673e-06,-9.97454903590295e-08,-1.08822592714690e-07,...
            0,-0.133139857979163,-0.00414256181683629,0.000262864571525132,-1.10676203510254e-05,3.29346821092389e-06,6.30578190511812e-07,...
           -1.10016435908046e-07,9.60123516706882e-09,-1.66687335795543e-10];
btm_amp1 = [0,0,-0.592517826962367,0,-0.756852538234236,-0.120471012915935,0,-0.497346406006019,0.0429731622959948,-0.00880166352826639,...
            0,-0.105989411947159,0.0152542625167179,-0.0206498286309134,-0.00201473084085135,0,0.234648591772909,0.0343236281787142,...
           -0.0140081026707057,-0.000649484855923335,7.37211564901966e-05,0,0.0584364637325225,0.0320638413485914,-0.00519346644829413,...
           -0.000214698511820207,0.000121884398089535,-5.91161973878564e-06,0,0.215455696359060,0.00758033986119054,-0.00299824087941529,...
           -0.000199159187535078,6.29304819201578e-05,-1.90095563059853e-06,-1.22584856321955e-06,0,0.0226844323534866,-0.00756392339089110,...
            0.000397596341060141,-0.000156791560413666,1.00253645396150e-05,-8.09887568218689e-07,-1.16058022421319e-08,-4.65246476549127e-08,...
            0,0.0445968172892046,0.000476993962347011,-4.12768934787218e-05,1.38867154043243e-05,8.72310948583619e-06,-1.26389076582070e-06,...
           -9.95271907834252e-08,-1.80893311278996e-10,4.28256505130100e-10];
       
btm_amp2 = [0,0,-0.0122479856685579,0,0.0201765042638813,0.0269100047817980,0,-0.0803803898641259,0.00846882091998009,-0.00525517731125233,...
            0,-0.0967410156603188,0.0151041216752409,-0.00633576939309317,-0.00280633072886889,0,0.0902200674692880,0.0123487649722207,...
           -0.00749518093486767,-1.00592455618512e-05,-0.000166141261450248,0,0.0433210031113057,0.0130063546683708,-0.00303008055221123,...
            9.61271494625296e-05,1.35171380724485e-05,-4.13238473468550e-06,0,-0.0298780663244048,-0.0107315294583263,-0.000197362578635673,...
           -0.000128792407635824,-6.93538884869542e-06,3.52999149471513e-06,6.99074114448022e-07,0,-0.0119903137251735,-0.00435138849078847,...
           -0.000432694881401380,-3.34792953917844e-05,-5.55435846832665e-06,2.06143714378024e-07,3.27163402770885e-08,1.54402008966990e-08,...
            0,-0.0223419072728562,-0.00453735096395125,-0.000497953074790828,-1.75672806810104e-05,9.94859655807274e-08,-2.00940417160627e-07,...
            6.95790417038359e-08,-1.35303487715157e-08,-4.67938959495910e-10];
        
btm_amp3 = [0,0,-0.133741319117053,0,-0.0886107845936190,0.0133733656872964,0,0.0724603512133792,-0.00204741614867426,-0.00343028352780000,...
            0,-0.0239990865845136,0.0106952193957979,0.000697982395401740,0.000351090948998333,0,0.0705176135374181,0.00132624179604311,...
           -0.000256421605349527,0.000442951910695917,0.000111839217995342,0,-0.0432908248371864,0.00129192164198252,-0.000770509459913362,...
            7.67133670909187e-05,3.15978344324762e-05,-2.93394558209848e-06,0,0.00854343128810993,-0.00434223503010926,0.000496337142121929,...
           -1.31251600980501e-05,8.65722000607495e-06,5.86638630464148e-07,6.12594488753280e-08,0,-0.0654466050275388,-0.00273729684245883,...
            3.05028125870384e-05,-1.37401587879031e-07,5.74989499675812e-07,-3.58609499372689e-07,-1.50030298656803e-07,3.45568665613545e-08,...
            0,-0.0121405408466051,-0.00163781893784001,-5.17485167430442e-05,5.50854763054154e-06,1.40595802940089e-06,-9.02712801219109e-08,...
           -2.43178518154195e-08,1.10182025929641e-09,5.52583995382292e-10];
       
btm_amp4 = [0,0,-0.153069712160729,0,-0.0594593279256315,-0.0148049746605976,0,0.000690578822448154,0.00737272003275804,-0.00465992615639334,...
            0,0.0123552841751027,0.00351252210495320,-0.00492345590863009,0.00118736621303254,0,-0.0163130386367319,0.000466546622785591,...
           -0.000254641911212074,0.000257788323526059,3.47389019706321e-05,0,0.00958575882905318,-0.00300223048719657,-0.000824329705257649,...
           -6.94163109021607e-05,1.81867749556230e-05,-5.55979315996809e-06,0,-0.0123973970235239,-0.00400925522496225,-5.15373771874713e-06,...
           -0.000107237568280748,-1.73117542186820e-06,-1.97738760565604e-06,5.65788127065801e-08,0,0.00140627852636228,-0.00247153171376130,...
           -0.000226377221830705,-5.70208797023784e-05,3.39698594789969e-06,3.64449063896633e-07,-9.77718383074134e-08,1.82201788126035e-08,...
            0,0.0104546503530950,0.000417919412278803,-0.000168919901607586,-8.19393039432030e-06,2.32808659060194e-07,2.25925124973257e-07,...
           -2.49175091377898e-08,3.50446433655499e-09,-1.94070436075581e-10];
       
btm_amp5 = [0,0,-0.318067710502201,0,-0.0579970550611717,-0.0163148672150660,0,-0.0247995310482872,0.00184058104759062,4.78971819911175e-05,...
            0,-0.0171073512965450,0.00365130610888561,0.00145783128202004,-0.000226881889946824,0,0.00755674647985190,0.000514957580615413,...
            0.000759533035273623,-7.75857263305171e-05,4.63062748297388e-06,0,0.0226092207361425,-0.000291965510756393,4.37391342224835e-05,...
           -2.84779565213827e-05,-1.65395117249092e-05,-4.77468452947909e-07,0,-0.00623867872781672,0.000271609462267271,3.20861851847977e-05,...
           -2.16697697948683e-05,-1.44145597656715e-06,1.24486068634153e-06,-6.28426546295651e-08,0,0.00227576760609467,-0.000145778705425263,...
           -6.75672730446659e-05,3.66978040337151e-06,4.83940851973245e-07,1.77643240321121e-08,4.55141728088402e-08,-1.38950422298613e-09,...
            0,0.000640225738785700,0.000109943699032360,-4.60170764447540e-05,2.72518830703328e-06,2.60096032033877e-07,2.24172342768408e-08,...
           -7.13834824083231e-09,-5.05421486538336e-10,4.98503970427622e-10];
       
btm_amp6 = [0,0,0.0436132188092647,0,-0.0246484999501010,-0.0178634979743666,0,-0.00966695333812173,-0.00205399617586783,-0.00324902800968307,...
            0,0.00980529908298968,-0.00109467788276475,-0.000793675066756176,-0.000749623746720517,0,-0.000156427285293542,-0.000460784039751379,...
            0.000181430188608916,0.000104733632507062,-7.30321479528749e-05,0,-0.00558075912172978,0.000767300353102789,-1.88154652652516e-05,...
           -1.14727220146119e-06,1.90605770654750e-05,-6.13104399536983e-07,0,0.00531032179297282,-0.000407391318433965,-0.000107814556874618,...
           -1.53926679728800e-05,-3.55625229117244e-07,5.03753366147334e-08,2.31271308250387e-08,0,-0.00101323304224313,-0.000614150635714083,...
            1.60698789132012e-05,-1.82764827340315e-06,-1.72271878392577e-06,1.43298729059203e-07,1.30566132758761e-08,9.46777156623161e-09,...
            0,0.00199571645844959,0.000821834795051578,2.23745005788397e-05,-4.30913839487089e-06,1.72890179447845e-07,-5.33871165314281e-08,...
           -2.92212234771776e-09,-5.73455493122187e-10,6.49711127030057e-10];
       
%****CONVERT LATITUDE & LONGITUDE COORDS IN DEGREES TO RADIAN
lat=deg2rad(lat); %Latitude coord in radian
lon=deg2rad(lon); %Longitude coord in radian

%******COMPUTE DAY OF YEAR(doy) 
%***CALL THE DayofYear.m FUNCTION
 doy=DayofYear(UTCtime);

%PARAMETERS t
t = sin(lat);

%****DEGREE n and ORDER m FOR LEGENDRE POLYNOMIAL
n = 9; 
m = 9; 
i = 0;
%dfac = 20;
dfac(1) = 1;

% DETERMINE n!  (faktorielle)  moved by 1
for i = 1:(2*n + 1)
    dfac(i+1) = dfac(i)*i;
end

%***DETERMINE LEGENDRE FUNCTIONS 
%REFERENCE:Heiskanen and Moritz,Physical Geodesy, 1967, eq. 1-62
P =[10, 10];
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
horth = h - undu;

%***SURFACE TEMPERATURE ON THE GEOID
atm = 0.0;
ata1 = 0.0;
ata2 = 0.0; 
ata3 = 0.0; 
ata4 = 0.0; 
ata5 = 0.0;
ata6 = 0.0;

for i = 1:55
    atm =  atm + (atm_mean(i)*aP(i) + btm_mean(i)*bP(i));
   ata1 = ata1 + (atm_amp1(i) *aP(i) + btm_amp1(i) *bP(i));  
   ata2 = ata2 + (atm_amp2(i) * aP(i) + btm_amp2(i) * bP(i));
   ata3 = ata3 + (atm_amp3(i) * aP(i) + btm_amp3(i) * bP(i));
   ata4 = ata4 + (atm_amp4(i) * aP(i) + btm_amp4(i) * bP(i));
   ata5 = ata5 + (atm_amp5(i) * aP(i) + btm_amp5(i) * bP(i));
   ata6 = ata6 + (atm_amp6(i) * aP(i) + btm_amp6(i) * bP(i));    
end

Tm0 =  atm + ata1*cos(doy/365.25*2*pi) + ata2 * sin((doy) / 365.25 * 2 * pi) + ata3 * cos((doy) / 365.25 * 4 * pi)  + ata4 * sin((doy) / 365.25 * 4 * pi)  + ata5 * cos((doy - floor(doy)) * 24 / 24 * 2 * pi) + ata6 * sin((doy - floor(doy)) * 24 / 24 * 2 * pi);

%****HEIGHT CORRECTION[SURFACE WEIGHTED MEAN TEMPERATURE]
Tm = Tm0 + alpha2 * horth;

%******************************************END OF GTmIII.m 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%C.******SUBROUTINE TO COMPUTE SURFACE TEMPERATURE(T0)BY GPT MODEL TO BE
%        USED TO COMPUTE SURFACE MEAN TEMPERATURE(Tm0) BY BEVIS MODEL
function [Tm,horth] = Bevis_GPT(Time,lat,lon,h)
%***************************************************************************
%***************************************************************************
%DESCRIPTION: 
%            This subroutine determines Global Surface Temperature(T0)
%            based on Spherical Harmonics up to degree and order 9 and
%            computes surface weighted-mean temperature by Bevis model
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%INPUT:                                                                    *
%1.     Time = Time of observation                                         *
% Time format:[year month day hour minute seconds] OR [year month day]     *                  
%2.      lat = Latitude of Station in degrees OR [D M S]                   *
%3.      lon = Longitude of Station in degrees OR [D M S]                  *
%4.      h = Station Height(Ellipsoidal) in meters                         *

%OUTPUT:                                                                   *
%       Tm =Surface Weighted Mean Temperature                              *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================
%REFERENCE: 
% J. Böhm, R. Heinkelmann, H. Schuh, Short Note: A Global Model of Pressure 
% and Temperature for Geodetic Applications, Journal of Geodesy, 
% doi:10.1007/s00190-007-0135-3, 2007.
%==========================================================================

%**********MODEL COEFFICIENTS
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

at_mean = [ ... 
+1.6257e+001,+2.1224e+000,+9.2569e-001,-2.5974e+001,+1.4510e+000, ...
+9.2468e-002,-5.3192e-001,+2.1094e-001,-6.9210e-002,-3.4060e-002, ...
-4.6569e+000,+2.6385e-001,-3.6093e-002,+1.0198e-002,-1.8783e-003, ...
+7.4983e-001,+1.1741e-001,+3.9940e-002,+5.1348e-003,+5.9111e-003, ...
+8.6133e-006,+6.3057e-001,+1.5203e-001,+3.9702e-002,+4.6334e-003, ...
+2.4406e-004,+1.5189e-004,+1.9581e-007,+5.4414e-001,+3.5722e-001, ...
+5.2763e-002,+4.1147e-003,-2.7239e-004,-5.9957e-005,+1.6394e-006, ...
-7.3045e-007,-2.9394e+000,+5.5579e-002,+1.8852e-002,+3.4272e-003, ...
-2.3193e-005,-2.9349e-005,+3.6397e-007,+2.0490e-006,-6.4719e-008, ...
-5.2225e-001,+2.0799e-001,+1.3477e-003,+3.1613e-004,-2.2285e-004, ...
-1.8137e-005,-1.5177e-007,+6.1343e-007,+7.8566e-008,+1.0749e-009]; 

bt_mean = [ ... 
+0.0000e+000,+0.0000e+000,+1.0210e+000,+0.0000e+000,+6.0194e-001, ...
+1.2292e-001,+0.0000e+000,-4.2184e-001,+1.8230e-001,+4.2329e-002, ...
+0.0000e+000,+9.3312e-002,+9.5346e-002,-1.9724e-003,+5.8776e-003, ...
+0.0000e+000,-2.0940e-001,+3.4199e-002,-5.7672e-003,-2.1590e-003, ...
+5.6815e-004,+0.0000e+000,+2.2858e-001,+1.2283e-002,-9.3679e-003, ...
-1.4233e-003,-1.5962e-004,+4.0160e-005,+0.0000e+000,+3.6353e-002, ...
-9.4263e-004,-3.6762e-003,+5.8608e-005,-2.6391e-005,+3.2095e-006, ...
-1.1605e-006,+0.0000e+000,+1.6306e-001,+1.3293e-002,-1.1395e-003, ...
+5.1097e-005,+3.3977e-005,+7.6449e-006,-1.7602e-007,-7.6558e-008, ...
+0.0000e+000,-4.5415e-002,-1.8027e-002,+3.6561e-004,-1.1274e-004, ...
+1.3047e-005,+2.0001e-006,-1.5152e-007,-2.7807e-008,+7.7491e-009]; 

at_amp = [ ... 
-1.8654e+000,-9.0041e+000,-1.2974e-001,-3.6053e+000,+2.0284e-002, ...
+2.1872e-001,-1.3015e+000,+4.0355e-001,+2.2216e-001,-4.0605e-003, ...
+1.9623e+000,+4.2887e-001,+2.1437e-001,-1.0061e-002,-1.1368e-003, ...
-6.9235e-002,+5.6758e-001,+1.1917e-001,-7.0765e-003,+3.0017e-004, ...
+3.0601e-004,+1.6559e+000,+2.0722e-001,+6.0013e-002,+1.7023e-004, ...
-9.2424e-004,+1.1269e-005,-6.9911e-006,-2.0886e+000,-6.7879e-002, ...
-8.5922e-004,-1.6087e-003,-4.5549e-005,+3.3178e-005,-6.1715e-006, ...
-1.4446e-006,-3.7210e-001,+1.5775e-001,-1.7827e-003,-4.4396e-004, ...
+2.2844e-004,-1.1215e-005,-2.1120e-006,-9.6421e-007,-1.4170e-008, ...
+7.8720e-001,-4.4238e-002,-1.5120e-003,-9.4119e-004,+4.0645e-006, ...
-4.9253e-006,-1.8656e-006,-4.0736e-007,-4.9594e-008,+1.6134e-009]; 

bt_amp = [ ... 
+0.0000e+000,+0.0000e+000,-8.9895e-001,+0.0000e+000,-1.0790e+000, ...
-1.2699e-001,+0.0000e+000,-5.9033e-001,+3.4865e-002,-3.2614e-002, ...
+0.0000e+000,-2.4310e-002,+1.5607e-002,-2.9833e-002,-5.9048e-003, ...
+0.0000e+000,+2.8383e-001,+4.0509e-002,-1.8834e-002,-1.2654e-003, ...
-1.3794e-004,+0.0000e+000,+1.3306e-001,+3.4960e-002,-3.6799e-003, ...
-3.5626e-004,+1.4814e-004,+3.7932e-006,+0.0000e+000,+2.0801e-001, ...
+6.5640e-003,-3.4893e-003,-2.7395e-004,+7.4296e-005,-7.9927e-006, ...
-1.0277e-006,+0.0000e+000,+3.6515e-002,-7.4319e-003,-6.2873e-004, ...
-8.2461e-005,+3.1095e-005,-5.3860e-007,-1.2055e-007,-1.1517e-007, ...
+0.0000e+000,+3.1404e-002,+1.5580e-002,-1.1428e-003,+3.3529e-005, ...
+1.0387e-005,-1.9378e-006,-2.7327e-007,+7.5833e-009,-9.2323e-009]; 

%****CONVERT LATITUDE & LONGITUDE COORDS IN DEGREES TO RADIAN
lat=deg2rad(lat); %Latitude coord in radian
lon=deg2rad(lon); %Longitude coord in radian

%********COMPUTE MODIFIED JULIAN DATE(MJD)
try
   dmjd=mjuliandate(Time);
   
catch
     %********COMPUTE MODIFIED JULIAN DATE(MJD)
     %GET UTC TIME COMPONENTs
     %1.******UTCdate
     Yr = Time(:,1);%get Hour
     Mn = Time(:,2);%get Month
    Day = Time(:,3);%get Day
     
     %******GET TIME[HOUR MINUTE SECONDs]
     if size(Time,2) == 6
        H    = Time(:,4);%get Hour
        MIN  = Time(:,5);%get Minute
        SECs = Time(:,6);%get Seconds
   
        %CHANGE SECs,MIN,HOUR WHOSE SEC==60
        MIN(SECs==60) = MIN(SECs==60)+1;
        H(MIN==60) = H(MIN==60)+1;
   
   
     elseif size(Time,2) == 3
           H    = zeros(size(Time,1),1);%create Hour of zeros
           MIN  = zeros(size(Time,1),1);%get Minute of zeros
           SECs = zeros(size(Time,1),1);%get Seconds of zeros    
     end 

%****MODIFIED JULIAN DAY
[~, dmjd ,~]=utc2JulianDay_DoY(Yr,Mn,Day,H,MIN,SECs);

end

%*****COMPUTE DAY OF YEAR
%REFERENCE DAY IS 28 JANUARY.THIS IS TAKEN FROM NIELL (1996) TO BE CONSISTENT
doy = dmjd  - 44239.d0 + 1 - 28;

%%***PARAMETERS t
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
            sum = sum + (-1)^k*dfac(2*i - 2*k + 1)/dfac(k + 1)/dfac(i - k + 1)/dfac(i - j - 2*k + 1)*t^(i - j - 2*k);
        end 
        %LEGENDRE FUNCTIONs moved by 1
        P(i + 1,j + 1) = 1.d0/2^i*sqrt((1 - t^2)^(j))*sum;
    end
end 

%*****SPHERICAL HARMONICS(formula 4.19 on page 77 of "Mathematische Methoden der Geowissenschaften VO"
i = 0;
for n = 0:9
    for m = 0:n
        i = i + 1;
        Cnm(i) = P(n+1,m+1)*cos(m*lon);
        Snm(i) = P(n+1,m+1)*sin(m*lon);
    end
end

%****GEOIDAL HEIGHT
undu = 0.d0;
for i = 1:55
   undu = undu + (a_geoid(i)*Cnm(i) + b_geoid(i)*Snm(i));
end

%***ORTHOMETRIC HEIGHT
horth = h - undu;

%****SURFACE TEMPERATURE ON THE GEOID 
atm = 0.d0;
ata = 0.d0;
for i = 1:55
    atm = atm + (at_mean(i)*Cnm(i) + bt_mean(i)*Snm(i));
    ata = ata + (at_amp(i) *Cnm(i) + bt_amp(i) *Snm(i));
end 
T0 =  atm + ata*cos(doy/365.25d0*2*pi);

%****HEIGHT CORRECTION FOR TEMPERATURE
T = T0 - 0.0065d0*horth; %[in degree Celcius(°C)]

%COMPUTE SURFACE MEAN TEMPERATURE(Tm0)-BEVIS et al(1992,1994) MODEL      
Tm=(T + 273.15).*0.72+70.2;

%******************************************END OF Bevis_GPT.m 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%              --------------------------------------------------------
%D.**************SUB-ROUTINE FOR CHECKING THE EXISTENCE OF TVGG GRID FILES
%              --------------------------------------------------------
function [GTVR_grid,folder] = SearchGTVRgrid()

%*********SEARCH FOR VMF GRID FILES FROM VARIOUS DIRECTORY
%1.SEARCH GPT GRIDs from 'GPT grids' folder
[GTVR_grid,Directory] = SearchGTVRGrid('../data/Tm files');

if isequal(GTVR_grid,0)%if VMF grids are not in 'VMF files' folder, try the 'data' folder
    
   [GTVR_grid,Directory] = SearchGTVRGrid('../data');%searching the 'data' folder
   
   if isequal(GTVR_grid,0) %if VMF grids are not in the 'data' folder, try the 'current' directory
       
      [GTVR_grid,Directory] = SearchGTVRGrid(pwd);%if VMF grids are not in the 'current' directory, then issue a message
      
      if isequal(GTVR_grid,0)
          
          beep %Give a beep sound
          errmsg2{1}=sprintf('No TVGG grid file found in goGPS directory.\n');
          errmsg2{2}=sprintf('Please Provide TVGG grid file & Try Again.\n');
          warndlg(errmsg2,'VMF grid file Error','modal') 
      
      else 
          folder = Directory;
      
      end 
      
   else 
       folder = Directory;
       
   end 
   
else 
    folder = Directory;
    
end 

%            ---------------------------------------------
%************IF A DIRECTORY/FOLDER WITH VMF GRIDs IS FOUND....
%            ---------------------------------------------
 if exist('folder','var') %If Folder is not Empty
    
    %SAVE DIRECTORY 
    setappdata(0,'directory',folder) 
 end
%===================================END OF Search FOR FOLDER WITH VMF grids
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


function [GTVR_grid,Folder] = SearchGTVRGrid(directory)

%ASSIGNMENT
folder = directory;  %folder/Directory suppose to be containing VMFG FILES in goGPS

%****FIRST CHECK IF FOLDER IS EMPTY OR NOT
SizeListing=size(dir(folder)); %Size of folder content

%********CHECK IF FOLDER IS EMPTY OR NOT
if any([SizeListing(1,1)==2,SizeListing(1,1)==0])% size of 2 means folder is empty & 0 means folder doesn't exist
                  
   GTVR_grid = 0; %Flag to indicate the absence of VMF grid file in directory
   Folder = []; %File directory is empty([]) if GTVR-Tm.txt is not found
   
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
          
       %*****GET INDEXEs FOR THE VARIOUS GPT GRID FILES(i.e.Integer Val /[])
       GTVRIndex = regexp(fileList,'\w*GTVR-Tm.txt'); %matches any words starting with gpt.  
          
       %************GET THE GRID FILES
       GTVRFiles  = fileList(~cellfun(@isempty,GTVRIndex));%GTVR files(eg:'GTVR-Tm.txt')
        
    else
        GTVRFiles = [];
        
        Folder = []; %File directory is empty([]) if GTVR-Tm.txt file is not found
        
    end  %if ~isempty(fileList)
    
    %IF NO GPT FILES EXIST IN THE CURRENT FOLDER, CHECK SUB-FOLDERS/              
    %SUB-DIRECTORIES FOR GPT GRID FILEs
    %                              +                        ==========
     if isempty(fileList)
              
        checkSubfolders = 1; %Flag to search Sub-Directories/folders
           
     elseif isempty(GTVRFiles)
                   
            checkSubfolders = 1;%Flag to search Sub-Directories/folders 
                         
     else
          if ~isempty(GTVRFiles)
          
             checkSubfolders=0;%Flag not to search Sub-Directories/folders
             
             GTVR_grid = 1; %Flag to indicate the presence of VMF grid file in directory
             
             %***CONSTRUCT THE FULL PATH [File directory [if GTVR-Tm.txt file is found)]      
             Folder = cellfun(@(x) fullfile(folder,x),GTVRFiles,'UniformOutput',false);%Prepend path to files  
             
             %CONVERT CELL ARRAY TO CHARACTER
             if iscell(Folder)
                Folder=char(Folder);
             end    
             
          end
           
     end  %\\if isempty(fileList)
                
     
     if checkSubfolders==1 %if this condition is true,then either files in the current
                           %Directory were not VMF grid files or the there were no files at all.        
                           %System has to look through Sub-Folders or Sub-directories for VMF Grid Files                                                              
                                                                                                                                                                                                                                                                                  
        %***OPEN ALL SUB-FOLDERS & GET VMF GRID FILES
                
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
            
                     j=1;%Loop index for the VMF grid file types
                     
                      for i = 1:length(filelists)%****LOOP OVER THE FILES & CHECK FOR VMF GRID FILES
                                                 %Identify the needed files and extract them separately
                                                 %into different file names(H0files,H6files,H12files & H18files )
                                                 %Note that dir also lists the directories,so you have to check for them. 
                                                         
                           %***GET NAME OF FILES FROM THE LIST OF FILES FOR EACH LOOP               
                           
                           FileName = filelists{i,1};
                           filepath = filepaths{i,1};

                           %********MATCH ANY WORD STARTING WITH GPT
                           if regexpi(FileName,'\w*GTVR-Tm.txt') 

                             %When Tested & the file is aGPT file
                             GTVRfiles{j,1}=FileName; %#ok<*NASGU>
                             GTVRpath{j,1}=filepath;%18 Hour file in cells
                             
                             j=j+1;%Update index
                                                                           
                           end  %//if any([regexpi(FileName,'\<gpt\w*') regexpi(FileName,'\w*gpt')])
                           
                      end %//for i = 1:length(fileLists) 
                                                                                                               
                  end   %//if ~isempty(fileLists)
                      
            else                          
                 %WRITE ERROR MESSAGE
                 erMsg{q,1}=sprintf('The Folder %s and  It''s Sub-Folder, %s , are empty / Contain No VMF Grid files .\n',folder,subDirs{iDir});
                         
                 q=q+1; %update index
                        
            end   %if Sizesublist(1,1)~=2
                   
        end    %for iDir = find(validIndex)
        
      if ~exist('GTVRfiles','var')
         
         GTVR_grid = 0;%Flag to indicate the absence of GTVR grid file in directory
         
         Folder = []; %File directory is empty([]) if GTVR-Tm.txt file is not found
         
         return
     
      else 
          GTVR_grid = 1; %Flag to indicate the presence of GTVR grid file in directory
          
          Folder = GTVRpath; %if GTVR-Tm.txt file is found
          
          %CONVERT CELL ARRAY TO CHARACTER
          if iscell(Folder)
             Folder=char(Folder);
          end   
                    
      end %//~exist('GTVRfiles','var')                                        
           
     end  %//if checkSubfolders==1
     
end %//if exist('listFiles','var')

%=========================================END OF SearchGTVRgrid.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%E.SUBROUTINE TO COMPUTE JULIAN,MODIFIED DAY(JD,MJD) &  DAY OF YEAR(DoY)
function [JD, MJD,DoY]=utc2JulianDay_DoY(Year,Month,Day,Hour,Minute,Seconds)

%**************************************************************************
%***DESCRIPTION:
%               The function "JulianDay_DoY" Computes Julian &  Modified 
%               Julian Day as well as Day of Year (DoY) from Calender date 
 
%USAGE: [JD, MJD,DoY]=JulianDay_DoY(Year,Month,Day,Hour,Minute,Seconds)
  
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
  %3.    DoY -------->Day of Year
 
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
          
if all(UT(:))~=0
    
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

%****COMPUTE DAY OF YEAR
%Call the "DayofYear function.m"
 DoY=DayofYear(Year,Month,Day);
%******************************************END OF utc2JulianDay_DoY.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%E.1*********SUBROUTINE TO COMPUTE DAY OF YEAR(DoY)
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

%E.2*********SUBROUTINE TO CONVERT TWO DIGITs TO YEAR 4 DIGITS
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