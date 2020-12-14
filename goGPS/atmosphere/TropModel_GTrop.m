%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *
%            "TropModel_GTrop" is a subroutine that Computes Zenith & Slant*
%             Tropospheric Delays(ZHD,ZWD,ZTD,SHD,SWD,STD) & Weighted-Mean * 
%             Temperature(Tm) for specific sites near the Earth Using the  * 
%             Global Tropospheric Model(GTrop).The GTrop model coefficients 
%             are stored in the grid points on a global scale based on a...* 
%             1°x 1° spatial grid resolution from the European Centre ...  *
%             for Medium-Range Weather Forecasts(ECMWF) ERA-Interim ...    *
%             reanalysis from 1979 to 2017 to compute the ZTDs  and Tm ... * 
%             globally/worldwide. Generally, the GTrop model uses site ... *
%             position[lat lon h],Year and Day of Year(DoY) to compute ... *
%             the ZTDs and Tm. The Slant delays(STD,SHD,SWD) are computed  * 
%             using the Neill Mapping Function(NMF).                       *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%USAGE:                                                                    *
%      [STD,SHD,SWD,ZTD,ZHD,ZWD,Tm] = TropModel_GTrop(Time,Lat,Lon,h,...   *
%                                                       GTrop_coeff,SatPos)*
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%1.     Time   = Time of observation                                       *
% Time format  :[year month day hour minute seconds] OR [year month day]   *                  
%2.       lat  = Latitude of Station in degrees  OR [D M S]                *
%3.       lon  = Longitude of Station in degrees OR [D M S]                *
%4.         h  = Height of Station in meters                               *                                               
%5.GTrop_coeff = GTrop model coefficients / grid values extracted from ... *  
%                GTropCoefficient.mat 1°x 1° external grid file            * 
%NOTE1:                                                                    *          
%      The GTrop model coefficients / grid values are extracted            *
%      using the "SearchGTropgrid.m" sub-routine                           *
%
%NOTE2:                                                                    *
%      The site positions can also be entered as a single input in the ff  *
%      order [lat lon h].In such case, GTrop_coeff will take the position  * 
%      Lon and SatPos the position of h.                                   *
%  e.g.:[STD,SHD,SWD,ZTD,ZHD,ZWD,Tm] = TropModel_GTrop(Time,[Lat,Lon,h],...
%                                                        GTrop_coeff,SatPos)
%  eg2: LET POS = [Lat,Lon,h]                                              *
%[STD,SHD,SWD,ZTD,ZHD,ZWD,Tm]=TropModel_GTrop(Time,POS,GTrop_coeff,SatPos) *
%      OR ---> {When SatPos is empty([]) or STDS NOT NEEDED}               *
% [STD,SHD,SWD,ZTD,ZHD,ZWD,Tm]=TropModel_GTrop(Time,POS,GTrop_coeff,[])    *
% [STD,SHD,SWD,ZTD,ZHD,ZWD,Tm]=TropModel_GTrop(Time,POS,GTrop_coeff)       *
%                                                        
%NOTE3:                                                                    *
%      If Slant delays are not needed, the SatPos input should be empty([])*
%      i.e. SatPos = [] or no entry at all.                                *
%e.g.:[~,~,~,ZTD,ZHD,ZWD,Tm]=TropModel_GTrop(Time,Lat,Lon,h,GTrop_coeff,[])*
%                OR                                                        *
%%e.g.:[~,~,~,ZTD,ZHD,ZWD,Tm] = TropModel_GTrop(Time,Lat,Lon,h,GTrop_coeff)*

%OTHER CONSIDERATIONs:
%1.*******Satellites Elevation Angle(satEL) can also be used in place of...*
%         Satellites Position (SatPos).  format should be (nx1)            * 
%         NOTE:                                                            *
%            Elevation Angles should be decimal degrees(eg:26.981.78.102)  *  
%NOTE:                                                                     *
%1.   a)Elevation Angles should be decimal degrees(eg:26.981,78.102)       *
%     b)HeightS are also given in meters                                   *
%2.   Receiver positions are in the ff formats for SINGLE INPUT:           *
%     a)For Geographic Coordinates (n x m):[Latitude Longitude height(h)]  *
%      eg:DMS:[6 40 21 -1 33 17 187.76]  OR                                *
%          DM:[6 40.35 -1 33.2833 187.76] OR                               *
%         Deg:[6.6725 -1.5547 187.76] Decimal degrees                      *

%For Western Longitudes place negative(-) infront of the value             *
% eg:[-1 33 34] or [-001 33 34] or [ 000 -24 16] or [ -1 33.2]OR [-1.3345] *

%For Southern Latitudes place negative(-) infront of the nunmber           *
%  eg:[-8 23 14] or [-008 23 14] or [ 000 -23 26] or [ -8 23.6]OR [-8.2315]*

%     b)For ECEF(XYZ) Coordinates (3 x n) OR (n x 3)matrix : This is also  *
%       applicable to Satellite Positions(SatPos).                         *
%     I.e.: 3xn=|[X eg:[6332942.597  6332976.932  6332957.890  6332977.582|*
%               | Y    -172955.641  -172805.004  -172878.972   -172804.786|*
%               | Z]    737935.003   737647.856   737824.057   737648.519]|*
%                ----------------------------------------------------------*
%           nx3=|[ X Y Z] eg:[6332942.597  -172955.641  737935.003 |       *
%               |             6332976.932  -172805.004  737647.856 |       *
%               |             6332957.890  -172878.972  737824.057 |       *
%               |             6332977.582  -172804.786  737648.519]|       *
%                --------------------------------------------------        *                                   *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%OUTPUT:
%1.     STD : Slant Total Tropospheric Delay in meters                     *
%2.     SHD : Slant Hydrostaic Tropospheric Delay in meters                *
%3.     SWD : slant Wet Tropospheric Delay  in meters                      *
%4.     ZTD : Zenith Total Tropospheric Delay in meters                    *
%5.     ZHD : Zenith Hydrostaic Tropospheric Delay in meters               *
%6.     ZWD : Zenith Wet Tropospheric Delay  in meters                     *
%7.     Tm  : Weighted Mean Temperature                                    *

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================
%REFERENCE:                                                                *
%          Sun, Z., Zhang, B., & Yao, Y. (2019). A Global Model for        * 
%          Estimating Tropospheric Delay and Weighted Mean Temperature     *
%          Developed with Atmospheric Reanalysis Data from 1979 to 2017.   * 
%          Remote Sensing, 11(16), 1893.https://doi.org/10.3390/rs11161893 *
%
%2.        The Origional Matlab source code of GTrop model is available at * 
%          https://github.com/sun1753814280/GTrop.                         *
%---------------------------------------------------------------------------
%3.        Niell, A.E.(1996).Global mapping functions for the atmosphere   * 
%          delay at radio wavelengths. J. Geophys. Res.,101, pp.3227-3246  *
%==========================================================================+                                 *
%Codes Modified by: 
%                  Samuel Osah,Msc Geomatic Engineering,(PhD Student,KNUST)*    
%                  Email: osahsamuel@yahoo.ca                              *
%                  Tel:+233 (0)246137410/+233 (0)509438484                 *
%==========================================================================
%**************************************************************************+
%***************************************************************************
function [STD,SHD,SWD,ZTD,ZHD,ZWD,Tm] = TropModel_GTrop(Time,Lat,Lon,h,...
                                                        GTrop_coeff,SatPos)

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
           
            if (any(nargin==[6,5])) %IF 6 or 5 inputs are provided
                
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
               
               if nargin == 5   %IF 5 INPUTS ARE PROVIDED
              
                  SatPos = [];%Assigning empty matrix ([]) to SatPos
                     
               end
               
            elseif  (any(nargin==[4,3,2])) %IF THREE INPUTS ARE PROVIDED
                   
                     if nargin == 4 %IF 4 INPUTS ARE PROVIDED
                   
                        %CHECK IF h POSITION IS NOT GTrop GRID VALUES(GTrop_coeff)
                        %OR GTrop_coeff & SatPos ASSIGNED TO Lon & h POSITIONS RESPECTIVELY
                        %NOTE:
                        %     THIS HAPPENS WHEN USER DECIDE TO PROVIDE LAT,LON & h
                        %     COORD AS A SINGLE INPUT,I.E.[Lat Lon h].IN THAT CASE,
                        %     GTrop_coeff WILL TAKE THE POSITION OF Lon & SatPos THAT OF h
                        
                        [nrow,ncol]   = size(h);%GET SIZE OF h
                        [nrow1,ncol1] = size(Lon);%GET SIZE OF Lon
                        
                        %***ROUND h,Lon TO THE NEAREST THOUSAND             
                        t1 = roundn(h,3);
                        t2 = roundn(Lon,3);
                        
                        if all([all([ndims(Lon) == 3,size(Lon,1) == 181,size(Lon,2) == 361,size(Lon,3) == 37]),...
                                any([any([isempty(h),isempty(t1)]),any([(nrow > 1 & ncol == 1),(nrow == 1 & ncol > 1),all(t1(:,:))==0]),...
                                                                   any([(ncol == 3 & (nrow > ncol | nrow < ncol | nrow == ncol)),...
                                                                        (nrow == 3 & (ncol > nrow | ncol < nrow | ncol == nrow))])])])
                           
                           %ASSIGNMENT
                           GTrop_coeff = Lon; %ASSUMING Lon POSITION TO BE GTrop_coeff
                           SatPos      = h;%ASSUMING h POSITION TO BE SatPos 
                           Rpos        = Lat; %ASSIGNING Lat TO Rpos(RECEIVER/SITE POSITON)
                             
                        elseif all([all([ndims(h) == 3,size(h,1) == 181,size(h,2) == 361,size(h,3) == 37]),...
                                    any([any([isempty(Lon),isempty(t2)]),any([(nrow1 > 1 & ncol1 == 1),(nrow1 == 1 & ncol1 > 1),all(t2(:,:))==0]),...
                                                                         any([(ncol1 == 3 & (nrow1 > ncol1 | nrow1 < ncol1 | nrow1 == ncol1)),...
                                                                              (nrow1 == 3 & (ncol1 > nrow1 | ncol1 < nrow1 | ncol1 == nrow1))])])])
                                
                                %ASSIGNMENT
                                GTrop_coeff = h; %ASSUMING h POSITION TO BE GTrop_coeff
                                SatPos      = Lon;%ASSUMING Lon POSITION TO BE SatPos 
                                Rpos        = Lat; %ASSIGNING Lat TO Rpos(RECEIVER/SITE POSITON)                        
                                                       
                        elseif all([all([ndims(h) == 3,size(h,1) == 181,size(h,2) == 361,size(h,3) == 37]),...
                                    ~any([any([isempty(Lon),isempty(t2)]),any([(nrow1 > 1 & ncol1 == 1),(nrow1 == 1 & ncol1 > 1),all(t2(:,:))==0]),...
                                                                          any([(ncol1 == 3 & (nrow1 > ncol1 | nrow1 < ncol1 | nrow1 == ncol1)),...
                                                                               (nrow1 == 3 & (ncol1 > nrow1 | ncol1 < nrow1 | ncol1 == nrow1))])])])                              
                                                                              
                                %ASSIGNMENT
                                GTrop_coeff = h; %ASSUMING h POSITION TO BE GTrop_coeff
                                SatPos      = [];%ASSUMING Lon POSITION TO BE SatPos 
                                
                                h = zeros(size(Lon,1),1);%Assigning zeros to  heights(h)
                                       
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
                                
                                
                        else
                             GTrop_coeff = []; %Assigning empty matrix ([]) to GTrop_coeff
                             SatPos      = [];%Assigning empty matrix ([]) to SatPos
                             
                             if all([isempty(Lon),isempty(h)])
                                
                                Rpos = Lat; %ASSIGNING Lat TO Rpos(RECEIVER/SITE POSITON) 
                               
                             elseif all([~isempty(Lon),isempty(h)])
                                
                                   if any(ncol1==[3,2,1])
                                     
                                      h = zeros(size(Lon,1),1);%Assigning zeros to  heights(h)
                                       
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
                                      
                                   else
                                       Rpos = Lat; %ASSIGNING Lat TO Rpos(RECEIVER/SITE POSITON)
                                   
                                   end %//if any(ncol1==[3,2,1])
                                      
                             else  
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
                                 
                             end %//if all([isempty(Lon),isempty(h)])
                          
                        end %//if all([all([ndims(Lon) == 3,size(Lon,1) == 181,size(Lon,2) == 361,size(Lon,3) == 37]),...
                            %          any([isempty(h),any([(nrow > 1 & ncol == 1),(nrow == 1 & ncol > 1)]),...
                            %                          any([(ncol == 3 & (nrow > ncol | nrow < ncol | nrow == ncol)),...
                            %                               (nrow == 3 & (ncol > nrow | ncol < nrow | ncol == nrow))])])
                            
                                
                     elseif nargin == 3 %IF 3 INPUTS ARE PROVIDED
                            
                            [nrow1,ncol1] = size(Lon);%GET SIZE OF Lon
                            
                            %ROUND Lon TO THE NEAREST THOUSAND
                            t2 = roundn(Lon,3);
                            
                            if all([ndims(Lon) == 3,size(Lon,1) == 181,size(Lon,2) == 361,size(Lon,3) == 37])
                          
                               Rpos        = Lat; %ASSIGNING Lat TO Rpos(RECEIVER/SITE POSITON) 
                               GTrop_coeff = Lon;%ASSUMING Lon POSITION TO BE GTrop_coeff
                               SatPos      = [];%Assigning empty matrix ([]) to SatPos 
                               
                            elseif any([any([(nrow1 > 1 & ncol1 == 1),(nrow1 == 1 & ncol1 > 1),all(t2(:,:))==0]),...
                                        any([(ncol1 == 3 & (nrow1 > ncol1 | nrow1 < ncol1 | nrow1 == ncol1)),...
                                             (nrow1 == 3 & (ncol1 > nrow1 | ncol1 < nrow1 | ncol1 == nrow1))])])
                                   
                                   Rpos        = Lat; %ASSIGNING Lat TO Rpos(RECEIVER/SITE POSITON)      
                                   GTrop_coeff = [];%Assigning empty matrix ([]) to GTrop_coeff      
                                   SatPos      = Lon;%ASSUMING Lon POSITION TO BE SatPos  
                                   
                            else
                                
                                 GTrop_coeff = []; %Assigning empty matrix ([]) to GTrop_coeff
                                 SatPos      = [];%Assigning empty matrix ([]) to SatPos 
                                 
                                 if ~isempty(Lon)
                                     
                                    if any(ncol1==[3,2,1])
                                     
                                       h = zeros(size(Lon,1),1);%Assigning zeros to  heights(h)
                                       
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
          
                                    else 
                                        Rpos    = Lat; %ASSIGNING Lat TO Rpos(RECEIVER/SITE POSITON) 
                                        
                                    end %//if any([ncol1==[3,2,1]])
                                    
                                 else
                                     Rpos    = Lat; %ASSIGNING Lat TO Rpos(RECEIVER/SITE POSITON) 
                                     
                                 end %//if ~isempty(Lon)
                                 
                            end %//if all([ndims(Lon) == 3,size(Lon,1) == 181,size(Lon,1) == 361,size(Lon,3) == 37])
                                       
                     elseif nargin == 2
                             
                            Rpos        = Lat; %ASSIGNING Lat TO Rpos(RECEIVER/SITE POSITON) 
                            GTrop_coeff = []; %Assigning empty matrix ([]) to GTrop_coeff
                            SatPos      = [];%Assigning empty matrix ([]) to SatPos 
                             
                     end   
                      
            end %if (any(nargin==[6,6]))
            
        end %//if (any(nargin==[6,5,4,3,2]))
        
        
        %CHECK IF INPUTS(SATELLITE POSITION IS NOT EMPTY([])
        %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          %1.SATELLITE POSITION [SatPos]
           if (~isempty(SatPos)| ~isnan(SatPos))
              
             %***ROUND TO THE NEAREST THOUSAND             
             t2=roundn(SatPos,3);
                  
             if  size(SatPos,2)==1 %if the number of columns in SatPos is 1 
                 
               if all(t2(:,1))==0 
                  Satpos = [];
                  satEL_deg = SatPos;%Assigning SatPos to satEL(Satellite Elevation) 
                  satEL_rad = satEL_deg.* pi / 180;%Elevation Angles converted to radian
               end   %if all(t2(:,1))==0 
               
             else
                  if all(t2(:,:))==0  
                     Satpos = []; 
                     satEL_deg = SatPos;%Assigning SatPos to satEL(Satellite Elevation) 
                     satEL_rad = satEL_deg.* pi / 180;%Elevation Angles converted to radian
                     
                  else  
                      Satpos=SatPos;%Satellite position in XYZ
                      
                  end %if all(t2(:,:))==0 
                
             end   %\\if size(SatPos,2)==1
             
           else
               Satpos = [];
               satpos_empty = 1; %flag to indicate that satpos is empty([])
               setappdata(0,'satpos_empty',satpos_empty) 
               
           end %\\if (~isempty(SatPos)| ~isnan(SatPos))
            
         %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
         %******************************************************************
         %         EXTRACT USER INPUT RECEIVER Coordinates,.....
         %*****************************************************************   
        if exist('Rpos','var') 
          
           %CHECK RECEIVER POSITION TYPE(XYZ/LAT LONG h) AS A SINGLE INPUT
           [Rrow,Rcol] = size(Rpos);%Get number of rows & columns of Rpos
          
          %***ROUND TO THE NEAREST THOUSAND
          t3=roundn(Rpos,3);
          
          if (Rcol>7 & Rrow==3)
             X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates 
             Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates 
             Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates
             
          elseif (Rcol==7 &(Rrow < 3 | Rrow > 3))
                lat=Rpos(:,1:3);%Assigning columns 1-3 to Latitude in DMS 
                lon=Rpos(:,4:6);%Assigning columns 4-6 to Longitude in DMS 
                h=Rpos(:,end);%Assigning end (7th) column to heights 
                
          elseif (Rcol==7 & Rrow==3)
                if all(t3(:,1:6))==0 %if columns 1-6 rounded to zeros(0's)
                                     %then inputs are latitude & longitude
                  lat=Rpos(:,1:3);%Assigning columns 1-3 to Latitude in DMS 
                  lon=Rpos(:,4:6);%Assigning columns 4-6 to Longitude in DMS 
                  h=Rpos(:,end);%Assigning end (7th) column to  heights 
                  
                elseif (all(t3(:,1:end))~=0 | all(t3(:,1:6))~=0 | ...
                                                        any(t3(:,1:6))~=0)
                      X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates 
                      Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates 
                      Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates 
                      
                end %if all(t3(:,1:6))==0
                                     
          elseif (Rcol==6 && (Rrow < 3 | Rrow > 3))
                lat=Rpos(:,1:3);%Assigning columns 1-3 to Latitude in DMS 
                lon=Rpos(:,4:6);%Assigning columns 4-6 to Longitude in DMS 
                 h=zeros(size(Rpos,1),1);%Assigning zeros to  heights
                          
          elseif (Rcol==6 & Rrow==3) 
                if all(t3(:,1:6))==0 %if columns 1-6 rounded to zeros(0's)
                                     %then inputs are latitude & longitude
                  lat=Rpos(:,1:3);%Assigning columns 1-3 to Latitude in DMS 
                  lon=Rpos(:,4:6);%Assigning columns 4-6 to Longitude in DMS 
                   h=zeros(size(Rpos,1),1);%Assigning zeros to  heights
                  
                elseif (all(t3(:,1:end))~=0 | any(t3(:,1:end))~=0)                                 
                      X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates 
                      Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates 
                      Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates        
                end %if all(t3(:,1:6))==0
                
          elseif (Rcol==5 && (Rrow < 3 | Rrow > 3))
                lat=Rpos(:,1:2);%Assigning columns 1-2 to Latitude in DM 
                lon=Rpos(:,3:4);%Assigning columns 3-4 to Longitude in DM 
                 h=Rpos(:,end);%Assigning end (5th) column to  heights
                 
          elseif  (Rcol==5 & Rrow==3) 
                if all(t3(:,1:4))==0 %if columns 1-4 rounded to zeros(0's)
                                     %then inputs are latitude & longitude
                  lat=Rpos(:,1:2);%Assigning columns 1-2 to Latitude in DM 
                  lon=Rpos(:,3:4);%Assigning columns 3-4 to Longitude in DM 
                   h=Rpos(:,end);%Assigning end (5th) column to  heights
                  
                elseif (all(t3(:,1:end))~=0 | all(t3(:,1:4))~=0 | ...
                                                        any(t3(:,1:4))~=0)                                 
                      X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates 
                      Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates 
                      Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates         
                end %if all(t3(:,1:4))==0
                
          elseif (Rcol==4 & (Rrow < 3 | Rrow > 3))
                lat=Rpos(:,1:2);%Assigning columns 1-2 to Latitude in DM 
                lon=Rpos(:,3:4);%Assigning columns 3-4 to Longitude in DM 
                 h=zeros(size(Rpos,1),1);%Assigning zeros to  heights
                 
          elseif (Rcol==4 & Rrow==3) 
                if all(t3(:,1:4))==0 %if columns 1-4 rounded to zeros(0's)
                                     %then inputs are latitude & longitude
                  lat=Rpos(:,1:2);%Assigning columns 1-2 to Latitude in DM 
                  lon=Rpos(:,3:4);%Assigning columns 3-4 to Longitude in DM 
                   h=zeros(size(Rpos,1),1);%Assigning zeros to  heights
                  
                elseif (all(t3(:,1:end))~=0 | all(t3(:,1:4))~=0 | ...
                                                        any(t3(:,1:4))~=0)                                 
                      X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates 
                      Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates 
                      Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates         
                end %if all(t3(:,1:4))==0
                
          elseif  Rcol==3 
              
                if all(t3(:,1:2))==0 %if columns 1-2 rounded to zeros(0's)
                                     %then inputs are latitude & longitude
                  lat=Rpos(:,1);%Assigning columns 1 to Latitude in dec.D 
                  lon=Rpos(:,2);%Assigning columns 2 to Longitude in dec.D 
                   h=Rpos(:,end);%Assigning end (3rd) column to  heights
                  
                elseif (all(t3(:,1:end))~=0 | all(t3(:,1:2))~=0 | ...
                                                        any(t3(:,1:2))~=0)
                      if (Rrow > 3 | Rrow < 3 )                            
                        X=Rpos(:,1);%Assigning 1st column to ECEF(X) Coordinates 
                        Y=Rpos(:,2);%Assigning 2nd column to ECEF(Y) Coordinates 
                        Z=Rpos(:,3);%Assigning 3rd column to row ECEF(Z) Coordinates
                        
                      elseif (Rrow == 3 )
                            if (all(Rpos(:,1)>0) & all(Rpos(:,2)<0) &...
                                                          all(Rpos(:,3)>0))
                              X=Rpos(:,1);%Assigning 1st column to ECEF(X) Coordinates 
                              Y=Rpos(:,2);%Assigning 2nd column to ECEF(Y) Coordinates 
                              Z=Rpos(:,3);%Assigning 3rd column to row ECEF(Z) Coordinates 
                              
                            elseif (all(Rpos(1,:)>0) & all(Rpos(2,:)<0) &...
                                                          all(Rpos(3,:)>0))
                                  X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates 
                                  Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates 
                                  Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates 
                                  
                            elseif  (all(Rpos(1,:)>0) & any(Rpos(2,:)<0) &...
                                                          all(Rpos(3,:)>0))
                                  X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates 
                                  Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates 
                                  Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates 
                     
                            elseif (all(Rpos(:,:)>0))
                                %IDENTIFYING DATA ARRANGEMENT IN Rpos if
                                %                    ALL XYZ INPUTs ARE +VE
                                mrow=mean(Rpos,2);%Find mean of Rpos along the rows
                                mcol=mean(Rpos,1);%Find mean of Rpos along the columns
                                %ROUND EACH MEAN VALUEs to 2 sifnificant
                                %                                   Figures
                                rmrow_2=round(mrow,2,'significant');
                                rmcol_2=round(mcol,2,'significant');
                                
                                if (strcmp(num2str(rmcol_2(1)),num2str(rmcol_2(2)))...
                                        | strcmp(num2str(rmcol_2(2)),num2str(rmcol_2(3)))...
                                        | strcmp(num2str(rmcol_2(1)),num2str(rmcol_2(3))))
                                    
                                   X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates 
                                   Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates 
                                   Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates
                                   
                                elseif (all(abs(mean(rmcol_2)- rmcol_2)==0) |...
                                               (rmcol_2(1)==rmcol_2(2) | ...
                                                rmcol_2(2)==rmcol_2(3) |... 
                                                rmcol_2(1)==rmcol_2(3)))
                                            
                                      X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates 
                                      Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates 
                                      Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates
                                      
                                elseif (~strcmp(num2str(rmrow_2(1)),num2str(rmrow_2(2)))...
                                        | ~strcmp(num2str(rmrow_2(2)),num2str(rmrow_2(3)))...
                                        | ~strcmp(num2str(rmrow_2(1)),num2str(rmrow_2(3))))
                                    
                                      X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates 
                                      Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates 
                                      Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates                                                                     

                                else                                 
                                    X=Rpos(:,1);%Assigning 1st column to ECEF(X) Coordinates 
                                    Y=Rpos(:,2);%Assigning 2nd column to ECEF(Y) Coordinates 
                                    Z=Rpos(:,3);%Assigning 3rd column to row ECEF(Z) Coordinates
                                    
                                end %if (strcmp(num2str(rmcol_2(1)),...
                 
                            end %if (all(Rpos(:,1)>0) || all(Rpos(:,2)<0)       

                      end %if (Rrow > 3 || Rrow < 3 )
                             
                end %if all(t3(:,1:4))==0
                
          elseif  (Rcol==2 &(Rrow < 3 | Rrow > 3))
                lat=Rpos(:,1);%Assigning columns 1 to Latitude in dec.D 
                lon=Rpos(:,2);%Assigning columns 2 to Longitude in dec.D 
                 h=zeros(size(Rpos,1),1);%Assigning zeros to  heights
                 
          elseif   (Rcol==2 & Rrow ==3) 
                if all(t3(:,1:2))==0 %if columns 1-2 rounded to zeros(0's)
                                     %then inputs are latitude & longitude
                  lat=Rpos(:,1);%Assigning columns 1 to Latitude in dec.D 
                  lon=Rpos(:,2);%Assigning columns 2 to Longitude in dec.D 
                   h=zeros(size(Rpos,1),1);%Assigning zeros to  heights     
                else
                    X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates 
                    Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates 
                    Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates      
                end %if all(t3(:,1:2))==0
                
          elseif ((Rcol==1 & Rrow ==3) & all(t3(:,1))~=0)
                X=Rpos(1,:)';%Assigning 1st row to ECEF(X) Coordinates 
                Y=Rpos(2,:)';%Assigning 2nd row to ECEF(Y) Coordinates 
                Z=Rpos(3,:)';%Assigning 3rd to row ECEF(Z) Coordinates
                
          end %if (Rcol>7 && Rrow==3)
          
        %*******CONVERT Lat,Long,h IF ANY TO XYZ   
        if exist('lat','var') || exist('lon','var')%Check if lat/lon exist
          
          %Call the 'geo2xyz.m' function 
          [X,Y,Z] = geo2xyz(lat,lon,h); 
          
          %*****CONVERT COORD IN DEGREEs TO RADIAN (deg2rad)
          degrad = pi/180.0; %deg2rad conversion factor
          
          if size(lat,2)== 3
             latD = dms2degrees(lat);%Convert Latitude in DMS to degrees
             lonD = dms2degrees(lon);%Convert Longitude in DMS to degrees
             
          elseif size(lat,2)== 2
                 latD = dm2degrees(lat);%Convert Latitude in DM to degrees
                 lonD = dm2degrees(lon);%Convert Longitude in DM to degrees 
                 
          else
              latD = lat;%Latitude in Decimal degrees
              lonD = lon;%Longitude in Decimal degrees 
                 
          end
          
          latR = latD.*degrad; %latitude in radian
          lonR = lonD.*degrad; %longitude in radian
          
          %SAVE COORDs for use
          setappdata(0,'latD',latD)
          setappdata(0,'lonD',lonD)
          setappdata(0,'latR',latR)
          setappdata(0,'lonR',lonR)
          setappdata(0,'ht',h)
                    
          %**CONVERT USER POSITION(XYZ) TO lat,long,h
        elseif (exist('X','var')|| exist('Y','var') || exist('Z','var'))
               
              %Call the 'xyz2LLH.m' function
              [latR,lonR,h,latD,lonD] = xyz2LLH(X,Y,Z);  %#ok<*ASGLU>
              
              %SAVE COORDs for use
              setappdata(0,'latD',latD)
              setappdata(0,'lonD',lonD)
              setappdata(0,'latR',latR)
              setappdata(0,'lonR',lonR)
              setappdata(0,'ht',h)
                  
        end %exist('lat','var') || exist('lon','var')
        
          %*******COMPUTE SATELLITE ELEVATION ANGLE IF NOT PROVIDED                   
        if ~exist('satEL_deg','var')%If Satellite Elevation is not given,then 
                                    %Compute Satellite Elevation with the given  
                                    %Receiver & Satellite positions
           
           if ~isempty(Satpos)                                  
             %Call the "satAzEl.m" function 
             [satEL_deg,satEL_rad]=satAzEl([X Y Z],Satpos);%#ok<*ASGLU> %compute satellite elevation angle             
           end
         
        end  %if ~exist('satEL','var'
        
        else
            latD = lat;%Latitude in Decimal degrees
            lonD = lon;%Longitude in Decimal degrees
            
        end  %if exist('Rpos','var') 
        
    otherwise   
             %ISSUE ERROR MESSAGE for 1 INPUT
              beep%Give a beep sound 
              errmsg{1}=sprintf('Insuficient Data Input / Wrong Data Input format .');
              errmsg{2}='';
              errmsg{3}='Please Check inputs / Data format & Try Again.';
              errordlg(errmsg,'Coordinate(s) Input Error','modal')  
              
              %Return empty ([]) outputs
              STD=[]; SHD=[]; SWD=[]; ZTD=[]; ZHD=[]; ZWD=[]; Tm=[];
              
              return                          
        
end  %switch nargin
%===============================END OF INPUT REFORMATTING
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%************TROPOHERIC DELAY MODELING USING UNB3m MODEL
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-          

%1.***************COMPUTE ZENITH TROPOSPHERIC DELAYs
%Call the "ZenithTropDelay_UNB3m.m" function
[ZHD,ZWD,ZTD,Tm] = getZTD_GTrop(UTCtime,latD,lonD,h,GTrop_coeff);

%(2)*********COMPUTE SLANT TROPOSPHERIC DELAY
if any([~exist('satpos_empty','var'),exist('satEL_deg','var'),~isempty(Satpos)])
    
   %(2.1)****COMPUTE MAPPING FUNCTIONs(MF)
   %Call the "NMF.m" function    
   [MFh,MFw] = NMF(UTCtime,latD,h,satEL_deg);                      
                            
   %Call the "SlantTropDelay.m" subroutine
   [STD,SHD,SWD] = SlantTropDelay(ZHD,ZWD,MFh,MFw);
   
else
     STD=zeros(size(ZHD));
     SHD=zeros(size(ZHD));
     SWD=zeros(size(ZHD));     
end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+


%**********VARIOUS SUB-ROUTINES TO COMPUTE TROPOSPHERIC DELAYS*************
%          ---------------------------------------------------
%
%A.SUBROUTINE TO COMPUTE SLANT TROPOSPHERIC DELAY (SATELLITE IN VIEW)
function [STD,SHD,SWD] =SlantTropDelay(ZHD,ZWD,MFh,MFw)

%**************************************************************************    
%**************************************************************************  
%***DESCRIPTION:
%              "SlantTropDelay" Computes Slant Tropospheric Delays Using   * 
%              ZHD & ZWD computed from Global Tropospheric Model(GTrop) &  *
%              MFh,MFw computed from the Niell's mapping functions(NMF)    *
%              Latitude and Day of Year (DOY) dependent Tropospheric model.*
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%USAGE:                                                                    *
%      [STD,SHD,SWD] =SlantTropDelay(ZHD,ZWD,MFh,MFw)                      *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%******INPUT:                                                              *                                                              *
%           ZHD : Zenith Hydrostatic Delay computed from the GTrop model   *
%           ZWD : Zenith Wet Delay computed from the GTrop model           *   
%           MFh : Hydrostatic Mapping Function computed from NMF model     *   
%           MFw : Wet Mapping Function computed from NMF model             * 
%****OUTPUT:                                                               *
%           SHD : Slant Hydrostatic Delay                                  *
%           SWD : Slant Wet Delay                                          *
%           STD : Slant Total Delay                                        *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+


%*********COMPUTE SLANT TROPOSPHERIC DELAY (SATELLITE IN VIEW)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%(1)*****INITIALIZING OUTPUT VARIABLEs

if any([iscell(MFh),iscell(MFw)])%if either of MFh,MFw is a cell array
    
  [nrow,ncol,dim] = size(MFh);%Get size of cell array
  SHD = cell(nrow,ncol,dim);
  [STD,SWD] = deal(SHD);%copy the contents of SHD to STD,SWD
  
  %(2)*****LOOP OVER TIME SETs,USER POSITIONs & SATELLITE ELEVATIONs
  
  for k = 1:dim %Loop over different sets of time
      
      %CONVERT CELL ARRAY TO NUMERIC ARRAY
      MFH = cell2mat(MFh(:,:,k) );
      
      MFW = cell2mat(MFw(:,:,k) );
                 
      for i=1:ncol %Loop over the Number of Receiver positions
    
          for j=1:nrow %Loop over the Number of Satellite Elevations
    
              %COMPUTE Slant Hydrostatic Delay(SHD)
              SHD{j,i,k}=ZHD(i,k).*MFH(j,i);

             %COMPUTE Slant Wet Delay(SWD)
             SWD{j,i,k}=ZWD(i,k).*MFW(j,i);

             %****COMPUTE TOTAL SLANT TROPOSPHERIC DELAY(STD)
             STD{j,i,k}=(ZHD(i,k).*MFH(j,i))+ (ZWD(i,k).*MFW(j,i));            
                  
          end %//for j=1:nrow
          
      end   %//for i=1:ncol
  
  end %//k = 1:dim
  
else %if Neither of MFh,MFw is a cell array
    
    SHD = zeros(size(MFh));%Assign zeros of nxm to SHD
    
    [STD,SWD] = deal(SHD);%copy the contents of SHD to  STD,SWD
    
    %(1)*****************LOOP OVER USER POSITIONs & SATELLITE ELEVATIONs
    %(1.1)FIND NUMBER OF ROWS & COLUMNs IN MFh
    [nrow,ncol]=size(MFh);

    for i=1:ncol %Loop over the Number of Receiver positions
    
        for j=1:nrow %Loop over the Number of Satellite Elevations
    
            %COMPUTE Slant Hydrostatic Delay(SHD)
            SHD(j,i)=ZHD(i,1).*MFh(j,i);

            %COMPUTE Slant Wet Delay(ZWD)
            SWD(j,i)=ZWD(i,1).*MFw(j,i);

            %****COMPUTE TOTAL SLANT TROPOSPHERIC DELAY(STD)
            STD(j,i)=SHD(j,i)+ SWD(j,i);            
                  
        end  %for j=1:nrow
          
    end     %for i=1:ncol

end %//if any([iscell(MFh),iscell(MFw)])
%***********************END OF SlantTropDelay.m ***********************    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%B.******SUBROUTINE TO COMPUTE ZENITH TROPOSPHERIC DELAYs
function [ZHD,ZWD,ZTD,Tm] = getZTD_GTrop(UTCtime,lat,lon,h,GTrop_Coeff)

%********COMPUTE ZTDs USING GPT2w MODEL
%GET SIZE OF USER INPUT TIME & POSITION
nrow_time=size(UTCtime,1);
nrow_pos=size(lat,1);

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
      ZHD = zeros(nrow_time,1);
      
      [ZWD,ZTD,Tm] = deal(ZHD);%create copy of ZHD in ZWD,ZTD,Tm
   
      for i=1:nrow_time
       
           %Call the "ZenithTropDelay.m" Function
           [ZHD(i,1),ZWD(i,1),ZTD(i,1),Tm(i,1)] = ZenithTropDelay(UTCtime(i,:),lat(i),lon(i),h(i),GTrop_Coeff); 
      end
      
   else 
       %*****INITIALIZE OUTPUTs 
       ZHD=zeros(nrow_pos,nrow_time);
      
      [ZWD,ZTD,Tm]=deal(ZHD);%create copy of ZHD in ZWD,ZTD,Tm
    
     for i=1:nrow_time %LOOP OVER TIME
         
        for j=1:nrow_pos %LOOP OVER POSITIONS
            
            %Call the "ZenithTropDelay.m" Function
           [ZHD(j,i),ZWD(j,i),ZTD(j,i),Tm(j,i)] = ZenithTropDelay(UTCtime(i,:),lat(j),lon(j),h(j),GTrop_Coeff); 
    
        end 
     end
    
   end %//if isequal(Identical,1)
    
else
    %*****INITIALIZE OUTPUTs 
    ZHD=zeros(nrow_pos,nrow_time);
    
    [ZWD,ZTD,Tm]=deal(ZHD);%create copy of ZHD in ZWD,ZTD,Tm
    
    for i=1:nrow_time %LOOP OVER TIME
        
        for j=1:nrow_pos %LOOP OVER POSITIONS
            
            %Call the "ZenithTropDelay.m" Function
            [ZHD(j,i),ZWD(j,i),ZTD(j,i),Tm(j,i)] = ZenithTropDelay(UTCtime(i,:),lat(j),lon(j),h(j),GTrop_Coeff);     
        end
    end
end
%***********************END OF getZTD_GTrop.m ***********************    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%(B.1)****SUBROUTINE TO COMPUTE ZENITH TROPOSPHERIC DELAY USING GTrop MODEL
function [ZHD,ZWD,ZTD,Tm] = ZenithTropDelay(UTCtime,lat,lon,h,GTrop_Coeff)

%***************************************************************************
%***************************************************************************
%DESCRIPTION:
%            "ZenithTropDelay" Computes Zenith Tropospheric Delays(ZTD,... * 
%            ZHD,ZWD) and Weighted Mean Temperature(Tm) for specific sites *
%            near the Earth Using the Global Tropospheric Model(GTrop).The *
%            GTrop model is based on a  1°x 1° spatial grid resolution from*
%            the European Centre for Medium-Range Weather Forecasts(ECMWF) *
%            ERA-Interim reanalysis from 1979 to 2017 to compute the ZTDs  * 
%            and Tm globally/worldwide. Generally, the GTrop model uses... *
%            site position[lat lon h],Year and Day of Year(DoY) to ...     *
%            compute the ZTDs and Tm.                                      *        
%%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%USAGE:                                                                    *                                                           *
%      [ZTD,ZHD,ZWD,Tm] = ZenithTropDelay(UTCtime,lat,lon,h,GTrop_Coeff)   * 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%****INPUTs:                                                               *
%1.       UTCtime:.........UTC time in [Year,Month,Day,Hour,Minute,Seconds]*                         
%2.           lat:.........Station geodetic latitude in [degrees]          *
%3.           lon:.........station ellipsoidal longitude in [degrees]      *
%4.           h  :.........Station ellipsoidal height in [km]              *
%5.  GTrop_Coeff :.........GTrop grid values or coefficients               *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%OUTPUTs:                                                                  *
%1.     ZTD : Zenith Total Tropospheric Delay in meters                    *
%2.     ZHD : Zenith Hydrostaic Tropospheric Delay in meters               *
%3.     ZWD : Zenith Wet Tropospheric Delay  in meters                     *
%4.     Tm  : Weighted Mean Temperature(Tm)  in Kelvin(K)                  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%REFERENCE:                                                                *
%          Sun, Z., Zhang, B., & Yao, Y. (2019). A Global Model for        * 
%          Estimating Tropospheric Delay and Weighted Mean Temperature     *
%          Developed with Atmospheric Reanalysis Data from 1979 to 2017.   * 
%          Remote Sensing, 11(16), 1893.https://doi.org/10.3390/rs11161893 *
%
%2.        The Origional Matlab source code of GTrop model is available at * 
%          https://github.com/sun1753814280/GTrop.                         *
%==========================================================================+                        
%Codes Modified by: 
%                  Samuel Osah,Msc Geomatic Engineering,(PhD Student,KNUST)*    
%                  Email: osahsamuel@yahoo.ca                              *
%                  Tel:+233 (0)246137410/+233 (0)509438484                 *
%==========================================================================
%**************************************************************************+

%GET UTC TIME COMPONENTs
%******DATE[YEAR MONTH DAY]
Yr  = UTCtime(:,1);%get Hour
Mn  = UTCtime(:,2);%get Month
Day = UTCtime(:,3);%get Day

%******TIME[HOUR MINUTE SECONDs]
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
doy = DayofYear(Yr,Mn,Day);

%CHECK IF GTrop_Coeff IS EMPTY([]) & GET GRID VALUES/COEFFICIENTS
if isempty(GTrop_Coeff)
   
   %GET SAVED GTrop GRID VALUES
   if ~isempty(getappdata(0,'gridV_GTrop'))
            
       GTrop_Coeff = getappdata(0,'gridV_GTrop');  
       
   else
       %SEARCH FOR GTrop model GRID FILE(GTropCoefficient.mat) & DIRECTORY & LOAD VALUES               
       %Call the "SearchGTropgrid.m" fxn
       [~,GTrop_Coeff] = SearchGTropgrid();%****LOADED GTrop MODEL COEFFICIENTS
       
   end
   
end
       
%****ASSIGNMENT
coefficient = GTrop_Coeff;%GTrop grid values/coefficients

%NOTE:---------------------------------------------------------------------
% The size of the coefficients 181*361*37, The first two dimensions represent 
% the latitude and the longitude, respectively, and the third dimension stores
% the model coefficients and the reference height at each grid point
% 
% For the third dimension:
% 
% 1-6   model coefficients of ZHD
% 7-12  model coefficients of ZWD
% 13-18 model coeffiicents of Tm
% 19-24 model coefficients of ZHD lapse rate
% 25-30 model coefficients of ZWD lapse rate
% 31-36 model coefficients of Tm lapse rate
% 37    height of grid point unit: km
%--------------------------------------------------------------------------

%CONVERT HEIGHTS IN METERS TO KM
h  = h./1000; %[Km]

%GET THE LATITUDE AND LONGITUDE OF THE FOUR CORNER POINTS OF THE GRID 
%POINT WHERE THE LATITUDE AND LONGITUDE OF THE TARGET POINT IS LOCATED
%{I.E.OBTAIN THE LONGITUDES & LATITUDES OF THE 4 GRID POINTS CLOSEST TO THE 
%     GNSS RECEIVER}
B1 = floor(lat);
B2 = ceil(lat);
L1 = floor(lon);
L2 = ceil(lon);

B = lat;
L = lon;

%THE CASE WHERE THE GNSS RECEIVER IS AT THE GRID POINT
if B1 ==B2 && L1==L2
    
    i            = B1 + 91;
    j            = L1 + 181;
    a(1,:)       = coefficient(i,j,:);
    [zhd,zwd,tm] = GTropGrid(h,Yr,doy,a);

    %THE CASE WHERE THE LATITUDE OF GNSS RECEIVER AND THAT OF THE GRID
    %POINT ARE THE SAME====================================================+
elseif B1==B2
    
    i               = B1 + 91;
    j               = L1 + 181;
    a(1,:)          = coefficient(i,j,:);
    [zhd1,zwd1,tm1] = GTropGrid(h,Yr,doy,a);
    
    % ---------------------------------------------------------------------
    i               = B1 + 91;
    j               = L2 + 181;
    a(1,:)          = coefficient(i,j,:);
    [zhd2,zwd2,tm2] = GTropGrid(h,Yr,doy,a);
    
    % ---------------------------------------------------------------------
    
    %LINEAR INTERPOLATION ALONG THE LONGITUDE DIRECTION
    zhd = ((L2-L)*zhd1+(L-L1)*zhd2)/(L2-L1);    
    zwd = ((L2-L)*zwd1+(L-L1)*zwd2)/(L2-L1);    
    tm  = ((L2-L)*tm1+(L-L1)*tm2)/(L2-L1);    
   
    %THE CASE WHERE THE LONGITUDE OF GNSS RECEIVER AND THAT OF THE GRID
    %POINT ARE THE SAME====================================================+   
elseif L1==L2
   
    i               = B1 + 91;
    j               = L1 + 181;
    a(1,:)          = coefficient(i,j,:);
    [zhd1,zwd1,tm1] = GTropGrid(h,Yr,doy,a);
    
    % ---------------------------------------------------------------------
    
    i               = B2 + 91;
    j               = L1 + 181;
    a(1,:)          = coefficient(i,j,:);
    [zhd2,zwd2,tm2] = GTropGrid(h,Yr,doy,a);
    
    % ---------------------------------------------------------------------
    
    %LINEAR INTERPOLATION ALONG THE LATITUDE DIRECTION
    zhd = ((B2-B)*zhd1+(B-B1)*zhd2)/(B2-B1);    
    zwd = ((B2-B)*zwd1+(B-B1)*zwd2)/(B2-B1);    
    tm  = ((B2-B)*tm1+(B-B1)*tm2)/(B2-B1); 
   
    %THE CASE WHERE THE GNSS RECEIVER IS NOT AT THE GRID POINT  
else
    i               = B1 + 91;
    j               = L1 + 181;
    a(1,:)          = coefficient(i,j,:);
    [zhd1,zwd1,tm1] = GTropGrid(h,Yr,doy,a);
    
    % ---------------------------------------------------------------------

    i               = B1 + 91;
    j               = L2 + 181;
    a(1,:)          = coefficient(i,j,:);
    [zhd2,zwd2,tm2] = GTropGrid(h,Yr,doy,a);
    
    % ---------------------------------------------------------------------
    
    i               = B2 + 91;
    j               = L1 + 181;
    a(1,:)          = coefficient(i,j,:);
    [zhd3,zwd3,tm3] = GTropGrid(h,Yr,doy,a);
    
    % ---------------------------------------------------------------------
    
    i               = B2 + 91;
    j               = L2 + 181;
    a(1,:)          = coefficient(i,j,:);
    [zhd4,zwd4,tm4] = GTropGrid(h,Yr,doy,a);
    
    % ---------------------------------------------------------------------
    
    %BILINEAR INTERPOLATION ALONG THE LATITUDE & LONGITUDEDIRECTION
    zhd = ((B2-B)*(L2-L)*zhd1+(B-B1)*(L-L1)*zhd4+(B-B1)*(L2-L)*zhd3+(B2-B)*(L-L1)*zhd2)/((B2-B1)*(L2-L1));   
    zwd = ((B2-B)*(L2-L)*zwd1+(B-B1)*(L-L1)*zwd4+(B-B1)*(L2-L)*zwd3+(B2-B)*(L-L1)*zwd2)/((B2-B1)*(L2-L1));    
    tm  = ((B2-B)*(L2-L)*tm1+(B-B1)*(L-L1)*tm4+(B-B1)*(L2-L)*tm3+(B2-B)*(L-L1)*tm2)/((B2-B1)*(L2-L1));    
    
end

%ASSIGNMENT
Tm = tm;

%CONVERT zhd & zwd to m
ZHD = zhd./1000;%[m]
ZWD = zwd./1000;%[m]

%Compute Zenith Total Delay(ZTD)
ZTD = ZHD + ZWD;%[m]
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%(B.1.1)****SUB-ROUTINE TO COMPUTE TROPOSPHERIC PARAMETER AT A GRID POINT
function [zhd,zwd,tm] = GTropGrid(h,year,doy,a)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:   
%            This function "GTropGrid" is used to calcualte tropospheric...*
%            parameters at a grid point
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%USAGE:
%       [zhd,zwd,tm] = GTropGrid(h,year,doy,a);
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%1.    h    : Ellipsoidal height GNSS Receiver in km                       *
%2.    year : Year of observation                                          *
%3.    doy  : Day of Year                                                  *
%4.    a    : Extracted coefficient from grid point                        *

%OUTPUT:                                                                   *
%1.     zhd : Zenith Hdrostatic Delay at grid point                        *
%2.     zwd : Zenith Wet Delay at grid point                               *
%3.     tm  : Weighted Mean Temperature at grid point                      *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

Tr = year + doy/365.25 - 1980;
C1 = cos(doy*2*pi/365.25);
S1 = sin(doy*2*pi/365.25);
C2 = cos(doy*4*pi/365.25);
S2 = sin(doy*4*pi/365.25);

% GET THE ZHD,ZWD,AND Tm AT THE REFERENCE LEVEL AND THEIR LAPSE RATES
zhd0 = a(1)*Tr + a(2) + a(3)*C1 + a(4)*S1 + a(5)*C2 + a(6)*S2;
zwd0 = a(7)*Tr + a(8) + a(9)*C1 + a(10)*S1 + a(11)*C2 + a(12)*S2;
tm0  = a(13)*Tr + a(14) + a(15)*C1 + a(16)*S1 + a(17)*C2 + a(18)*S2;
zhdh = a(19)*Tr + a(20) + a(21)*C1 + a(22)*S1 + a(23)*C2 + a(24)*S2;
zwdh = a(25)*Tr + a(26) + a(27)*C1 + a(28)*S1 + a(29)*C2 + a(30)*S2;
tmh  = a(31)*Tr + a(32) + a(33)*C1 + a(34)*S1 + a(35)*C2 + a(36)*S2;
h0   = a(37);

% *****APPLY HEIGHT CORRECTION 
zhd = zhd0*(1 - zhdh*(h - h0))^5.225;
zwd = zwd0*(1 - zwdh*(h - h0))^5.225;
tm  = tm0 - tmh*(h - h0);
%=========================================END OF GTropGrid.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 



%              -----------------------------------------------------------
%(B.1.2)********SUB-ROUTINE FOR CHECKING THE EXISTENCE OF GTrop GRID FILE
%              ------------------------------------------------------------
function [GTrop_grid,gridV_GTrop,folder] = SearchGTropgrid()
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *                                              
%            The "SearchGTropgrid" subroutine searches Global Tropospheric +
%            Model(GTrop) Grid file in a given Directory/Folder as...      +
%            indicated in the sub-routine. Situation where GTrop file ...  + 
%            is not found in the default directory/folder,it searches ...  +   
%            recursively through all Sub-Directories /Folders of the given +  
%            Directory.GTrop file is extracted by looping through all ...  + 
%            the listed files in the provided folder or sub-folders.       +
%            Finally,if GTrop file is still not found in the default ...   + 
%            directory/folder and its sub-folders,the search for GTrop ... +
%            file is extended to the current directory.                    +           
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%%USAGE:                                                                   +
%       [GTrop_grid,gridV_GTrop,folder] = SearchGTropgrid()                +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%NOTE:                                                                     *
%     THE FF FOLDERS/DIRECTORIES ARE SEARCHED FOR GTrop GRID FILE           *
%1.   'GTrop file' folder; main folder for GTrop GRID                     *
%2.   'TropoGRIDS' folder; main folder for Tropospheric & Tm GRIDs         *
%3.   'data' folder; folder for goGPS data                                 *
%4.   'pwd'; Current folder/Directory. In goGPS, it is the 'goGPS' folder  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%--------------------------------------------------------------------------
%OUTPUTs                                                                   *
%1. GTrop_grid : Flag to indicate presence of GTrop grid file.1 to indicate*
%                 there is file in the folder & 0 to mean absence of file  *
%2.gridV_GTrop :  Extracted GTrop Grid Values                              *
%3.    folder  : The found directory or folder(File Path) for GTrop grid..* 
%                 file.e.g.:                                               *
%                 'C:...data\TropoGRIDS\GTrop file\GTropcoefficient.mat    *       
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
% WRITTEN BY: 
%            OSAH SAMUEL, MSC GEOMATIC ENGINEERING (PhD STUDENT)           +
%            Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com                +
%            Phone:+233(0)246137410 / +233(0)509438484                     +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%**************************************************************************+
    
%         -----------------------------------------------------
%*********SEARCH FOR GTrop GRID FILES FROM THE VARIOUS DIRECTORY
%         -----------------------------------------------------

%1.SEARCH GTrop GRIDs from 'GTrop file' folder
[GTrop_grid,GTropgDir] = SearchGTropGrid('../data/TropoGRIDS/GTrop file');

if isequal(GTrop_grid,0)%if GTrop grids are not in 'GTrop file' folder, try the 'TropoGRIDS'
    
   %2.SEARCH GTrop GRIDs from 'TropoGRIDS' folder 
   [GTrop_grid,GTropgDir] = SearchGTropGrid('../data/TropoGRIDS');

   if isequal(GTrop_grid,0)%if GTrop grids are not in 'TropoGRIDS' folder, try the 'data' folder
    
      [GTrop_grid,GTropgDir] = SearchGTropGrid('../data');%searching the 'data' folder
   
      if isequal(GTrop_grid,0) %if GTrop grids are not in the 'data' folder, try the 'current' directory
       
         [GTrop_grid,GTropgDir] = SearchGTropGrid(pwd);
      
         if isequal(GTrop_grid,0)%if GTrop grids are not in the 'current' directory, then issue a message
          
            beep %Give a beep sound
            errmsg2{1}=sprintf('No GTrop grid file found in goGPS directory.\n');
            errmsg2{2}=sprintf('Please Provide GTrop grid file & Try Again.\n');
            warndlg(errmsg2,'GTrop grid file Error','modal') 
            
            %ASSIGN EMPTY([]) MATRIX
            gridV_GTrop = [];
            folder      = [];
            
            return
      
         else  
             folder = GTropgDir;%GTrop grid directory/full path [EX.C:\Users\...data\TropoGRIDS\GTrop file\GTropcoefficient.mat]
      
         end  
      
      else  
          folder = GTropgDir; %GTrop grid directory/full path
       
      end  
   
   else   
       folder = GTropgDir; %GTrop grid directory/full path
       
   end 
       
else 
     folder = GTropgDir; %GTrop grid directory/full path  
    
end  

%            -------------------------------------------------
%************IF A DIRECTORY/FOLDER WITH GTrop GRID IS FOUND....
%            -------------------------------------------------
 if ~isempty(folder) %If Folder is not Empty
     
    %****LOAD TVGG-Tm MODEL COEFFICIENTS
    GridV_GTrop = load(folder);%181×361×37 ERA-INTERIM GRIDS
    gridV_GTrop = GridV_GTrop.coefficient;  
 end

 %SAVE DIRECTORY & GTrop GRID VALUES
 setappdata(0,'gridV_GTrop',gridV_GTrop)
 setappdata(0,'GTropgDir',folder)
 %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
 
%           -------------------------------------------------------------
%(B.1.2.1)**MAIN SUB-ROUTINE FOR CHECKING THE EXISTENCE OF GTrop GRID FILE
%           --------------------------------------------------------------

function [GTrop_grid,Folder] = SearchGTropGrid(directory)

%ASSIGNMENT
folder = directory;  %folder/Directory suppose to be containing GTrop FILE in goGPS

%****FIRST CHECK IF FOLDER IS EMPTY OR NOT
SizeListing=size(dir(folder)); %Size of folder content

%********CHECK IF FOLDER IS EMPTY OR NOT
if any([SizeListing(1,1)==2,SizeListing(1,1)==0])% size of 2 means folder is empty & 0 means folder doesn't exist
                  
   GTrop_grid = 0; %Flag to indicate the absence of GTrop grid file in directory
   Folder = []; %File directory is empty([]) if GTropCoefficient.mat is not found
   
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
          
       %*****GET INDEXEs FOR THE GTrop GRID FILE(i.e.Integer Val /[])
       GTropIndex = regexp(fileList,'\w*GTropCoefficient.mat'); %matches any words with GTropCoefficient.mat.  
          
       %************GET THE GRID FILES
       GTropFiles  = fileList(~cellfun(@isempty,GTropIndex));%GTrop files(eg:'GTropCoefficient.mat')
       
    else
        GTropFiles = [];
           
    end  %if ~isempty(fileList)
    
    %IF NO GTrop FILES EXIST IN THE CURRENT FOLDER, CHECK SUB-FOLDERS/              
    %SUB-DIRECTORIES FOR GTrop GRID FILEs
    %                              +                        ==========
     if isempty(fileList)
              
        checkSubfolders = 1; %Flag to search Sub-Directories/folders
           
     elseif isempty(GTropFiles)
                   
            checkSubfolders = 1;%Flag to search Sub-Directories/folders 
                         
     else
          if ~isempty(GTropFiles)
          
             checkSubfolders=0;%Flag not to search Sub-Directories/folders
             
             GTrop_grid = 1; %Flag to indicate the presence of GTrop grid file in directory
             
             %***CONSTRUCT THE FULL PATH [File directory [if GTropCoefficient.mat file is found)]       
             Folder = cellfun(@(x) fullfile(folder,x),GTropFiles,'UniformOutput',false);%Prepend path to files [File directory [if GTropCoefficient.mat file is found)] 
             
             %SORT FILES & REMOVE DUPLICATE FILES
             [GTropFiles,i_GTrop] = unique(GTropFiles);%GET GTrop FILES with no repetitions & INDEX(i_GTrop)
             Folder =(Folder(i_GTrop));%FILES WITH PATH
             
             %CONVERT CELL ARRAY TO CHARACTER
             if iscell(Folder)
                Folder=char(Folder);
             end    
             
             if iscell(GTropFiles)
                GTropFiles = char(GTropFiles);
             end    
             
          end
           
     end  %\\if isempty(fileList)
                
     
     if checkSubfolders==1 %if this condition is true,then either files in the current
                           %Directory were not GTrop grid file or there were no files at all.        
                           %System has to look through Sub-Folders or Sub-directories for GTrop Grid File                                                              
                                                                                                                                                                                                                                                                                  
        %***OPEN ALL SUB-FOLDERS & GET GTrop GRID FILE
                
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
            
                     j=1;%Loop index for the GTrop grid file type
                     
                      for i = 1:length(filelists)%****LOOP OVER THE FILES & CHECK FOR GTrop GRID FILE
                                                 %Identify the needed files and extract them separately
                                                 %into different file name(GTropfiles)
                                                 %Note that dir also lists the directories,so you have to check for them. 
                                                         
                           %***GET NAME OF FILES FROM THE LIST OF FILES FOR EACH LOOP               
                           
                           FileName = filelists{i,1};
                           filepath = filepaths{i,1};

                           %********MATCH ANY WORD WITH GTropCoefficient.mat
                           if regexpi(FileName,'\w*GTropCoefficient.mat') 

                             %When Tested & the file is a GTrop file
                             GTropfiles{j,1}=FileName; %#ok<*NASGU>
                             GTroppath{j,1}=filepath;%GTrop file Path in cells
                             
                             j=j+1;%Update index
                                                                           
                           end  %//if regexpi(FileName,'\w*GTropCoefficient.mat') 
                           
                      end %//for i = 1:length(fileLists) 
                                                                                                               
                  end   %//if ~isempty(fileLists)
                      
            else                          
                 %WRITE ERROR MESSAGE
                 erMsg{q,1}=sprintf('The Folder %s and  It''s Sub-Folder, %s , are empty / Contain No GTrop Grid file .\n',folder,subDirs{iDir});
                         
                 q=q+1; %update index
                        
            end   %if Sizesublist(1,1)~=2
                   
        end    %for iDir = find(validIndex)
        
      if ~exist('GTropfiles','var')
         
         GTrop_grid = 0;%Flag to indicate the absence of GTrop grid file in directory
         
         Folder = []; %File directory is empty([]) if GTropCoefficient.mat file is not found
         
         return
     
      else 
          GTrop_grid = 1; %Flag to indicate the presence of GTrop grid file in directory
          
          Folder = GTroppath; %if GTropcoefficient.mat file is found
          
          %SORT FILES & REMOVE DUPLICATE FILES
          [GTropFiles,i_GTrop] = unique(GTropfiles);%GET GTrop FILES with no repetitions & INDEX(i_GTrop)
           Folder =(Folder(i_GTrop));%FILES WITH PATH
         
          %CONVERT CELL ARRAY TO CHARACTER
          if iscell(Folder)
             Folder=char(Folder);
          end 
          
          if iscell(GTropFiles)
             GTropFiles = char(GTropFiles);
          end     
                    
      end %//~exist('GTropfiles','var')                                        
           
     end  %//if checkSubfolders==1
     
end %//if exist('listFiles','var')

%=========================================END OF SearchGTropGrid.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%                      END OF SearchGTropgrid
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%(C.0)*******SUBROUTINE TO COMPUTE NEILL MAPPING FUNCTION(NMF)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
function [NMFh,NMFw] = NMF(UTCtime,lat,hgt,satELEV)

%**************************************************************************
%***DESCRIPTION:
%              This subroutine determines the Global Mapping Functions GMF *
%              given receiver/station position in [lat long hgt],...       * 
%              Satellites Elevation Angle(el),and Modified Julian Day(mjd).*

%             Original source codesis by Böhm et al (2006), modified by ...* 
%             Osah Samuel.                                                 *                                                              
%USAGE:                                                                    *
%      General:                                                            *
%              [NMFh,NMFw] = UTCtime(UTCtime,lat,hgt,satEL)                *       
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%****INPUTs:                                                               *
%1.       UTCtime:.........UTC time in [Year,Month,Day,Hour,Minute,Seconds]*                    
%2.           lat:.........Station geodetic latitude in [Degrees]          *
%3.           hgt:.........Station ellipsoidal height in [m]               *
%4.         satEL:.........Satellites Elevation Angle in [degrees]         *
% 
%***OUTPUTs:                                                               *
%1.         NMFh:..........Neill Hydrostatic mapping function              *
%2.         NMFw:..........Neill Wet mapping function                      *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%GET SIZE OF USER INPUT TIME & POSITION
nrow_time=size(UTCtime,1);
nrow_pos=size(lat,1);

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
       
      %1.*****INITIALIZING OUTPUT VARIABLEs
      satel=zeros(size(satELEV,1),size(satELEV,2));%Assign zeros of nxm to satel
      [NMFh,NMFw]=deal(satel);%copy the contents of satel to NMFh,NMFw
      
       %FIND NUMBER OF ROWS & COLUMNs IN satELEV
       [nrow,ncol]=size(satELEV);
       
          for i=1:ncol %Loop over Station
    
             for j=1:nrow %Loop over Satellites Elevations

                %CONVERT ELEVATION ANGLEs TO RADIAN  
                %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
                %First, Check if the elevation is less than 0, set it to .1 deg
                %REMARKs:
                %       Valid elevations are between -pi/2 and pi/2. Elevations below 
                %       0.1 will have the a delay mapped to 0.1 degree.No mapping  
                %       below zero degree elevation is modeled.
						
                EL_zero = find(satELEV(j,i) < .0017);
                
                if ~isempty(EL_zero)
                    
                   satELEV(EL_zero) = ones(size(satELEV(EL_zero)))*.0017;
                   
                end  % if ~isempty(EL_zero)
      
                %CONVERT ELEVATION ANGLES IN DEGREES TO RADIANS
                satEL =  abs(satELEV(j,i)) * pi / 180; 
                
                %Call the "nmf.m" function
                [NMFh(j,i),NMFw(j,i)] = nmf(UTCtime(i,:),lat(i),hgt(i),satEL);
                
             end
             
          end
          
   else %for different times
       
      [nrow,ncol]=size(satELEV);%Get size of satELEV
       NMFh = cell(nrow,ncol,nrow_time);%Create a cell array of NMFh Output
       NMFw = deal(NMFh);%Create copy of NMFh 
       
      for k =1:nrow_time %Loop over UTCtime for different sets of time
          
         for i=1:ncol %Loop over Station
    
             for j=1:nrow %Loop over Satellites Elevations
            
            %CONVERT ELEVATION ANGLEs TO RADIAN  
            %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
            %First, Check if the elevation is less than 0, set it to .1 deg
            %REMARKs:
            %       Valid elevations are between -pi/2 and pi/2. Elevations below 
            %       0.1 will have the a delay mapped to 0.1 degree.No mapping  
            %       below zero degree elevation is modeled.
						
                EL_zero = find(satELEV(j,i) < .0017);
                
                if ~isempty(EL_zero)
                    
                   satELEV(EL_zero) = ones(size(satELEV(EL_zero)))*.0017;
                   
                end  % if ~isempty(EL_zero)
      
                %CONVERT ELEVATION ANGLES IN DEGREES TO RADIANS
                satEL =  abs(satELEV(j,i)) * pi / 180; 
                
                 %Call the "nmf.m" function
                [NMFh{j,i,k},NMFw{j,i,k}] = nmf(UTCtime(k,:),lat(i),hgt(i),satEL);
                                                             
             end
         end
      end
      
   end %//if isequal(Identical,1) 
   
else  
     %FIND NUMBER OF ROWS & COLUMNs IN satELEV
     [nrow,ncol]=size(satELEV);%Get size of satELEV  
     
     nrow_time=size(UTCtime,1); %Get the # of Rows in UTCtime 
      
     if nrow_time > 1 %(Indication for different sets of time)
          
       NMFh = cell(nrow,ncol,nrow_time);%Create a cell array of NMFh Output
       NMFw = deal(NMFh);%Create copy of NMFh 
       
       for k =1:nrow_time %Loop over sets if time
          
          for i=1:ncol %Loop over Station
    
          for j=1:nrow %Loop over Satellites Elevations
            
              %CONVERT ELEVATION ANGLEs TO RADIAN  
              %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
              %First, Check if the elevation is less than 0, set it to .1 deg
              %REMARKs:
              %       Valid elevations are between -pi/2 and pi/2. Elevations below 
              %       0.1 will have the a delay mapped to 0.1 degree.No mapping  
              %       below zero degree elevation is modeled.
						
              EL_zero = find(satELEV(j,i) < .0017);
                
              if ~isempty(EL_zero)
                    
                 satELEV(EL_zero) = ones(size(satELEV(EL_zero)))*.0017;
                   
              end   % if ~isempty(EL_zero)
      
              %CONVERT ELEVATION ANGLES IN DEGREES TO RADIANS
              satEL =  abs(satELEV(j,i)) * pi / 180; 
              
              %Call the "nmf.m" function 
              [NMFh{j,i,k},NMFw{j,i,k}] = nmf(UTCtime(k,:),lat(i),hgt(i),satEL); 
                
          end                             
          
          end 
          
       end 
       
     else %for a single time info
         %INITIALIZE OUTPUT
         satel=zeros(size(satELEV,1),size(satELEV,2));%Assign zeros of nxm to satel
         
         [NMFh,NMFw]=deal(satel);%copy the contents of satel to NMFh,NMFw

            
        for i=1:ncol %Loop over Station
    
           for j=1:nrow %Loop over Satellites 
              
              %CONVERT ELEVATION ANGLEs TO RADIAN  
              %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
              %First, Check if the elevation is less than 0, set it to .1 deg
              %REMARKs:
              %       Valid elevations are between -pi/2 and pi/2. Elevations below 
              %       0.1 will have the a delay mapped to 0.1 degree.No mapping  
              %       below zero degree elevation is modeled.
						
              EL_zero = find(satELEV(j,i) < .0017);
              if ~isempty(EL_zero)
                satELEV(EL_zero) = ones(size(satELEV(EL_zero)))*.0017;
              end   % if ~isempty(EL_zero)
      
             %CONVERT ELEVATION ANGLES IN DEGREES TO RADIANS 
             satEL =  abs(satELEV(j,i)) * pi / 180;  
             
             %Call the "nmf.m" function
             [NMFh(j,i),NMFw(j,i)] = nmf(UTCtime,lat(i),hgt(i),satEL);
                                                                            
           end
           
        end
        
     end %//if nrow_time > 1
                
end %//if isequal(nrow_time,nrow_pos)

%******************************************END OF NMF.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%(C.1)*******MAIN SUB-ROUTINE TO COMPUTE NEILL MAPPING FUNCTION(NMF)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
function [nmfh,nmfw] = nmf(UTCtime,lat,hgt,satEL) 

%***************************************************************************
%DESCRIPTION:                                                              * 
%            "nmf" is a subroutine that Computes the wet and dry ...       *
%             Niell mapping function[NMF] 1996                             *
%             The coefficients of the hydrostatic mapping function depend  * 
%             on the latitude and height above sea level of the receiver   *
%             station, and on the day of the year(DOY). The wet mapping    * 
%             function on the other hand, depends only on latitude.        *
%USAGE:                                                                    *
%      General:                                                            *
%              [nmfh,nmfw] = nmf(UTCtime,lat,hgt,satEL)                    *       
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%****INPUTs:                                                               *
%1.       UTCtime:.........UTC time in [Year,Month,Day,Hour,Minute,Seconds]*                         
%2.           lat:.........Station geodetic latitude in [degrees]          *
%4.           hgt:.........Station ellipsoidal height in [meters]          *
%5.         satEL:.........Satellites Elevation Angle in [radians]         *
% 
%***OUTPUTs:                                                               *
%1.         nmfh:..........Neill Hydrostatic mapping function              *
%2.         nmfw:..........Neill Wet mapping function                      *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-           
%==========================================================================+
%REFERENCE:                                                                *
%          Niell, A.E.(1996).Global mapping functions for the atmosphere   * 
%          delay at radio wavelengths. J. Geophys. Res.,101, pp.3227-3246  *
%==========================================================================+
%Written by: Samuel Osah,Msc Geomatic Engineering ,2016                    *    
%             Email: osahsamuel@yahoo.ca                                   *
%             Tel:+233 (0)246137410/+233 (0)509438484                      *
%==========================================================================
%**************************************************************************+
%GET UTC TIME COMPONENTs
%1.******UTCdate
  Yr = UTCtime(:,1);%get Hour
  Mn = UTCtime(:,2);%get Month
 Day = UTCtime(:,3);%get Day
%2.UTCtime
   H = UTCtime(:,4);%get Hour
 MIN = UTCtime(:,5);%get Minute
SECs = UTCtime(:,6);%get Seconds

%****DEFINE pi
pi = 3.14159265359d0;

%Coefficient for changing DOY to radian 
doy2rad = 2*pi/365.25d0; 
 
%LATITUDE ARRAY FOR THE HYDROSTATIC & WET MAPPING FUNCTION COEFFICIENTs
TableLat = [15,30,45,60,75]';

%*********************Initialize NEILL look-up table
%1.*******HYDROSTATIC MAPPING COEFFICIENTs
%AVERAGE OF COEFFICIENTs a,b and c CORRESPONDING 
%             TO THE GIVEN LATITUDE
abc_avg = [1.2769934  2.9153695  62.610505  %15
           1.2683230  2.9152299  62.837393  %30
           1.2465397  2.9288445  63.721774  %45
           1.2196049  2.9022565  63.824265  %60
           1.2045996  2.9024912  64.258455]*1.0e-3 ;%75
           %    a          b         c  

%SEASONAL VARIATION/AMPLITUDE OF COEFFICIENTs a,b and c CORRESPONDING 
%                  TO THE GIVEN LATITUDE
abc_amp =  [   0          0          0              %15
           1.2709626  2.1414979  9.0128400          %30 
           2.6523662  3.0160779  4.3497037          %45
           3.4000452  7.2562722  84.795348          %60
           4.1202191  11.723375  170.37206]*1.0e-5; %75
           %    a          b         c           
 
%HEIGHT CORRECTION COEFFICIENTs   
a_ht = 2.53d-5;
b_ht = 5.49d-3;
c_ht = 1.14d-3;
                                        
%2.***************WET MAPPING COEFFICIENTs
%******AVERAGE OF COEFFICIENTs a,b and c CORRESPONDING 
%                    TO THE GIVEN LATITUDE
abc_w2po = [5.8021897*1.0e-4  1.4275268*1.0e-3  4.3472961*1.0e-2    %15
            5.6794847*1.0e-4  1.5138625*1.0e-3  4.6729510*1.0e-2    %30 
            5.8118019*1.0e-4  1.4572752*1.0e-3  4.3908931*1.0e-2    %45
            5.9727542*1.0e-4  1.5007428*1.0e-3  4.4626982*1.0e-2    %60 
            6.1641693*1.0e-4  1.7599082*1.0e-3  5.4736038*1.0e-2];  %75 
            %    a                   b               c              
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%********INTERPOLATE AVERAGE & AMPLITUDE COEFFICIENT VALUEs 
%Call the "int_Neill_Mapping_COE.m" function
[avg,amp,abc_wet] = int_Neill_Mapping_COE(lat,TableLat,abc_avg,abc_amp,abc_w2po);

%***************COMPUTING MAPPING FUNCTIONs (NMFh,NMFw)
%GET Wet Mapping Coefficients
aw=abc_wet(:,1);%Assigning all rows,1st column to aw
bw=abc_wet(:,2);%Assigning all rows,2nd column to bw
cw=abc_wet(:,3);%Assigning all rows,3rd column to cw

%1.*****INITIALIZING OUTPUT VARIABLEs
%Mapping Functions
satel=zeros(size(satEL,1),size(satEL,2));%Assign zeros of nxm to satel
[nmfh,nmfw]=deal(satel);%copy the contents of satel to MFh,MFw

%FIND NUMBER OF ROWS & COLUMNs IN satEL
[nrow,ncol]=size(satEL);

%**************COMPUTE DAY OF YEAR(DOY)   
%Call the "utc2JulianDay_DoY.m" function
[~, ~,DOY]=utc2JulianDay_DoY(Yr,Mn,Day,H,MIN,SECs);
           
for i=1:ncol %Loop over Stations
    
   for j=1:nrow %Loop over Satellites 
      
   %*************COMPUTE HYDROSTATIC MAPPING COEFFICIENTs (a,b,c)     
    %************DEAL WITH SOUTHERN HEMISPHERE & YEARLY VARIATION
    %NOTE:  
         %Since no data from the southern hemisphere were used in developing
         %these mapping functions, the inversion of the seasons has been
         %accounted for simply by adding half a year to the phase for
         %southern latitudes. 
    if (lat(i) < 0) %if Southern Hemisphere
      DOY = DOY + 365.25/2; %Add half a year to the phase
    end 
    %***COMPUTE COS OF PHASE ANGLE
    COSphase = cos((DOY - 28) * doy2rad );

    %********NOW, COMPUTE HYDROSTATIC COFFICIENT VALUEs(ah,bh,ch)
    ah = avg(i,1) - amp(i,1) * COSphase;
    bh = avg(i,2) - amp(i,2) * COSphase;
    ch = avg(i,3) - amp(i,3) * COSphase;

    %*****COMPUTE sine ELEVATION ANGLE
    sine = sin(satEL(j,i));

    %********COMPUTE NEILL HYDROSTATIC MAPPING FUNCTION value
    beta   = bh/( sine + ch  );
    gamma  = ah/( sine + beta);
    topcon = (1+ ah/(1 + bh/(1 + ch)));
    nmfh   = topcon/(sine+gamma); 
    
    %*********APPLY HEIGHT CORRECTION     
    hs_km  = hgt(i)/1000.d0;%convert height to km
    beta         = b_ht/( sine + c_ht);
    gamma        = a_ht/( sine + beta);
    topcon       = (1 + a_ht/(1 + b_ht/(1 + c_ht)));
    ht_corr_coef = 1/sine - topcon/(sine + gamma);
    ht_corr      = ht_corr_coef * hs_km;
    nmfh(j,i)        = nmfh + ht_corr; 
    
    %********COMPUTE NEILL WET MAPPING FUNCTION value
    beta   = bw(i)/( sine + cw(i) );
    gamma  = aw(i)/( sine + beta);
    topcon = (1 + aw(i)/(1.d0 + bw(i)/(1 + cw(i))));
    nmfw(j,i)   = topcon/(sine+gamma);
      
   end  %for j=1:nrow
   
end %for i=1:nrow
%******************************************END OF nmf.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%(C.2)********SUBROUTINE TO INTERPOLATE NEILL MAPPING COEFFICIENTs
function [avg,amp,abc_wet] = int_Neill_Mapping_COE(UserLat,TableLat,HydavgValues...
                                                   ,ampValues,WetavgValues)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            The function "int_Neill_Mapping_CO" Interpolates and Computes * 
%            Hydrostatic & Wet mapping function coefficients(a,b,c) from...*  
%            averaged  and seasonal variation/Amplitude values used ...    * 
%            for tropospheric delay prediction given receiver latitude,... * 
%            based on Neill Mapping Function(NMF)Model.Parameters above... *    
%            |lat|<=15° and |lat|>=75° are extracted directly while ...    * 
%            Parameters for latitudes 15°<|lat|<75° are linearly...        *
%            interpolated between values of two  closest latitudes.        *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%       UserLat: Station/User geodetic latitude(s) in degrees(deg) (n x 1) *
%      TableLat: Latitude array for the Hydrostatic & Wet mapping function *
%                coefficients in degrees(deg).                             *
%  HydavgValues: Hydrostatic Average Values                                *
%     ampValues: Amplitude Values for Hydrostatic Mapping function         *
%  WetavgValues: Wet Average Values                                        *

%AVERAGE OF COEFFICIENTs a,b and c CORRESPONDING TO THE GIVEN LATITUDE     *
%HydavgValues =[1.2769934  2.9153695  62.610505         %15                *
%               1.2683230  2.9152299  62.837393         %30                *
%               1.2465397  2.9288445  63.721774         %45                *
%               1.2196049  2.9022565  63.824265         %60                *
%               1.2045996  2.9024912  64.258455]*1.0e-3;%75                *
%                   a          b         c                                 *
%AMPLITUDE OF COEFFICIENTs a,b and c CORRESPONDING TO THE GIVEN LATITUDE   *
%ampValues =  [   0          0          0               %15                *
%              1.2709626  2.1414979  9.0128400          %30                *
%              2.6523662  3.0160779  4.3497037          %45                *
%              3.4000452  7.2562722  84.795348          %60                *
%              4.1202191  11.723375  170.37206]*1.0e-5; %75                *
%                  a          b         c                                  *         
%2.********WET MAPPING FUNCTION                                            *
%WetavgValues =[5.8021897*1.0e-4  1.4275268*1.0e-3  4.3472961*1.0e-2  %15  *
%               5.6794847*1.0e-4  1.5138625*1.0e-3  4.6729510*1.0e-2  %30  *
%               5.8118019*1.0e-4  1.4572752*1.0e-3  4.3908931*1.0e-2  %45  *
%               5.9727542*1.0e-4  1.5007428*1.0e-3  4.4626982*1.0e-2  %60  *
%               6.1641693*1.0e-4  1.7599082*1.0e-3  5.4736038*1.0e-2];%75  *
%                    a                  b                c               
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%OUTPUT:                                                                   *
%       avg - Average(a,b,c) hydrostatic MF coefficients (nx3)             *
%       amp - seasonal  variation/Amplitude coefficients of a, b and c(nx3)*
%       abc_wet-coefficients of a, b and c(nx3)                            *
%WHERE:
%      nx3-->Means n Rows,three(3) Columns (e.g. 4x3)
%      n-->Indicates the Number of User/Receiver Input Coordinates(Latitude)
%          That is, each row in (avg,amp,abc_wet) represents Receiver Coord*
%     The Columns arrangements are as follows:
%          1       2       3       
%     avg:[a       b       c  ]
%     amp:[a       b       c  ]
% abc_wet:[a       b       c  ]
% =========================================================================
%REFERENCE:                                                                *
%          A.E.Niell, Global mapping functions for the atmosphere delay at * 
%          radio wavelengths, 1996                                         *
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************
%***************************************************************************
% Latitude=[15;16;25;30;35;45;60;75;80];
%*********ASSIGNMENT
abc_avg=HydavgValues;
abc_amp=ampValues;
abc_w2po=WetavgValues;

%****************INTERPOLATIONS  
%******INITIALIZE OUTPUT
 
 %***HYDROSTATIC
 avg=zeros(length(UserLat),size(abc_avg,2));%Assigning zeros to avg 
 amp=deal(avg);%copy the contents avg to amp
 %***WET
 abc_wet=deal(avg);%copy the contents avg to abc_wet
 
 for i=1:length(UserLat) %Loop through user given Latitude coordinates
     
     for j=1:length(TableLat) %Loop through Look-up Table
         
        if (j==1 & UserLat(i)<= TableLat(j))
           
           %***HYDROSTATIC
           avg(i,:) = abc_avg(j,:); %Average(a,b,c) hydrostatic MF coefficients
           amp(i,:) = abc_amp(j,:); %Amplitude of coefficients a, b and c
           
           %***WET
           abc_wet(i,:) = abc_w2po(j,:); %Wet coefficients of a, b and c 
           
        elseif (j==length(TableLat) & UserLat(i)>= TableLat(j))
               
               %***HYDROSTATIC
               avg(i,:) = abc_avg(j,:); %Average(a,b,c) hydrostatic MF coefficients
               amp(i,:) = abc_amp(j,:); %Amplitude of coefficients a, b and c
               
               %***WET
               abc_wet(i,:) = abc_w2po(j,:); %Wet coefficients of a, b and c 
        else     
            if (UserLat(i)== TableLat(j))
              
              %***HYDROSTATIC  
              avg(i,:) = abc_avg(j,:); %Average a hydrostatic MF coefficients
              amp(i,:) = abc_amp(j,:); %Amplitude of coefficients a, b and c
              
              %***WET
              abc_wet(i,:) = abc_w2po(j,:); %Wet coefficients of a, b and c 
              
            elseif (UserLat(i)>TableLat(j,1) && TableLat(j,1)<TableLat(j+1,1))
                  %********HYDROSTATIC
                  diff_avg= abc_avg(j+1,:)-abc_avg(j,:);%Difference in average Met parameters between two closest latitudes
                  diff_amp= abc_amp(j+1,:)-abc_amp(j,:);%Difference in seasonal variation of Met parameters between two closest latitudes
                  diff_latRMF=UserLat(i)-TableLat(j,1);%Difference in latitude values between user Lat and closest latitude of Met parameters
                  diff_latMF=TableLat(j+1,1)-TableLat(j,1);%Difference in latitude values between user Lat and closest latitude of Met parameters
                 
                  %********WET
                  diff_wet= abc_w2po(j+1,:)-abc_w2po(j,:);
                  
                  %***********PERFORM INTERPOLATION
                  %***HYDROSTATIC
                  avg(i,:) = ((diff_latRMF/diff_latMF).*diff_avg)+abc_avg(j,:); %Average Meteorological parameters
                  amp(i,:) = ((diff_latRMF/diff_latMF).*diff_amp)+ abc_amp(j,:);%Seasonal variation Meteorological parameters          
                 
                  %***WET
                  abc_wet(i,:) = ((diff_latRMF/diff_latMF).*diff_wet)+ abc_w2po(j,:);
                  
            end %if (latR(i)== lat_MF(j))
        end %if (j==1 & latR(i)<= lat_MF(j))
        
       if (j+1)>length(TableLat)
             
          break
       end  %if j+1>length(lat_MF) 
        
     end %for j=1:length(lat_MF)
     
 end %for i=1:length(latR)
%******************************END OF int_Neill_Mapping_COE.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%                 *END OF NEILL MAPPING FUNCTION*
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%(D) SUBROUTINE TO COMPUTE JULIAN, MODIFIED DAY(JD,MJD) &  DAY OF YEAR(DoY)
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
          
if any(UT(:) ~=0)
    
  try
     JD = floor( (365.25 * Y) ) + floor( (30.6001 * (M+1)) ) + D + ( UT / 24 ) + 1720981.5;
  catch      
       JD = floor( (365.25 * Y) ) + floor( (30.6001 * (M+1)) )- 15 + 1720996.5 + D + ( UT / 24 );    
  end   %try
                   
  %****CALCULATION OF MODIFIED JULIAN DATE(MJD)
  MJD = JD-2400000.5;%Modified Julian date

else
    JD=juliandate(Yr, Mn, D);
    MJD=mjuliandate(Yr, Mn, D);
end

%****COMPUTE DAY OF YEAR
%Call the "DayofYear function.m"
 DoY=DayofYear(Year,Month,Day);
%******************************************END OF utc2JulianDay_DoY.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
 
%(D.1)*********SUBROUTINE TO COMPUTE DAY OF YEAR(DoY)
function DoY=DayofYear(Year,Month,Day)

switch nargin    
    case 1
          if size(Year,2)==3
             Yr=Year(:,1);%Assigning 1st column to years
             Mn=Year(:,2);%Assigning 2ND column to months
             D=Year(:,3);%Assigning 3RD column to days             
          end          
    otherwise        
             Yr=Year;%Assigning Year to Yr
             Mn=Month;%Assigning Month to Mn
             D=Day;%Assigning Day to D
             
end %switch nargin
                          
%******CALCULATION OF DAY OF YEAR(DoY)
try
   I = ones(length(Yr),1);
 DoY = datenum([Yr Mn D]) - datenum([Yr I I]) + 1;
catch    
     %Compute Julian Date
     jDay  = juliandate(Yr,Mn,D,0,0,0);%Midnight Morning of this Day
  jDayJan1 = juliandate(Yr,1,1,0,0,0);%Midnight Morning of January 1st
 dayNumber = jDay - jDayJan1 + 1;
       DoY = dayNumber;%Day of the Year     
end %try
%******************************************END OF DayofYear.m 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%(D.2)*********SUBROUTINE TO CONVERT TWO DIGITs TO YEAR 4 DIGITS
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
%1.******SUBROUTINE TO CONVERT XYZ COORDs TO GEOGRAPHIC COORDs
function [X,Y,Z] = geo2xyz(latitude,longitude,height,RefEllipsoid)

% *************************************************************************
% *************************************************************************
%DESCRIPTION:                                                              *
%           geo2xyz  calculate Cartesian coordinates(XYZ) given geodetic 
%           coordinates Latitude(degree/dms,Longitude(degree/dms,and ...   * 
%           height(m)above reference ellipsoid(RefEllipsoid) given by the  * 
%                                                                     user.* 

%The function Calls the following function(s):                             *
%1. [a,finv] = Elipsoidpara(RefEllipsoid);                                 *

%INPUT:                                                                    *
%1      latitude -> latitude of point in decimal degrees or DMS            * 
%2.     longitude -> longitude of point in decimal degrees or DMS          *
%3.     height -> ellipsoidal height in meters(m)                          * 
%4.     RefEllipsoid-Reference Ellipsoid in single quote(eg:'WGS84',...    *
%                                                       'GRS80','Pz-90')   *
%NOTE:                                                                     * 
%  if there is no height put zero(0)                                       *
% =========================================================================
%OUTPUT:                                                                   *
%1.    X --> 3D Cartesian X coordinates in meters                          *
%2.    Y --> 3D Cartesian Y coordinates in meters                          *
%3.    Z-->  3D Cartesian Z coordinates in meters                          *
% =========================================================================
%REFERENCE:*
%1         B.Hofmann-Wellenhof, H.Lichtenegger and J.Collins: GPS Theory   *
%          and practice. 2001. Fifth revised edition. Springer,Wien,New ... 
%                                                           York.p.280-282 *
%2        GILBERT STRANG AND KAI BORRE: Linear Algebra,Geodesy,and GPS     *
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%**************************************************************************
%**************************************************************************

%****FORMAT DATA INPUT
switch nargin 
    
    case 1 
         RefEllipsoid='WGS84'; 
         lat=latitude;%Assigning latitude to lat
         
         %***GET THE NUMBER OF COLUMNS IN lat
         clat=size(lat,2);%No. of columns of lat matrix

         if clat==7
           latDeg  = dms2degrees(lat(:,1:3));%Latitude in DMS to decimal degrees
           longDeg = dms2degrees(lat(:,4:6));%Longitude in DMS to decimal degrees
           h       = lat(:,7);%Ellipsoidal heights in the 7TH Column
         elseif  clat==6   
               latDeg  = dms2degrees(lat(:,1:3));%Latitude in DMS to decimal degrees
               longDeg = dms2degrees(lat(:,4:6));%Longitude in DMS to decimal degrees
               h       = zeros(size(latDeg,1),1);%%Assigning zeros to Ellipsoidal height(h) 
         elseif  clat==5    
               latDeg  = dm2degrees(lat(:,1:2));%Latitude in DM to decimal degrees
               longDeg = dm2degrees(lat(:,3:4));%Longitude in DM to decimal degrees
               h       = lat(:,5);%Ellipsoidal heights in the 5TH Column
         elseif  clat==4  
               latDeg  = dm2degrees(lat(:,1:2));%Latitude in DM to decimal degrees
               longDeg = dm2degrees(lat(:,3:4));%Longitude in DM to decimal degrees
               h       = zeros(size(latDeg,1),1);%%Assigning zeros to Ellipsoidal height(h)    
         elseif  clat==3   
               latDeg  = lat(:,1);%Latitude in decimal degrees in the 1ST Column
               longDeg = lat(:,2);%Longitude in decimal degrees in the 2ND Column
               h       = lat(:,3);%Ellipsoidal heights in the 3RD Column
         elseif  clat==2  
               latDeg  = lat(:,1);%Latitude in decimal degrees in the 1ST Column
               longDeg = lat(:,2);%Longitude in decimal degrees in the 2ND Column
               h       = zeros(size(latDeg,1),1);%%Assigning zeros to Ellipsoidal height(h)
         end    %if clat==7

    case {4,3,2} 
         lat  = latitude;%Assigning latitude to latDeg
         long = longitude;%Assigning longitude to lonitude
         
         %%***GET THE NUMBER OF COLUMNS IN lat & long 
         clat = size(lat,2);%No. of columns of latDeg matrix
         clong= size(long,2);%No. of columns of longDeg matrix
         
         if clat==3 %if No. of columns is 2 then data is in DMS   
           latDeg=dms2degrees(lat);%Latitude in DMS to decimal degrees
         elseif   clat==2 %if No. of columns is 2 then data is in DM    
               latDeg=dm2degrees(lat);%Latitude in DM to decimal degrees
         else   %if No. of columns is 1 then data is in decimal degrees
             latDeg=lat;   
         end     %if clat==3  
    
           if clong==3 %if No. of columns is 2 then data is in DMS  
             longDeg=dms2degrees(long);%Latitude in DMS to decimal degrees
    
           elseif   clong==2 %if No. of columns is 2 then data is in DM   
                 longDeg=dm2degrees(long);%Latitude in DM to decimal degrees

           else   %if No. of columns is 1 then data is in decimal degrees 
               longDeg=long;   
           end   %if clong==3 
           
           if (any(nargin==[4,3,2]))
            
          if nargin ==4 
            h = height;%Assigning height to h
          elseif nargin ==3
                RefEllipsoid='WGS84';
                h = height;%Assigning height to h
          elseif nargin ==2
                RefEllipsoid='WGS84';
                h = zeros(size(latitude,1),1);%Assigning zeros to ...
                                              %       Ellipsoidal height(h)
          end %if nargin ==4
          
           end %if (any(nargin==[4,3,2]))
           
           %CHECK IF SIZE OF h CONFORMS WITH SIZE OF lat/long
           if (size(h,1)==1 && (size(lat,1) && size(long,1))>1)         
            if h==0       
              h=repmat(h,size(lat,1),1);%Replicating height to conform with 
                                        %              size of lat or long   
            else
                h=[h;repmat(zeros(size(lat,1)-1,1),1)];%Replicating height to 
                                                      %conform with size of 
                                                      %         lat or long     
            end %if h==0
            
           elseif (size(h,1)>1 && (size(lat,1) && size(long,1))>1)
               
                 if size(h,1)<size(lat,1)
                   h=[h;repmat(zeros(size(lat,1)-size(h,1),1),1)];%Replicating height  
                                                                 %to conform  
                                                                 %with size of 
                                                                 %lat or long
                 end %if size(h,1)<size(lat,1)
                                                                                                                                                         
           end   %if (size(h,1)==1 && (size(lat,1)&& size(long,1))>1)
       
end   % switch nargin

%ISSUE ERROR MESSAGE IF SIZE OF lat & long ARE UNEQUAL

if size(lat,1)>size(long,1) %Check sizes
   beep%Give a beep sound 
   errmsg{1}=sprintf('Size of Latitude Coordinates  >  Size of Longitude Coordinates   %d > %d .',size(lat,1),size(long,1));
   errmsg{2}='';
   errmsg{3}='Please ensure that size of Latitude & Longitude Coordinates are the same .';
   errordlg(errmsg,'Coordinate(s) Input Error','modal')   
   return 
   
elseif size(lat,1)<size(long,1) %Check sizes 
      beep%Give a beep sound 
      errmsg{1}=sprintf('Size of Latitude Coordinates  <  Size of Longitude Coordinates  %d > %d .',size(lat,1),size(long,1));
      errmsg{2}='';
      errmsg{3}='Please ensure that size of Latitude & Longitude Coordinates are the same .';
      errordlg(errmsg,'Coordinate(s) Input Error','modal')   
      return 
      
end %if size(lat,1)>size(long,1) %Check sizes

%******CALCULATION 
%            
%INITIALIZE OUTPUT
X=zeros(size(latDeg,1),1);
Y=deal(X);
Z=deal(X);

%****GET ELLIPSOID PARAMETERS
%      
[a,finv] = Ellipsoidpara(RefEllipsoid);

f=1/finv;%flattening

%***compute square of eccentricity
e2 = (2.*f)-((f).^2);% first eccentricity

for i=1:size(latDeg,1)
    
%***CONVERT LatDeg & LongDeg to Radians 
latRad=deg2rad(latDeg(i));%converting latitude values in degrees to radian
longRad=deg2rad(longDeg(i));%converting longitude values in degrees to radian

%***COMPUTE RADIUS OF PRIME VERTICAL(N)
N = a./sqrt(1-(e2.*(sin(latRad)).^2));%prime vertical radius

%***COMPUTE DISTANCE FROM Z-AXIS(P) 
P = (N + h(i)).*cos(latRad);

%COMPUTE XYZ COORDINATES
Z(i,1) = roundn((((1-e2).*N + h(i)) .* sin(latRad)),-3); %cartesian Z coordinates in meters
X(i,1) = roundn(P.*cos(longRad),-3); %cartesian X coordinates in meters
Y(i,1) = roundn(P.*sin(longRad),-3); %cartesian Y coordinates in meters

% fprintf('\n X =%12.3f\n Y =%12.3f\n Z =%12.3f\n',X(i),Y(i),Z(i));

end
%%%%%%%%%%%%%%%%%%%%%%%%END OF geo2xyz.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

%2.******SUBROUTINE TO CONVERT XYZ COORDs TO GEOGRAPHIC COORDs
function [latRAD,longRAD,h,latDEC,longDEC] = xyz2LLH(X,Y,Z,RefEllipsoid)
% *************************************************************************+
%DESCRIPTION:                                                              +
%           xyz2LLH Converts Cartesian(XYZ)Coordinates to Geographic/...
%                                                    Ellipsoidal coords    +
%           (lat,lon,alt) on WGS-84 according to a non iterative method...
%                                                         (Bowring method) +

%INPUT:                                                                    +
%1       X                                                                 + 
%2.      Y  > vectors of Cartesian(ECEF) Coordinates  (m)                  +
%3.      Z                                                                 + 
%4.      RefEllipsoid-Reference Ellipsoid in single quote(eg:'WGS84','GRS80+
%                                                                ','Pz-90')  
% =========================================================================+
%OUTPUT:                                                                   +
%1.    latRAD  --> latitude in radians                                     +
%2.    longRAD --> longitude in radians                                    +
%3.          h --> ellipsoidal height in units like the input              +
%4.    latDEC  --> latitude in decimal degrees                             +
%5.    longDEC --> longitude in decimal degrees                            +
% =========================================================================+
%REFERENCE:
%1         B.Hofmann-Wellenhof, H.Lichtenegger and J.Collins: GPS Theory   +
%          and practice. 2001. Fifth revised edition. Springer,Wien,...    +
%                                                       New York.p.280-282 +
%2         GILBERT STRANG AND KAI BORRE: Linear Algebra,Geodesy,and GPS    
% =========================================================================+
%  Written by OSAH SAMUEL, Msc Geomatic Engineering ,2016                                 
%       Email: osahsamuel@yahoo.ca                                                       
%         Tel:+233 (0)246137410/+233 (0)509438484                                         
%==========================================================================+
%**************************************************************************+
%**************************************************************************+
%****FORMAT DATA INPUT

switch nargin      
    case 1
          RefEllipsoid='WGS84';
          XYZ=X  ;       
         [nrow,ncol]=size(XYZ);

         if (ncol==3 && (nrow>ncol || nrow<ncol ))
           X=XYZ(:,1);
           Y=XYZ(:,2);
           Z=XYZ(:,3) ;   
         elseif  ((nrow==3) && (ncol>nrow || ncol<nrow ))    
               X=(XYZ(1,:))';
               Y=(XYZ(2,:))';
               Z=(XYZ(3,:))'; 
         elseif  ((nrow==3) && (ncol==3 ))   
              %***Round X to 2 digit numbers
              X_vpa=vpa(roundn(XYZ,5),2);%Variable precision arithmetic(vpa).
              Xm=mean(X_vpa);%Mean
              if (find(X_vpa(:,1)-Xm(1)==0)| find(X_vpa(:,2)-Xm(2)==0)| find(X_vpa(:,3)-Xm(3)==0))   
                X=XYZ(:,1);
                Y=XYZ(:,2);
                Z=XYZ(:,3);
              elseif (find(X_vpa(:,1)-Xm(1)~=0)| find(X_vpa(:,2)-Xm(2)~=0)| find(X_vpa(:,3)-Xm(3)~=0))    
                    X=(XYZ(1,:))';
                    Y=(XYZ(2,:))';
                    Z=(XYZ(3,:))';  
              end             
         end 
 case 2
RefEllipsoid=Y;     

%***FORMAT DATA
 XYZ=X  ;       
[nrow,ncol]=size(XYZ);
if (ncol==3 && (nrow>ncol || nrow<ncol ))
X=XYZ(:,1);
Y=XYZ(:,2);
Z=XYZ(:,3) ;   
elseif ((nrow==3) && (ncol>nrow || ncol<nrow ))     
X=(XYZ(1,:))';
Y=(XYZ(2,:))';
Z=(XYZ(3,:))';      
elseif ((nrow==3) && (ncol==3 ))   
    
%***Round X to 2 digit numbers
X_vpa=vpa(roundn(XYZ,5),2);%Variable precision arithmetic(vpa).
Xm=mean(X_vpa);%Mean
if (find(X_vpa(:,1)-Xm(1)==0)| find(X_vpa(:,2)-Xm(2)==0)| find(X_vpa(:,3)-Xm(3)==0))   
X=XYZ(:,1);
Y=XYZ(:,2);
Z=XYZ(:,3);
elseif (find(X_vpa(:,1)-Xm(1)~=0)| find(X_vpa(:,2)-Xm(2)~=0)| find(X_vpa(:,3)-Xm(3)~=0))    
X=(XYZ(1,:))';
Y=(XYZ(2,:))';
Z=(XYZ(3,:))';  
end         
end
 case 3
     
RefEllipsoid='WGS84'; 
 
end   
        
%******CALCULATION 
%            
%INITIALIZE OUTPUT
latRAD=zeros(size(X,1),1);
longRAD=deal(latRAD);
latDEC=deal(latRAD);
longDEC=deal(latRAD);
h=deal(latRAD);

%****GET ELLIPSOID PARAMETERS
%      
[a,finv] = Ellipsoidpara(RefEllipsoid);

f=1/finv;% flattening
b=a*(1-f);%radus of semi-minor axis
e2=(2*f)-f.^2; %eccentrcity squared
se2=e2./(1-e2) ; % second eccentricity squared

for i=1:size(X,1)
%***COMPUTE DISTANCE IN XY PLANE FROM GEOCENTER 
p=sqrt((X(i)).^2+(Y(i)).^2);
 
try
theta=atan((a.*Z(i))./(b.*p));
catch
    theta = atan2(Z(i).*a,p.*b);
end

%**COMPUTE LONGITUDE OF POINT IN RADIANS
longRAD(i,1)=atan2(Y(i),X(i));%longitude of points
% Format the longitude value
if (X(i) < 0) && (Y(i) > 0)
    longRAD(i,1)= pi + longRAD(i,1);
elseif (X(i) < 0) && (Y(i) < 0)
    longRAD(i,1) = longRAD(i,1) - pi;
end

%**COMPUTE LATITUDE OF POINT IN RADIANS
latRAD(i,1)=atan((Z(i)+(se2.*b.*((sin(theta)).^3)))./(p-(e2*a*((cos(theta)).^3))));%latitude of points

%****COMPUTE ELLIPSOIDAL HEIGHT(h)
V = a./sqrt(1-(e2.*(sin(latRAD(i,1))).^2)); % prime vertical radius of curvature(v)
h(i,1)=(p./cos(latRAD(i,1)))-V;% ellipsoidal height

%****CONVERT LATITUDE & LONGITUDE IN RADIAN TO DECIMAL DEGREES
latDEC(i,1) = rad2deg(latRAD(i,1));
longDEC(i,1)= rad2deg(longRAD(i,1));

end
%%%%%%%%%%%%%%%%%%%%%%%%END OF xyz2LLH.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%3.******SUBROUTINE TO CONVERT XYZ COORDs TO GEOGRAPHIC COORDs
function [a,finv,RefEllipsoid] = Ellipsoidpara(RefEllipsoid)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *
%            Elipsoidpara Returns the radius of the Semi Major Axis(a) and * 
%            inverse flattening of a given Ellipsoid                       * 
%INPUT:                                                                    *
%      RefEllipsoid -The name of the Ellipsoid in single quote(eg:'WGS84') *

%OUTPUT:                                                                   *
%       a    = Semi Major Axis of the given ellipsoid in meters            *
%       finv = inverse flattening of the given ellipsoid                   *

%REFERENCE:                                                                *
%         Department of Defence World Geodetic System 1984                 *
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *                                       
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************
%***************************************************************************

if ischar(RefEllipsoid)
 
  refE=RefEllipsoid;%Refernce Ellipsoid

 %***GPS REFERENCE FRAME
if (strcmpi(refE,'WGS84')|| (strcmpi(refE,'WGS1984')||strcmpi(refE,'GPS')))    
RefEllipsoid='WGS84';
a=6378137.0;%radius of major axis on the WGS 84 ellipsoid
finv=298.257223563;%inverse flattening on the WGS 84 ellipsoid

elseif (strcmpi(refE,'WGS72')|| strcmpi(refE,'WGS1972') )     
RefEllipsoid='WGS72';
a=6378135.0;%radius of major axis on the WGS 84 ellipsoid
finv=298.26;%inverse flattening on the WGS 84 ellipsoid

%***Geodetic Reference System 1980
elseif (strcmpi(refE,'GRS80')|| strcmpi(refE,'GRS1980') )     
RefEllipsoid='GRS80'; 
a=6378137.0;%radius of major axis on the WGS 84 ellipsoid
finv=298.257222101;%inverse flattening on the WGS 84 ellipsoid

%Geodetic Reference System 1967
elseif (strcmpi(refE,'GRS67') || strcmpi(refE,'GRS1967'))
RefEllipsoid='GRS67';     
a=6378160.0;%radius of major axis
finv=298.247167427;%inverse flattening
 
%GLONASS REFERENC FRAME
elseif (strcmpi(refE,'Pz-90')||(strcmpi(refE,'Pz-90.02')||strcmpi(refE,'GLONASS')))     
RefEllipsoid='Pz-90'; 
a=6378136.0;%radius of major axis on the WGS 84 ellipsoid
finv=298.257839303;%inverse flattening on the WGS 84 ellipsoid

%***BEIDOU REFERENC FRAME
elseif (strcmpi(refE,'CGCS2000')|| strcmpi(refE,'BEIDOU') )     
RefEllipsoid='CGCS2000'; 
a=6378137.0;%radius of major axis on the WGS 84 ellipsoid
finv=298.257222101;%inverse flattening on the WGS 84 ellipsoid

%***GHANA
%OLD
elseif (strcmpi(refE,'WAROFFICE') || strcmpi(refE,'WO')) 
RefEllipsoid='WAROFFICE';    
a=6378299.996;%radius of major axis
finv=296;%inverse flattening

elseif (strcmpi(refE,'ACCRADATUM') || strcmpi(refE,'AD'))
RefEllipsoid='WAROFFICE';  
a=6378299.996;%radius of major axis
finv=296;%inverse flattening

elseif (strcmpi(refE,'LEIGONDATUM') || strcmpi(refE,'LD'))
RefEllipsoid='Clarke1880';  
a=20926202 ;%radius of major axis
finv=293.465;%inverse flattening

%NEW DATUM(Ghana Geodetic Datum)    
elseif (strcmpi(refE,'GGD') || strcmpi(refE,'Ghana Geodetic Datum'))
RefEllipsoid='GRS80'; %Geodetic Reference System 1980
a=6378137.0;%radius of major axis on the WGS 84 ellipsoid
finv=298.257222101;%inverse flattening on the WGS 84 ellipsoid
    
elseif (strcmpi(refE,'Clarke1866') || strcmpi(refE,'Clarke66'))
RefEllipsoid='Clarke66';      
a=6378206.4;%radius of major axis
finv=294.9786982; %inverse flattening  

elseif (strcmpi(refE,'Clarke1880') || strcmpi(refE,'Clarke80'))
RefEllipsoid='Clarke80';      
a=6378249.145;%radius of major axis
finv=293.465;%inverse flattening    

elseif (strcmpi(refE,'Airy1830') || strcmpi(refE,'Airy30'))
RefEllipsoid='Airy30';    
a=6377563.396;%radius of major axis
finv=299.3249646;%inverse flattening

elseif (strcmpi(refE,'Everest30') || strcmpi(refE,'Everest1830'))
RefEllipsoid='Everest30';    
a=6377276.0;%radius of major axis
finv=300.8;%inverse flattening
  
elseif (strcmpi(refE,'Bessel41') || strcmpi(refE,'Bessel1841'))
RefEllipsoid='Bessel41';     
a=6377397.0;%radius of major axis
finv=299.15;%inverse flattening
end

else
RefEllipsoid='WGS84';
a=6378137.0;%radius of major axis on the WGS 84 ellipsoid
finv=298.257223563;%inverse flattening on the WGS 84 ellipsoid    
end   
%%%%%%%%%%%%%%%%%%%%%%%%END OF Elipsoidpara.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%4.******SUBROUTINE TO COMPUTE SATELLITE AZIMUTH & ELEVATION ANGLE

function [satEL_deg,satEL_rad,satAZ_deg,satAZ_rad,Range]=satAzEl(UserXYZ,SatXYZ,...
                                                              RefEllipsoid)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *
%            "satAzEl" Computes Satellites Azimuth(AZ),Elevation Angle(EL) *
%             and Range(D),relative to observer position/Ground Station /..*
%             from User to Satellite given User and Satellite Position(s)  *   
  
%The function Calls the following Subroutine(s):                           *
%1. [a,finv] = Elipsoidpara(RefEllipsoid);                                 *
%2. [latRAD,longRAD,h] = xyz2LLH(X,Y,Z,RefEllipsoid);                      *

%INPUT:                                                                    *
%      UserXYZ - Observer position/Ground station coordinates(ECEF)in(m) ..*
%                                                  (3 x n) / (n x 3)matrix *
%      SatXYZ  - ECEF satellite coordinates (m), (3 x n) / (n x 3) matrix  *

%     I.e.: 3xn=|[X eg:[6332942.597  6332976.932  6332957.890  6332977.582|*
%               | Y    -172955.641  -172805.004  -172878.972   -172804.786|*
%               | Z]    737935.003   737647.856   737824.057   737648.519]|*
%                ----------------------------------------------------------*
%           nx3=|[ X Y Z] eg:[6332942.597  -172955.641  737935.003 |       *
%               |             6332976.932  -172805.004  737647.856 |       *
%               |             6332957.890  -172878.972  737824.057 |       *
%               |             6332977.582  -172804.786  737648.519]|       *
%                --------------------------------------------------        *
%      RefEllipsoid-Reference Ellipsoid in single quote(eg:'WGS84','GRS80',*
%                                                                   'Pz-90)*
%NOTE:
%      If RefEllipsoid is not indicated,system uses 'WGS84' as default     *

%OUTPUT:                                                                   *
%       satAZ_rad = Azimuth (radians ClockWise from North)                 *
%		satEL_rad = Elevation Angle (radians)                              *
%       satAZ_deg = Azimuth (degrees ClockWise from North)                 *
%		satEL_deg = Elevation Angle (degrees)                              *
%		    Range = Receiver-Satellite Distance in units like the input    *

% REMARKS:                                                                 *
%         All the catch routins in this m-file are alternative solutions to*
%         their corresponding try(s).all give the same results.                               *

% =========================================================================
%REFERENCE:                                                                *
%1         B.Hofmann-Wellenhof, H.Lichtenegger and J.Collins: GPS Theory   *
%          and practice. 2001. Fifth revised edition. Springer, Wien, New..*
%                                                            York.p.280-282*
%2         GPS Theory and application",edited by B.Parkinson,J.Spilker     *
%3.        Alfred Leick, GPS Satellite Surveying, 2nd ed.,...              *
%          Wiley-Interscience,John Wiley & Sons, New York, 1995.           *                     * 
%4.        GILBERT STRANG AND KAI BORRE: Linear Algebra,Geodesy,and GPS    *
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *                                      
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************
%***************************************************************************

%****FORMAT DATA INPUT

switch nargin
    
    case 2
               
RefEllipsoid='WGS84';
          
end

[nrow,ncol]=size(UserXYZ);

if(ncol==3 && (nrow>ncol || nrow<ncol ))

UserXYZ=UserXYZ';

end
[nrow,ncol]=size(SatXYZ);

if(ncol==3 && (nrow>ncol || nrow<ncol ))

SatXYZ=SatXYZ';

end

%INITIALIZING OUTPUT VARIABLEs
satAZ_rad=zeros(size(SatXYZ,2),size(UserXYZ,2));%Assign zeros of nxm to satAZ_rad
[satEL_rad,satAZ_deg,satEL_deg,Range]=deal(satAZ_rad);%copy the contents of satAZ_rad  
                                                  %to all the requested outputs
Az_rad=zeros(size(SatXYZ,2),1);%Assign zeros of nx1 to Az_rad
[El_rad,Az_deg,El_deg,hdist,range]=deal(Az_rad);%copy the contents of Az_rad  
                                                %to all the requested outputs
for j=1:size(UserXYZ,2)%LOOP THROUGH ALL USER POSITIONs
    
for i=1:size(SatXYZ,2)%LOOP THROUGH ALL SATELLITE POSITIONs

%**CONVERT USER POSITION(XYZ) TO latitude and longitude 
[lat,long] = xyz2LLH(UserXYZ(:,j),RefEllipsoid);

%FINDING Baseline differences
dxyz = SatXYZ(:,i)-UserXYZ(:,j);

try

%*****COMPUTING Rotation matrix
%REMARKS:
%      R is direction cosine matrix to transform position vector components 
%      from geocentric equatorial frame into the topocentric horizon frame.

try
R = [-sin(long)  -sin(lat)*cos(long)  cos(lat)*cos(long);
      cos(long)  -sin(lat)*sin(long)  cos(lat)*sin(long);
       0              cos(lat)            sin(lat)     ];
 %**COMPUTE UNIT VECTORS   
  vENU = R'*dxyz; 
   
catch
   
 R = [   -sin(long)                        cos(long)         0     ;  
              -sin(lat)*cos(long)  -sin(lat)*sin(long)  cos(lat);               
              cos(lat)*cos(long)   cos(lat)*sin(long)  sin(lat)]; 
          
%**COMPUTE UNIT VECTORS            
vENU = R*dxyz ;         
        
end

%***EXTRACTING INDIVIDUAL LOCAL COORDINATES
E = vENU(1);% 'East'-coordinate relative to local origin (meters)
N = vENU(2);% 'North'-coordinate relative to local origin (meters)
U = vENU(3);% Up-coordinate relative to local origin (meters)

%***COMPUTING AZIMUTH(Az) AND ELEVATION ANGLE(El)FROM NEU 
hdist(i,1) = sqrt(E.^2+N.^2);%Horizontal Distance
if hdist(i,1) < 1.e-20
   Az_rad(i,1) = 0; % Radians
   El_rad(i,1) = pi/2;% Radians
else
    
 %***COMPUTE AZIMUTH(Az);   
try
   Az_rad(i,1) = atan2(E,N);% Radians
catch
   Az_rad(i,1)=rad2deg(atan2(E/norm(vENU),N/norm(vENU)));% Radians
end
if Az_rad(i,1) < 0
   Az_rad(i,1) = Az_rad(i,1)+ 2*pi;% Radians
end  
%****COMPUTE ELEVATION ANGLE(El)
try
   El_rad(i,1) = atan2(U,hdist(i,1));% Radians
catch
   El_rad(i,1)=asin(U/norm(vENU));% Radians
end

end

catch
%REMARK:
%       This an alternative way to compute satellite Azimuth and Elevation
%       Angle.Both(the try and catch) are suppose to produce same results.

%**COMPUTE UNIT VECTOR FROM OBSERVATION STATION(USER) TO SATELLITE POSITION*
r=sqrt((dxyz(1)).^2+(dxyz(2)).^2+(dxyz(3)).^2);
Dxyz=dxyz./r;
dx = Dxyz(1);
dy = Dxyz(2);
dz = Dxyz(3);

%****COMPUTE ROTATED UNIT VECTORS IN ENU FORM OBSERVATION STATION TO ...   *
%                                                        SATELLITE POSITION*
try
%*****COMPUTING Rotation matrix
%REMARKS:
%       R is direction cosine matrix to transform position vector components 
%       from geocentric equatorial frame into the topocentric horizon frame.    

R = [ -sin(long)                        cos(long)         0     ;  
      -sin(lat)*cos(long)  -sin(lat)*sin(long)  cos(lat);               
      cos(lat)*cos(long)   cos(lat)*sin(long)  sin(lat)]; 

%**COMPUTE UNIT VECTORS  
vENU = R*Dxyz;  

%***EXTRACTING INDIVIDUAL LOCAL COORDINATES
E = vENU(1);% 'East'-coordinate relative to local origin (meters)
N = vENU(2);% 'North'-coordinate relative to local origin (meters)
U = vENU(3);% Up-coordinate relative to local origin (meters)
   
catch
    
N =  sin(lat).*(-dx.*cos(long) - dy.*sin(long))+dz.*cos(lat);
E = - dx.* sin(long)+ dy.*cos(long) ;
U = cos(lat).*(dx.*cos(long)+dy.*sin(long)) + dz.*sin(lat);
               
end

%***COMPUTING AZIMUTH(Az) AND ELEVATION ANGLE(El)FROM NEU 
        
%****COMPUTE ELEVATION ANGLE(El)
%ElRAD=asin(U/norm(vENU));% Radians
try
	El_rad(i,1)= (pi / 2 - acos(U));% Radians
    
catch 
   LE = sqrt(E.^2+N.^2);
   El_rad(i,1) =(pi/2) - atan2(LE,U);% Radians
end
   El_deg(i,1)=rad2deg(El_rad(i,1));%Degrees
    
%***COMPUTE AZIMUTH(Az);
try
Az_rad(i,1) = atan2(E,N);% Radians
catch
Az_rad(i,1)=atan2(E/norm(vENU),N/norm(vENU));% Radians
end
%***CHECK FOR NEGATIVE ANGLES
if (Az_rad(i,1) < 0);
Az_rad(i,1) = Az_rad(i,1) + 2 * pi;
end
Az_deg(i,1) = rad2deg(Az_rad(i,1));	% degrees 
end    

%***COMPUTE RANGE
range(i,1) = sqrt(dxyz(1)^2+dxyz(2)^2+dxyz(3)^2);

%***CONVERT ANGLES IN RADIANS TO DEGREES
Az_deg(i,1) = rad2deg(Az_rad(i,1));%Degrees
El_deg(i,1) = rad2deg(El_rad(i,1));%Degrees
    
end

%***ARRANGE DATA FOR EACH USER POSITION/STATION
satAZ_rad(:,j)=Az_rad;  %Satellite Azimuth in radian
satEL_rad(:,j)=El_rad; %Satellite Elevation in radian
satAZ_deg(:,j)=Az_deg ; %Satellite Azimuth in decimal degrees
satEL_deg(:,j)=El_deg; %Satellite Elevation in decimal degrees
Range(:,j)=range;% Range(Distance)b/n User Position & Satellite

end %for j=1:size(UserXYZ,2)

%%%%%%%%%%%%%%%%%%%%%%%%END OF satAzimuthElevation.m %%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-