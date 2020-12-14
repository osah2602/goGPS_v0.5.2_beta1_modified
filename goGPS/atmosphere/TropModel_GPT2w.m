%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%*************"TropModel_GPT2w" is a routine that determines the           *
%              Meteorological Parameters and then  Computes the wet,dry,... *
%              and Total Tropospheric Delays Using Hydrostatic Delay Model *  
%              by Saatamoinen as refined by Davis et al(1985),  Askne ...  * 
%              and Nordius (1987)Wet Model and Vienna Mapping function1(VMF1).    
%                                             
%*************The routine determines pressure, temperature, temperature    *
%             lapse rate, mean temperature of the water vapor, water vapor * 
%             pressure, hydrostatic and wet mapping function coefficients  * 
%             ah and aw, water vapour decrease factor and geoid undulation * 
%             for specific sites near the Earth surface It is based on:    *                                             
%             1.a 1 x 1 degree external grid file ('gpt2_1w.grd')          *
%             2.a 5 x 5 degree external grid file ('gpt2_5w.grd') with mean*
%             values as well as sine and cosine amplitudes for the annual  * 
%             and semiannual variation of the coefficients.                *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
%NOTE:  This file goes together with "readGPTgrid.m", which reads the grid *
%       and saves it in cell structures.                                   *
%       +The "readGPTgrid.m" should be ran first to provide the grid input *
%       file (grid)                                                        *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
%                                                                    
%USAGE:
%      [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_GPT2w(Time,ReceiverPos,SatPos,grid,grid_res,Timevar)
%                                                         
%Others:                                                                   *
%      [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_GPT2w(Time,ReceiverPos,ElevationAngle,grid,grid_res,Timevar)
%                                                   
%      [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_GPT2w(Time,ReceiverPos,SatPos,grid,grid_res) 
%                                                         
%      [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_GPT2w(Time,ReceiverPos,[],grid,grid_res) 
%                                                           
%*******EXAMPLE:                                                           +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%1.EXTRACT GRID FILE USING "readGPTgrid.m"                                 +
%USAGE:                                                                    +
%      [gpt_grid] = readGPTgrid(GRIDfile,GPTmodel,grid_res)                +
%WHERE:                                                                    +
%a.    GRIDfile : GPT model grid file in single quote (e.g:'gpt2_1w.grd',  +
%                 'gpt2_5w.grd', or 'gpt2_1w.mat','gpt2_5w.mat')           +   
%b.    GPTmodel : GPT model type. indicate by 'GPT2w'                      +
%c.    grid_res : Grid resolution. 1 for 1x1 and 5 for 5x5 resolutions resp.

%e.g.: grid = readGPTgrid('gpt2_5w.grd','GPT2w',5) ;                       +
%           OR                                                             +
%      grid = readGPTgrid('gpt2_5w.mat','GPT2w',5) ;                       +

%2. NOW RUN THE "TropModel_GPT2w.m" TO COMPUTE ZENITH DELAYS               +
%INPUTS(eg):                                                               +
%1.Time        = [2019 1 5]                                                +
%2.ReceiverPos = [6 40 3.7449 -1 34  0.7263  260.2410]                     +
%3.SatPos      = [];                                                       +
%4.grid        = grid (i.e.extracted from "readGPTgrid.m")                 +
%5.grid_res    = 5                                                         +
%6.Timevar     = 0 i.e. 0 means GPT2w with time variation, 1 means static  +

%==> [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_GPT2w(Time,ReceiverPos,SatPos,grid,grid_res,Timevar)
%
%i.e.[STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_GPT2w([2019 1 5],[6 40 3.7449 -1 34  0.7263  260.2410],[],grid,5,0)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%**THE FUNCTION CALLS THE FOLLOWING SUBROUTINE(s):                         *
%1.[X,Y,Z] = geo2xyz(latitude,longitude,height,RefEllipsoid);              *
%2.[latRAD,longRAD,h,latDEC,longDEC] = xyz2LLH(X,Y,Z,RefEllipsoid);        *
%3.[a,finv] = Elipsoidpara(RefEllipsoid);                                  *
%4.[Az_rad,El_rad,Az_deg,El_deg,D]=satAzimuthElevation(UserXYZ,SatXYZ,...  *
%                                                            RefEllipsoid);*
%5.[P,T,dT,Tm,es,ah,aw,lambda,undu] = GPT2w_5(UTCtime,lat,long,h,Timevar); *
%6.[FileName, FilePath] = searchGRID(varargin)                             *
%7.[VMFh,VMFw] = VMF(ah,aw,UTCtime,lat,hgt,satELEV)                        *
%8.ZHD =ZenithTropDelay(P,es,Tm,lambda )                                   *
%9.SHD =SlantTropDelay(P,es,Tm,lambda,MFh,MFw)                             *
%10.[JD, MJD,DoY]=utc2JulianDay_DoY(Year,Month,Day,Hour,Minute,Seconds)    *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%      Generally, the function accepts six(6) sets of inputs:              *
%1.           Time: Observation time in[year,Month,Day,Hour,Min,Sec]       *
%2.    ReceiverPos: Receiver position in either Latitude,Longitude &...    *
%                     height or ECEF(XYZ) Coordinates                      *
%3.         SatPos: Satellite position(s) in ECEF(XYZ) Coordinates         *
%4            grid: grid values in cells extracted from the grid file      *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%NOTE: To get the grid values(grid) the "readGPTgrid.m" subroutine has to  *
%      be ran first.
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%5.       grid_res: grid resolution (°) (possible: 1 or 5).indicate 1 for  *
%               1°x1° external grid file and 5 for 5°x5° external grid file*

%6.        Timevar: is an indicator of time variation to either compute Met*
%                   parameters in staic mode or in annual and semiannual   * 
%                   terms. Indicate by specifying 1 or 0.                  *
%              case 1: no time variation but static quantities             *
%              case 0: with time variation (annual and semiannual terms)   *

%OTHER CONSIDERATIONs:                                                     *
%1.*******For only zenith tropospheric delays, set SatPos to empty([])     *

%2.*******Satellites Elevation Angle(satEL) can also be used in place of...*
%         Satellites Position (SatPos).  format should be (nx1)or (n x 3)  * 
%         NOTE:                                                            *
%              Elevation Angles should be decimal degrees(eg:26.981.78.102)*
%3.*******Where SatPos is empty ([]) or Not a number(nan), the subroutine  *          
%         Computes the delays only in the zenith direction                 * 
%4.*******Where only Ellipsoidal Height is provided as input for receiver  *      
%         position(ReceiverPos), the subroutine returns an error message...*        
%5.*******Again, for empty ([]) or not a number (nan) entry/input for ...  *
%         ReceiverPos and Time, the subroutine returns an error message    *
%6.*******For only one(1) input, error message is given and the program    *
%         Terminates / stop working given empty matrices([]) or zeros as   * 
%         outputs                                                          *                                                         
%NOTE:                                                                     *
%1.   a)Elevation Angles should be decimal degrees(eg:26.981,78.102)       *
%     b)HeightS are also given in meters                                   *
%2.   Receiver positions are in the ff formats:                            *
%     a)For Geographic Coordinates (n x m):[Latitude Longitude height(h)]  *
%      eg:DMS:[6 40 21 -1 33 17 187.76]  OR                                *
%          DM:[6 40.35 -1 33.2833 187.76] OR                               *
%         Deg:[6.6725 -1.5547 187.76] Decimal degrees                      *

%For Western Longitudes place negative(-) infront of the value             *
% eg:[-1 33 34] or [-001 33 34] or [ 000 -24 16] or [ -1 33.2]OR [-1.3345] *

%For Southern Latitudes place negative(-) infront of the nunmber           *
%  eg:[-8 23 14] or [-008 23 14] or [ 000 -23 26] or[ -8 23.6]OR [-8.2315] *

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
%                --------------------------------------------------        *                             
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%OUTPUT:                                                                   *
%1.     STD => Slant Total Tropospheric Delay in meters                    *
%2.     SHD => Slant Hydrostaic Tropospheric Delay in meters               *
%3.     SWD => slant Wet Tropospheric Delay  in meters                     *
%4.     ZTD => Zenith Total Tropospheric Delay in meters                   *
%5.     ZHD => Zenith Hydrostaic Tropospheric Delay in meters              *
%6.     ZWD => Zenith Wet Tropospheric Delay  in meters                    *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================
%%REFERENCE:                                                               *
%           J. Böhm, G. Möller, M. Schindelegger, G. Pain, R. Weber,...    *
%           Development of an improved blind model for slant delays in the *
%           troposphere (GPT2w),GPS Solutions, 2014,...                    *
%           doi:10.1007/s10291-014-0403-7                                  *

%Original codes by Böhm et al 2014                                         *
%Modified by: Samuel Osah,Msc Geomatic Engineering ,2016                   *    
%             Email: osahsamuel@yahoo.ca                                   *
%             Tel:+233 (0)246137410/+233 (0)509438484                      *  
%==========================================================================+
%**************************************************************************+
%**************************************************************************+

function [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_GPT2w(Time,ReceiverPos,SatPos,grid,grid_res,Timevar)
                                                                                                        
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%**********CHECK INPUTs & REFORMAT INPUT DATA
switch nargin
    
    case {6,5,4,3,2} %When all inputs are provided        
        
        if (any(nargin==[6,5,4,3,2]))
            
          if nargin ==6 %for six(6) or all inputs...  
             
             %******Assignments 
             if isempty(grid)%grid data input is empty                
              %Retrieve stored grid 
              grid =getappdata(0,'gpt_grid' ) ;
            end
                        
            if any([isempty(grid_res),isnan(grid_res)])%grid resolution input is empty                
               grid_res = 5;%set grid resolution to 5
            end
            
            %******Assignments  
            if any([isempty(Timevar),isnan(Timevar)])%Timevar input is empty                
               Timevar = 0;%set time variation incator to 0
            end
            
          elseif nargin == 5 %for ONLY five(5) inputs...
              
                 %******Assignments 
                 
                 if isempty(grid)%grid data input is empty                
                    %Retrieve stored grid 
                    grid =getappdata(0,'gpt_grid' ) ;
                 end 
            
                 if any([isempty(grid_res),isnan(grid_res)])%grid resolution input is empty                
                    grid_res = 5;%set grid resolution to 5
                 end 
                 
                  Timevar = 0;%set time variation incator to 0
            
          elseif nargin ==4 %for ONLY four(4) inputs... 
              
                 %******Assignments
                 
                 if isempty(grid)%grid data input is empty                
                    %Retrieve stored grid 
                    grid =getappdata(0,'gpt_grid' ) ;
                 end 
                 
                  grid_res = 5;%set grid resolution to 5
                  Timevar = 0;%set time variation incator to 0
                 
          elseif  nargin ==3 %for ONLY three(3) inputs... 
              
                  %******Assignments                            
                  grid =getappdata(0,'gpt_grid' ) ; %Retrieve stored grid 
                  grid_res = 5;%set grid resolution to 5
                  Timevar = 0;%set time variation incator to 0
                 
          elseif nargin == 2
              
                 %******Assignments
                 SatPos  = [];%Assign empty matrix([]) to SatPos
                 grid =getappdata(0,'gpt_grid' ) ; %Retrieve stored grid 
                 grid_res = 5;%set grid resolution to 5
                 Timevar = 0;%set time variation incator to 0
                   
          end %if nargin ==6
          
          %*******************CHECK FOR EMPTY([]) INPUTs
          
         if isempty(Time)
            
           %ISSUE ERROR MESSAGE  
           beep%Give a beep sound 
           errmsg0{1}=sprintf('Observation / Reception Time is not provided i.e it''s Empty ([ ]).\n');
           errmsg0{2}='Please provide Observation / reception Time & Try again.';
           errordlg(errmsg0,'Time Input Error','modal')  
             
           %RETURN EMPTY([]) DELAYS
           STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           return
        end
        
        if isempty(ReceiverPos)
            
           %ISSUE ERROR MESSAGE 
            beep%Give a beep sound 
            errmsg10{1}=sprintf('Reciever / Site Position is not provided i.e it''s Empty ([ ]).\n');
            errmsg10{2}='Please provide Receiver / Site position & Try again.';
            errordlg(errmsg10,'Reciever Position(s) Input Error','modal')  
           
            %RETURN EMPTY([]) DELAYS
            STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
            return
            
        elseif ~isempty(ReceiverPos)
               %IF RECEIVER POSITION IS NOT EMPTY,BUT ONLY HEIGHTS ARE PROVIDED
               
               if size(ReceiverPos,2)==1 
                  
                   %ISSUE ERROR MESSAGE for HEIGHT INPUT 
                    beep%Give a beep sound 
                    errmsg50{1}=sprintf('Insuficient Input for Receiver / Station Position .\n');
                    errmsg50{2}=sprintf('Latitude & Longitude Coordinates have not been provided.\n');
                    errmsg50{3}=sprintf('Please provide Latitude & Longitude Coordinates & Try again.');
                    errordlg(errmsg50,'Coordinate(s) Input Error','modal')  
                    
                    %RETURN EMPTY([]) PARAMETERs
                    STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
            
                    return
                    
               else 
                    Rpos=ReceiverPos;%Receiver position(XYZ/LAT LONG h) 
               end
               
        end
          
          
       if isempty(grid)
                              
          %ISSUE ERROR MESSAGE  
          beep%Give a beep sound 
          errmsg21{1}=sprintf('No input for VMF grid file, i.e. input is empty ( [ ] ).\n');
          errmsg21{2}='Please provide VMF grid file & Try again.';
          errordlg(errmsg21,'Grid file Input Error','modal') 
                     
          %RETURN EMPTY([]) DELAYS
          STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           
          return  
          
       end
                         
       if isempty(grid_res)
                      
          %ISSUE ERROR MESSAGE  
          beep%Give a beep sound 
          errmsg20{1}=sprintf('No input for VMF grid resolution, i.e. input is empty ( [ ] ) .\n');
          errmsg20{2}='Please provide VMF grid resolution & Try again.';
          errordlg(errmsg20,'Grid resolution Input Error','modal') 
                     
          %RETURN EMPTY([]) DELAYS
          STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           
          return 
                         
       end %//if isempty(grid_res)
            
          %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-          
          
          %*****CHECK IF SATELLITE ELEVATION ANGLE IS ENTERED INSTEAD
          %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
                                                                                          
          %CHECK IF INPUTS(RECEIVER & SATELLITE POSITIONs ARE NOT EMPTY([])
          %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          %1.SATELLITE POSITION [SatPos]
           if (~isempty(SatPos)| ~isnan(SatPos))
              
             %***ROUND TO THE NEAREST THOUSAND             
             t2=roundn(SatPos,3);
                  
             if  size(SatPos,2)==1 %if the number of columns in SatPos is 1                                                      
               if  all(t2(:,1))==0                       
                 satEL_deg=SatPos;%Assigning SatPos to satEL(Satellite Elevation) 
                 satEL_rad=satEL_deg.* pi / 180;%Elevation Angles converted to radian
               end   %if all(t2(:,1))==0                  
             
             else 
                 if all(t2(:,:))==0 %if all rows and columns rounded to zeros(0's)
                                     %then inputs are elevation angles
                    satEL_deg=SatPos;%Assigning SatPos to satEL(Satellite Elevation) 
                    satEL_rad=satEL_deg.* pi / 180;%Elevation Angles converted to radian 
                    
                 else                          
                     Satpos=SatPos;%Satellite position in XYZ                      
                 end
                 
             end   %\\if size(SatPos,2)==1
             
                                    
           end %\\if (~isempty(SatPos)| ~isnan(SatPos))
                                    
          %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          %******CHECK FOR TIME INPUTs
          if any(~isempty(Time(:,:)))| any(~isnan(Time(:,:)))
            
            %CHECK THE TIME FORMAT
            ncol=size(Time,2);% finding Number of columns 
            
            switch ncol              
                case {6,5,4,3}
                  
                  if (any(ncol==[6,5,4,3]))
                    
                     Date=Time(:,1:3);%Assigning column 1 to 3 of time to date                     
                     %***Date
                     Year=Date(:,1);%Assigning 1st column to years
                     Month=Date(:,2);%Assigning 2ND column to months
                     Day=Date(:,3);%Assigning 3RD column to days
                     
                    if ncol==6                                                 
                      time=Time(:,4:6);%Assigning column 4 to 6 of tyme to time                                            
                      %***Time
                      Hour=time(:,1);%Assigning 1st column to hours    
                      Minute=time(:,2);%Assigning 2nd column to min  
                      Seconds=time(:,3);%Assigning 3rd column to sec                       
                    elseif ncol==5                         
                          time=Time(:,4:5);%Assigning column 4 to 5 of to tyme                          
                          %Time
                          Hour=time(:,1);%Assigning 1st column to hours    
                          Minute=time(:,2);%Assigning 2nd column to min               
                          Seconds=zeros(size(Hour,1),1);%Assigning zeros to Seconds                 
                    elseif  ncol==4 
                          time=Time(:,4);%Assigning column 4 of time to tyme
                          Hour=time(:,1);%Assigning 1st column to hours 
                          Minute=zeros(size(Hour,1),1);%Assigning zeros to Minute
                          Seconds=zeros(size(Hour,1),1);%Assigning zeros to Seconds                           
                    elseif  ncol==3
                          Hour=zeros(size(Year,1),1);%Assigning zeros to Hour
                          Minute=zeros(size(Year,1),1);%Assigning zeros to Minute
                          Seconds=zeros(size(Year,1),1);%Assigning zeros to Seconds
                                            
                    end   %if ncol==6
                    
                    %*****CHANGE SECs, MIN, HOUR WHOSE SEC == 60
                     Minute(Seconds==60) = Minute(Seconds==60)+1;
                     Seconds(Seconds==60) = 0;
                     Hour(Minute==60) = Hour(Minute==60)+1;
                     Minute(Minute==60)=0;
                     
                    UTCtime=[Year Month Day Hour Minute Seconds];%UTC time
                    
                  end   %if (any(ncol==[6,5,4,3])) 
                  
                     %******CHECK THAT YEAR,MONTH & DAY ARE INTEGERS   
                     if any(any( (Date-fix(Date)) ~= 0 ))
             
                       %****DISPLAY ERROR MESSAGE          
                       beep %Give a beep sound 
                       errmsg{1}=sprintf('Year,Month & Day must be integers \n');                       
                       errmsg{3}=sprintf('Please repeat with valid values/Check Date. \n');
                       errordlg(errmsg,'Date Input Error','modal')            
                       return;
                            
                       %******CHECK THE VALIDITY OF MONTH, DAY & TIME 
                     elseif  (any(Month(:)>12  | Month(:)<1 ))
             
                           %****DISPLAY ERROR MESSAGE          
                           beep %Give a beep sound  
                           errmsg{1}=sprintf('The Month input must be a number in the range [1-12] \n');                           
                           errmsg{2}=sprintf('Please repeat with valid values/Check Date. \n');
                           errordlg(errmsg,'Date(Month) Input Error','modal')
                           return; 
               
                     elseif  (any(Day(:)>31 | Day(:)<1 ))
             
                           %****DISPLAY ERROR MESSAGE          
                           beep %Give a beep sound  
                           errmsg{1}=sprintf('The Day input must be a number in the range [1-31] \n');                          
                           errmsg{2}=sprintf('Please repeat with valid values/Check Date. \n');
                           errordlg(errmsg,'Date(Day) Input Error','modal') 
                           return
               
                           %*****************CHECK TIME VALIDITY
                     elseif  (any(Hour(:)>24 | Hour(:)<0))
             
                           %****DISPLAY ERROR MESSAGE          
                           beep %Give a beep sound  
                           errmsg{1}=sprintf('The Hour time input must be a number in the range [0-24] \n');                        
                           errmsg{2}=sprintf('Please repeat with valid values/Check Time. \n');
                           errordlg(errmsg,'Time(Hour) Input Error','modal')
                           return;
               
                     elseif  (any(Minute(:)>60 | Minute(:)<0))
             
                           %****DISPLAY ERROR MESSAGE          
                           beep %Give a beep sound  
                           errmsg{1}=sprintf('The Minute time input must be a number in the range [0-60] \n');                          
                           errmsg{2}=sprintf('Please repeat with valid values/Check Time. \n');
                           errordlg(errmsg,'Time(Minute) Input Error','modal')
                           return;
               
                     elseif  (any(Seconds(:)>60 | Seconds(:)<0))
             
                           %****DISPLAY ERROR MESSAGE          
                           beep %Give a beep sound  
                           errmsg{1}=sprintf('The Seconds time input must be a number in the range [0-60] \n');                           
                           errmsg{2}=sprintf('Please repeat with valid values/Check Time. \n');
                           errordlg(errmsg,'Time(Seconds) Input Error','modal')
                           return;
                          
                     end  %if any(any( (date-fix(date)) ~= 0 ))                       
                               
            end   %switch ncol
            
          else %if time is not provided
               %****DISPLAY ERROR MESSAGE          
               beep %Give a beep sound  
               errmsg{1}=sprintf('Time input is empty ( [ ] ) \n');              
               errmsg{2}=sprintf('Please Provide Time input and try again. \n');
               errordlg(errmsg,'Time Input Error','modal')
                                    
          end %if any([~isempty(Time),~isnan(Time)])
          
        end  %if (any(nargin==[4,3,2]))
         
        %*******IF USER INPUTs are Coordinates,.....
        if exist('Rpos','var')
          
           %CHECK RECEIVER POSITION TYPE(XYZ/LAT LONG h)
           [Rrow,Rcol]=size(Rpos);%Get number of rows & columns of Rpos
          
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
             latD=dms2degrees(lat);%Convert Latitude in DMS to degrees
             lonD=dms2degrees(lon);%Convert Longitude in DMS to degrees
             
          elseif size(lat,2)== 2
                 latD=dm2degrees(lat);%Convert Latitude in DM to degrees
                 lonD=dm2degrees(lon);%Convert Longitude in DM to degrees 
                 
          else
               latD = lat;%Latitude in [degrees]
               lonD = lon;% Longitude in [degrees] 
                 
          end
          
          latR=latD.*degrad; %latitude in radian
          lonR=lonD.*degrad; %longitude in radian
          
          %SAVE COORDs for use
          setappdata(0,'lat',latR)
          setappdata(0,'lon',lonR)
          setappdata(0,'ht',h)
                    
          %**CONVERT USER POSITION(XYZ) TO lat,long,h
        elseif (exist('X','var')|| exist('Y','var') || exist('Z','var'))
               
              %Call the 'xyz2LLH.m' function
              [latR,lonR,h,latD] = xyz2LLH(X,Y,Z); 
              
              %SAVE COORDs for use
              setappdata(0,'lat',latR)
              setappdata(0,'lon',lonR)
              setappdata(0,'ht',h)
                  
        end %exist('lat','var') || exist('lon','var')
        
          %*******COMPUTE SATELLITE ELEVATION ANGLE IF NOT PROVIDED                   
        if ~exist('satEL_deg','var')%If Satellite Elevation is not given,then 
                                    %Compute Satellite Elevation with the given  
                                    %Receiver & Satellite positions
           
           if exist('Satpos','var')                        
         
             %Call the "satAzEl.m" function 
             [satEL_deg,satEL_rad]=satAzEl([X Y Z],Satpos);%compute satellite elevation angle
             
           end
         
        end %if ~exist('satEL','var'
        
        end %if exist('Rpos','var')
                                                                       
    otherwise
             %ISSUE ERROR MESSAGE for 4,3,2,1 INPUTs
              beep%Give a beep sound 
              errmsg{1}=sprintf('Insuficient Data Input / Wrong Data Input format .');
              errmsg{2}='';
              errmsg{3}='Please Check inputs / Data format & Try Again.';
              errordlg(errmsg,'Coordinate(s) Input Error','modal')  
              
              %Return empty ([]) outputs
              STD=[]; SHD=[]; SWD=[]; ZTD=[]; ZHD=[]; ZWD=[];
              
              return                          
end %switch nargin

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%************TROPOHERIC DELAY MODELING USING GPT2w 5x5 degree MODEL
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
%1.********COMPUTE ZENITH TROPOSPHERIC DELAYs
%Call the "getZTD_GPT2w.m" Function
[ZHD, ZWD, ZTD,ah,aw,T] = getZTD_GPT2w(UTCtime,latR,lonR,h,grid,grid_res,Timevar); %#ok<*ASGLU>
 
%2.**********COMPUTE SLANT TROPOSPHERIC DELAYs
if exist('satEL_deg','var')
    
   %(2.1)****COMPUTE MAPPING FUNCTIONs(MF)BY VMF1
   %Call the "VMF.m" function
   [VMFh,VMFw]=VMF1(ah,aw,UTCtime,latR,h,satEL_deg);   
     
   %(2.2)Call the "SlantTropDelay.m" function
   [STD, SHD, SWD] =SlantTropDelay(ZHD,ZWD,VMFh,VMFw);
else
     %Create zero matrix
     STD=zeros(size(ZHD));
     SHD=zeros(size(ZHD));
     SWD=zeros(size(ZHD));
end
 
%SAVE TEMPERATURE FOR USE
setappdata(0,'T_GPT',T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TropModel_GPT2w_5_x_5.m 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%A.SUBROUTINE TO COMPUTE SLANT TROPOSPHERIC DELAY (SATELLITE IN VIEW)
function [STD, SHD, SWD] =SlantTropDelay(ZHD,ZWD,MFh,MFw)                                                 
                                             
%DESCRIPTION:
%            "SlantTropDelay" Computes Slant Tropospheric Delay Using ...  * 
%             Davis et al(1985), and Askne and Nordius (1987) models for   *
%             Zenith Hydrostaic and Wet Tropospheric Delays Compuation and *
%             Vienna  Mapping function 1 (VMF1) for the Slant  delays.     *            
%******SYNTAX:                                                             *
%             SlantTropDelay(ZHD,ZWD,MFh,MFw)                              *
%******INPUT:                                                              *
%           ZHD = Zenith Hydrostatic Delay computed from  Davis et al model*
%           ZWD = Zenith Wet Delay computed from Askne and Nordius model   *   
%           MFh = Hydrostatic Mapping Function computed from VMF1 model    *
%           MFw = Wet Mapping Function computed from VMF1 model            *
   
%****OUTPUT:
%           SHD = Slant Hydrostatic Delay                                  *
%           SWD = Slant Wet Delay                                          *
%           STD = Slant Total Delay                                        *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%*********COMPUTE SLANT TROPOSPHERIC DELAY (SATELLITE IN VIEW)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%(1)*****INITIALIZING OUTPUT VARIABLEs

if any([iscell(MFh),iscell(MFw)])%if either of MFh,MFw is a cell array
    
  [nrow,ncol,dim] = size(MFh);%Get size of cell array
  SHD = cell(nrow,ncol,dim);
  [STD,SWD] = deal(SHD);%copy the contents of SHD to all the requested outputs
  
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
    
    [STD,SWD] = deal(SHD);%copy the contents of SHD to all the requested outputs
    
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
function [ZHD,ZWD,ZTD,ah,aw,T] = getZTD_GPT2w(UTCtime,lat,lon,hell,grid,grid_res,Timevar)

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
      ZHD=zeros(nrow_time,1);
      
      [ZWD,ZTD,ah,aw,T]=deal(ZHD);%create copy of ZHD in ZWD,ZTD,ah,aw,T
   
      for i=1:nrow_time
       
          if grid_res == 1 %if grid resolution is 1°x 1°
              
             %Call the "GPT2w_1_x_1.m" Function
             [ZHD(i,1),ZWD(i,1),ZTD(i,1),ah(i,1),aw(i,1),T(i,1)] = GPT2w_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid,Timevar); 
          
          elseif grid_res == 5 %if grid resolution is 5°x 5°
              
                 %Call the "GPT2w_5_x_5.m" Function
                [ZHD(i,1),ZWD(i,1),ZTD(i,1),ah(i,1),aw(i,1),T(i,1)] = GPT2w_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid,Timevar);
          end 
      end
      
   else 
       %*****INITIALIZE OUTPUTs 
       ZHD=zeros(nrow_pos,nrow_time);
      
      [ZWD,ZTD,ah,aw,T]=deal(ZHD);%create copy of ZHD in ZWD,ZTD,ah,aw,T
    
     for i=1:nrow_time %LOOP OVER TIME
         
        for j=1:nrow_pos %LOOP OVER POSITIONS
            
            %****CHECK GRID RESOLUTION
            
            if grid_res == 1 %if grid resolution is 1°x 1°
                
               %Call the the 1 degree grid version; "GPT2w_1_x_1.m"  
               [ZHD(j,i),ZWD(j,i),ZTD(j,i),ah(j,i),aw(j,i),T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);
           
            elseif grid_res == 5 %if grid resolution is 5°x 5°
                
                   %Call the the 5 degree grid version; "GPT2w_5_x_5.m"  
                   [ZHD(j,i),ZWD(j,i),ZTD(j,i),ah(j,i),aw(j,i),T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);   
            end
        end
     end
    
   end %//if isequal(Identical,1)
    
else
    %*****INITIALIZE OUTPUTs 
    ZHD=zeros(nrow_pos,nrow_time);
    
    [ZWD,ZTD,ah,aw,T]=deal(ZHD);%create copy of ZHD in ZWD,ZTD,ah,aw,T
    
    for i=1:nrow_time %LOOP OVER TIME
        
        for j=1:nrow_pos %LOOP OVER POSITIONS
            
            %****CHECK GRID RESOLUTION
            
            if grid_res == 1 %if grid resolution is 1°x 1°
                
               %Call the the 1 degree grid version; "GPT2w_1_x_1.m"  
               [ZHD(j,i),ZWD(j,i),ZTD(j,i),ah(j,i),aw(j,i),T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);
           
            elseif grid_res == 5 %if grid resolution is 5°x 5°
                
                   %Call the the 5 degree grid version; "GPT2w_5_x_5.m"  
                   [ZHD(j,i),ZWD(j,i),ZTD(j,i),ah(j,i),aw(j,i),T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);   
            end
        end
    end
end
%***********************END OF getZTD_GPT2w.m ***********************    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%B_1.SUBROUTINE TO COMPUTE METEOROLOGICAL PARAMETERS FOR ZTDs FROM 5°x 5° GRID MODEL
function [ZHD,ZWD,ZTD,ah,aw,T] = GPT2w_5_x_5(UTCtime,lat,lon,hell,grid,Timevar)

%**************************************************************************
%DESCRIPTION:
%            This subroutine determines pressure, temperature, temperature *
%            lapse rate, mean temperature of the water vapor, water vapor  * 
%            pressure, hydrostatic and wet mapping function coefficients   * 
%            ah and aw, water vapour decrease factor and geoid undulation  * 
%            AND FINALLY computes Zenith tropospheric delays for specific  *
%            sites near the Earth surface. It is based on a  5°x 5°        *
%            external grid file ('gpt2_5w.grd') with mean values as well as*
%            sine and cosine amplitudes for the annual and semiannual      * 
%            variation of the coefficients.                                *
%USAGE:
%      [ZHD,ZWD,ZTD,ah,aw] = GPT2w_5(UTCtime,lat,lon,hell,Timevar)         *
%INPUTs:                                                                   *
%1.   UTCtime:  UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *                                                                
%2.       lat:  Station geodetic latitude in radians [-pi/2:+pi/2] (n x 1) *
%3.       lon:  Station geodetic longitude in radians [-pi:pi] or [0:2pi]  *
%4.      hell:  ellipsoidal height in [meters] (n x 1)                     *
%5    Timevar:  case 1: no time variation but static quantities            *
%               case 0: with time variation (annual and semiannual terms)  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%OUTPUTs:
%1.     ZHD: Zenith Hydrostatic Delay                                      *
%2.     ZWD: Zenith Wet Delay                                              *
%3.     ZTD: Zenith Total Delay                                            *
%4.      ah: hydrostatic mapping function coefficient at zero height (VMF1)* 
%5.      aw: wet mapping function coefficient (VMF1) (vector of length nstat)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%NOTE:                                                                     *     
%    The hydrostatic mapping function coefficients have to be used with the*
%     height dependent Vienna Mapping Function 1  because the              *
%     coefficients refer to zero height.                                   *
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

%--------------------------------------------------------------------------
%(1)***********************EXTRACT METEOROLOGICAL PARAMETERS
%--------------------------------------------------------------------------

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
p    =  zeros([nstat, 1]);
T    =  zeros([nstat, 1]);
dT   = zeros([nstat, 1]);
Tm   = zeros([nstat, 1]);
e    =  zeros([nstat, 1]);
ah   = zeros([nstat, 1]);
aw   = zeros([nstat, 1]);
la   = zeros([nstat, 1]);
undu = zeros([nstat, 1]);

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
        undu(k) = u_grid(ix);
        hgt = hell(k)-undu(k);
            
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
        redh = hgt - Hs_grid(ix);

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
            hgt = hell(k)-undul(l);
        
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
            Hs1 = Hs_grid(indx(l));
            redh = hgt - Hs1;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF MET PARAMETERS EXTRACTION

%--------------------------------------------------------------------------
%(2)***COMPUTE ZENITH TROPOSPHERIC DELAYS USING EXTRACTED MET PARAMETERS
%--------------------------------------------------------------------------
%NOTE:
%      DAVIS ET AL DRY MODEL REQUIRES +VE ORTOMETRIC HEIGHT TO PERFORM
%--------------------------------------------------------------------------

%COMPUTE CORRECTION FACTOR FOR LOCAL GRAVITY ACCELERATION 
   
dgref=1-0.00266.*cos(2.*lat)-0.00028e-3.*(hell-undu);

gm   = 9.784 .* dgref;%Acceleration Gravity at the Atmospheric column in m/s^2

%**************COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
%              =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%USE OF DAVIS ET AL (1985) DRY MODEL   

 ZHD=(0.0022768 .* p)./dgref;  %Zenith Hydrostatic or dry delay [m]    

%**************COMPUTE ZENITH WET DELAY(ZHD)
%              =-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%USE OF  Askne and Nordius (1987) WET MODEL      
%*****DEFINE CONSTANTs
k1 = 77.604; 				   % K/hPa
k2 = 64.79; 				   % K/hPa
Mw = 18.0152;%Molar Mass of Water in [g/mol]
Md = 28.9644; %Molar Mass of Dry Air  in [g/mol]
k2p= k2 - k1*Mw/Md;          % K/hPa  [called k2 prime]
k3 = 377600; 				   % KK/hPa          
dMtr = 28.965*10^-3;%molar mass of dry air in kg/mol    
R = 8.3143;%universal gas constant in J/K/mol         
Rd = R/dMtr ; %specific gas constant for dry air

%ASSIGNMENT
lambda = la;

%*****COMPUTE ZENITH WET DELAY(ZWD)
ZWD = ((1e-6.*(k2p + k3./Tm).*Rd)./(gm.*(lambda + 1))).*e;%[m]
                
%****FINALLY, COMPUTE ZENITH TOTAL DELAY(ZTD)
ZTD = ZHD + ZWD;%Zenith Total Delay [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GPT2w_5 x 5.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%B_2.******SUBROUTINE TO COMPUTE METEOROLOGICAL PARAMETERS FOR ZTDs FROM 1°x 1° GRID MODEL
function [ZHD,ZWD,ZTD,ah,aw,T] = GPT2w_1_x_1(UTCtime,lat,lon,hell,grid,Timevar)

%**************************************************************************
%DESCRIPTION:
%            This subroutine determines pressure, temperature, temperature *
%            lapse rate, mean temperature of the water vapor, water vapor  * 
%            pressure, hydrostatic and wet mapping function coefficients   * 
%            ah and aw, water vapour decrease factor and geoid undulation  * 
%            AND FINALLY computes Zenith tropospheric delays for specific  *
%            sites near the Earth surface. It is based on a  1°x 1°        *
%            external grid file ('gpt2_1w.grd') with mean values as well as*
%            sine and cosine amplitudes for the annual and semiannual      * 
%            variation of the coefficients.                                *
%USAGE:
%      [ZHD,ZWD,ZTD,ah,aw] = GPT2w_5(UTCtime,lat,lon,hell,Timevar)         *
%INPUTs:                                                                   *
%1.   UTCtime:  UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *                                                                
%2.       lat:  Station geodetic latitude in radians [-pi/2:+pi/2] (n x 1) *
%3.       lon:  Station geodetic longitude in radians [-pi:pi] or [0:2pi]  *
%4.      hell:  ellipsoidal height in [meters] (n x 1)                     *
%5    Timevar:  case 1: no time variation but static quantities            *
%               case 0: with time variation (annual and semiannual terms)  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%OUTPUTs:
%1.     ZHD: Zenith Hydrostatic Delay                                      *
%2.     ZWD: Zenith Wet Delay                                              *
%3.     ZTD: Zenith Total Delay                                            *
%4.      ah: hydrostatic mapping function coefficient at zero height (VMF1)* 
%5.      aw: wet mapping function coefficient (VMF1) (vector of length nstat)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%NOTE:                                                                     *     
%    The hydrostatic mapping function coefficients have to be used with the*
%     height dependent Vienna Mapping Function 1  because the              *
%     coefficients refer to zero height.                                   *
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
  Yr = UTCtime(:,1);%get Hour
  Mn = UTCtime(:,2);%get Month
 Day = UTCtime(:,3);%get Day
%2.UTCtime
   H = UTCtime(:,4);%get Hour
 MIN = UTCtime(:,5);%get Minute
SECs = UTCtime(:,6);%get Seconds

%--------------------------------------------------------------------------
%(1)***********************EXTRACT METEOROLOGICAL PARAMETERS
%--------------------------------------------------------------------------

%****MODIFIED JULIAN DAT
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
p    =  zeros([nstat, 1]);
T    =  zeros([nstat, 1]);
dT   = zeros([nstat, 1]);
Tm   = zeros([nstat, 1]);
e    =  zeros([nstat, 1]);
ah   = zeros([nstat, 1]);
aw   = zeros([nstat, 1]);
la   = zeros([nstat, 1]);
undu = zeros([nstat, 1]);

%****LOOP OVER STATIONs
for k = 1:nstat
    
    %Only +VE longitude in [degrees]
    if lon(k) < 0
        plon = (lon(k) + 2*pi)*180/pi;
    else
        plon = lon(k)*180/pi;
    end
    %Transform to polar distance in [degrees]
    ppod = (-lat(k) + pi/2)*180/pi; 

    %find the index (line in the grid file) of the nearest point
	%changed for the 1 degree grid (GP)
    ipod = floor((ppod+1)); 
    ilon = floor((plon+1));
    
    %Normalized (to one) differences, can be +VE or -VE
	%changed for the 1 degree grid (GP)
    diffpod = (ppod - (ipod - 0.5));
    difflon = (plon - (ilon - 0.5));
    
    %Added by HCY
	%Changed for the 1 degree grid (GP)
    if ipod == 181
        ipod = 180;
    end
	%Added by GP
    if ilon == 361
		ilon = 1;
    end
    if ilon == 0
		ilon = 360;
    end 

    %Get the number of the corresponding line
	%changed for the 1 degree grid (GP)
    indx(1) = (ipod - 1)*360 + ilon;
    
    %Near the poles: nearest neighbour interpolation, otherwise: bilinear
	%with the 1 degree grid the limits are lower and upper (GP)
    bilinear = 0;
    if ppod > 0.5 && ppod < 179.5 
           bilinear = 1;          
    end          
    
    %Case of nearest neighborhood
    if bilinear == 0

        ix = indx(1);
        
        %Transforming ellipsoidal height to orthometric height
        undu(k) = u_grid(ix);
        hgt = hell(k)-undu(k);
            
        %Pressure, Temperature at the height of the grid
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
        redh = hgt - Hs_grid(ix);

        %Temperature at station height in [Celsius]
        T(k) = T0 + dT(k)*redh - 273.15;
        
        %Temperature lapse rate in [degrees / km]
        dT(k) = dT(k)*1000;

        %virtual temperature in [Kelvin]
        Tv = T0*(1+0.6077*Q);       
        c = gm*dMtr/(Rg*Tv);
        
        %Pressure in [hPa]
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

		%Mean temperature of the water vapor Tm - added by GP
        Tm(k) = Tm_grid(ix,1) + ...
                Tm_grid(ix,2)*cosfy + Tm_grid(ix,3)*sinfy + ...
                Tm_grid(ix,4)*coshy + Tm_grid(ix,5)*sinhy;
		
		%Water vapor pressure in hPa - changed by GP
		e0 = Q*p0/(0.622+0.378*Q)/100; % on the grid
		e(k) = e0*(100*p(k)/p0)^(la(k)+1); % on the station height - (14) Askne and Nordius, 1987
		                    
     else % bilinear interpolation
        
        ipod1 = ipod + sign(diffpod);
        ilon1 = ilon + sign(difflon);
        
		% changed for the 1 degree grid (GP)
        if ilon1 == 361
            ilon1 = 1;
        end
        if ilon1 == 0
            ilon1 = 360;
        end
        
        %Get the number of the line
		%changed for the 1 degree grid (GP)
        indx(2) = (ipod1 - 1)*360 + ilon;  % along same longitude
        indx(3) = (ipod  - 1)*360 + ilon1; % along same polar distance
        indx(4) = (ipod1 - 1)*360 + ilon1; % diagonal
        
        for l = 1:4
              
            %Transforming ellipsoidal height to orthometric height:
            %Hortho = -N + Hell
            undul(l) = u_grid(indx(l)); %#ok<*AGROW>
            hgt = hell(k)-undul(l);
        
            %Pressure, Temperature at the height of the grid
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
            Hs1 = Hs_grid(indx(l));
            redh = hgt - Hs1;

            %lapse rate of the temperature in [degree / m]
            dTl(l) = dT_grid(indx(l),1) + ...
                     dT_grid(indx(l),2)*cosfy + dT_grid(indx(l),3)*sinfy + ...
                     dT_grid(indx(l),4)*coshy + dT_grid(indx(l),5)*sinhy; 

            %Temperature reduction to station height
            Tl(l) = T0 + dTl(l)*redh - 273.15;

            %Virtual temperature
            Tv = T0*(1+0.6077*Ql(l));  
            c = gm*dMtr/(Rg*Tv);
            
            %Pressure in [hPa]
            pl(l) = (p0*exp(-c*redh))/100;
            
            %Hydrostatic coefficient ah
            ahl(l) = ah_grid(indx(l),1) + ...
                     ah_grid(indx(l),2)*cosfy + ah_grid(indx(l),3)*sinfy + ...
                     ah_grid(indx(l),4)*coshy + ah_grid(indx(l),5)*sinhy;
            
            %Wet coefficient aw
            awl(l) = aw_grid(indx(l),1) + ...
                     aw_grid(indx(l),2)*cosfy + aw_grid(indx(l),3)*sinfy + ...
                     aw_grid(indx(l),4)*coshy + aw_grid(indx(l),5)*sinhy;
					 
			%Water vapor decrease factor la - added by GP
			lal(l) = la_grid(indx(l),1) + ...
					 la_grid(indx(l),2)*cosfy + la_grid(indx(l),3)*sinfy + ...
					 la_grid(indx(l),4)*coshy + la_grid(indx(l),5)*sinhy;
					 
			% mean temperature of the water vapor Tm - added by GP
			Tml(l) = Tm_grid(indx(l),1) + ...
					 Tm_grid(indx(l),2)*cosfy + Tm_grid(indx(l),3)*sinfy + ...
					 Tm_grid(indx(l),4)*coshy + Tm_grid(indx(l),5)*sinhy;
					 		 
			% water vapor pressure in hPa - changed by GP
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
        
        %Undulation
        R1 = dnpod2*undul(1)+dnpod1*undul(2);
        R2 = dnpod2*undul(3)+dnpod1*undul(4);
        undu(k) = dnlon2*R1+dnlon1*R2;
		
		%Water vapor decrease factor la - added by GP
        R1 = dnpod2*lal(1)+dnpod1*lal(2);
        R2 = dnpod2*lal(3)+dnpod1*lal(4);
        la(k) = dnlon2*R1+dnlon1*R2;
		
		%Mean temperature of the water vapor Tm - added by GP
        R1 = dnpod2*Tml(1)+dnpod1*Tml(2);
        R2 = dnpod2*Tml(3)+dnpod1*Tml(4);
        Tm(k) = dnlon2*R1+dnlon1*R2;
                    
    end %//if bilinear == 0
    
end%//for k = 1:nstat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF MET PARAMETERS EXTRACTION

%--------------------------------------------------------------------------
%(2)***COMPUTE ZENITH TROPOSPHERIC DELAYS USING EXTRACTED MET PARAMETERS
%--------------------------------------------------------------------------
%NOTE:
%      DAVIS ET AL DRY MODEL REQUIRES +VE ORTOMETRIC HEIGHT TO PERFORM
%--------------------------------------------------------------------------

%COMPUTE CORRECTION FACTOR FOR LOCAL GRAVITY ACCELERATION    
dgref = 1-0.00266.*cos(2.*lat)-0.00028e-3.*(hell-undu);

gm    = 9.784 .* dgref;%Acceleration Gravity at the Atmospheric column in m/s^2

%**************COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
%              =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%USE OF DAVIS ET AL (1985) DRY MODEL   

 ZHD =(0.0022768 .* p)./dgref;  %Zenith Hydrostatic or dry delay [m]    

%**************COMPUTE ZENITH WET DELAY(ZHD)
%              =-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%USE OF  Askne and Nordius (1987) WET MODEL      
%*****DEFINE CONSTANTs
k1 = 77.604; 				   % K/hPa
k2 = 64.79; 				   % K/hPa
Mw = 18.0152;%Molar Mass of Water in [g/mol]
Md = 28.9644; %Molar Mass of Dry Air  in [g/mol]
k2p= k2 - k1*Mw/Md;          % K/hPa  [called k2 prime]
k3 = 377600; 				   % KK/hPa          
dMtr = 28.965*10^-3;%molar mass of dry air in kg/mol    
R = 8.3143;%universal gas constant in J/K/mol         
Rd = R/dMtr ; %specific gas constant for dry air

%ASSIGNMENT
lambda = la;

%*****COMPUTE ZENITH WET DELAY(ZWD)
ZWD = ((1e-6.*(k2p + k3./Tm).*Rd)./(gm.*(lambda + 1))).*e;%[m]
                
%****FINALLY, COMPUTE ZENITH TOTAL DELAY(ZTD)
ZTD = ZHD + ZWD;%Zenith Total Delay [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GPT2w_1_x_1.m  %%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%D.SUBROUTINE TO COMPUTE VIENNA MAPPING FUNCTION(VMF)
function [VMFh,VMFw] = VMF1(ah,aw,UTCtime,lat,hell,satELEV)

%GET SIZE OF USER INPUT TIME & POSITION
nrow_time=size(UTCtime,1);
nrow_pos=size(lat,1);

%1.******UTCdate
Yr  = UTCtime(:,1);%get Year
Mn  = UTCtime(:,2);%get Month
Day = UTCtime(:,3);%get Day

%ASSIGNMENT
Yr1  = Yr;
Mn1  = Mn;
Day1 = Day;


if isequal(nrow_time,nrow_pos)
      
   if isequal(Yr,Yr1) && isequal(Mn,Mn1) && isequal(Day,Day1)
       
      %1.*****INITIALIZING OUTPUT VARIABLEs
      satel=zeros(size(satELEV,1),size(satELEV,2));%Assign zeros of nxm to satel
      [VMFh,VMFw]=deal(satel);%copy the contents of satel to MFh,MFw
      
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
      
                satEL = satELEV(j,i) * pi / 180; %Elevation Angles converted to radian
                
                %Call the "vmf1.m" function
                [VMFh(j,i),VMFw(j,i)] = vmf1(ah(i),aw(i),UTCtime(i,:),lat(i),hell(i),satEL);
                
             end
             
          end
          
   else %for different times
       
      [nrow,ncol]=size(satELEV);%Get size of satELEV
       ncol_ah=size(ah,2); %Get the # of Columns in ah (indication of diff sets of time)
       VMFh = cell(nrow,ncol,ncol_ah);%Create a cell array of VMFh Output
       VMFw = deal(VMFh);%Create copy of VMFh 
       
      for k =1:ncol_ah %Loop over VMF3 mapping factors for different sets of time
          
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
      
                satEL = satELEV(j,i) * pi / 180; %Elevation Angles converted to radian
                
                 %Call the "vmf1.m" function
                [VMFh{j,i,k},VMFw{j,i,k}] = vmf1(ah(i,k),aw(i,k),UTCtime(k,:),lat(i),hell(i),satEL);
                
             end
         end
      end
      
   end %//if isequal(Yr,Yr1) && isequal(Mn,Mn1) && isequal(Day,Day1)
   
else  
     %FIND NUMBER OF ROWS & COLUMNs IN satELEV
     [nrow,ncol]=size(satELEV);%Get size of satELEV  
     
     ncol_ah=size(ah,2); %Get the # of Columns in ah 
      
     if ncol_ah > 1 %(Indication for different sets of time)
          
       VMFh = cell(nrow,ncol,ncol_ah);%Create a cell array of VMFh Output
       VMFw = deal(VMFh);%Create copy of VMFh 
       
       for k =1:ncol_ah %Loop over sets if time
          
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
      
              satEL = satELEV(j,i) * pi / 180; %Elevation Angles converted to radian
              
              %Call the "vmf1.m" function 
              [VMFh{j,i,k},VMFw{j,i,k}] = vmf1(ah(i,k),aw(i,k),UTCtime(k,:),lat(i),hell(i),satEL);
                
          end 
          
          end 
          
       end 
       
     else %for a single time info
         %INITIALIZE OUTPUT
         satel=zeros(size(satELEV,1),size(satELEV,2));%Assign zeros of nxm to satel
         
         [VMFh,VMFw]=deal(satel);%copy the contents of satel to MFh,MFw

            
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
      
             satEL = satELEV(j,i) * pi / 180; %Elevation Angles converted to radian 
             
             %Call the "vmf1.m" function
             [VMFh(j,i),VMFw(j,i)] = vmf1(ah(i),aw(i),UTCtime,lat(i),hell(i),satEL);
         
         
           end
           
        end
        
     end %//if ncol_ah > 1
                
end %//if isequal(nrow_time,nrow_pos)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF VMF1.m  %%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%D_1.SUBROUTINE TO COMPUTE VIENNA MAPPING FUNCTION1(VMF1)
function [VMFh,VMFw] = vmf1(ah,aw,UTCtime,lat,hgt,satELEV)

%**************************************************************************
%***DESCRIPTION:
%              This subroutine determines the Tropospheric Mapping         *
%              Function using the Vienna Mapping Functions 1(VMF1) for     * 
%              Specific site given receiver/station position in [lat hgt], *
%              Satellites Elevation Angle(el),and receiver reception time  *
%              in the utc time format                                      *                                       
%USAGE:                                                                    *
%      General:                                                            *
%              [VMFh,VMFw] = vmf1(ah,aw,UTCtime,lat,hgt,satELEV)           *       

%****INPUTs:                                                               *
%1.         ah: hydrostatic coefficient a (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)
%2.         aw: wet coefficient a         (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)                           *
%3.    UTCtime: UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *
%4.        lat: station ellipsoidal latitude in [radians]                  *
%5.        hgt: station ellipsoidal height in [meters]                     *
%6.    satELEV: Satellites Elevation Angle in [degrees]                    *
% 
%***OUTPUTs:                                                               *
%         VMFh:Vienna Hydrostatic mapping function                         *
%         VMFw:Vienna Wet mapping function                                 *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% REFERENCE:                                                               *
%           Troposphere mapping functions for GPS and very long baseline   *  
%           interferometry from European Centre for Medium-Range Weather   * 
%           Forecasts operational analysis data,J. Geoph. Res., Vol. 111,..* 
%           B02406, doi:10.1029/2005JB003629.                              *

%Original codes by Boehm et al (2006)                                      *
%Modified by: Samuel Osah,Msc Geomatic Engineering ,2016                   *    
%             Email: osahsamuel@yahoo.ca                                   *
%             Tel:+233 (0)246137410/+233 (0)509438484                      *        
%==========================================================================+
%**************************************************************************+
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

%********COMPUTE MODIFIED JULIAN DATE(MJD)
%Call the "utc2JulianDay_DoY.m" function
[~, dmjd,~]=utc2JulianDay_DoY(Yr,Mn,Day,H,MIN,SECs);
    
%*********COMPUTE DAY OF YEAR(doy)
%NOTE:Reference day is 28 January. This is taken from Niell (1996) 
%to be consistent     
doy = dmjd  - 44239.d0 + 1 - 28;    

    
%(1)***********COMPUTE HYDROSTATIC MAPPING FUNCTION(VMFh)
%*****DEFINE MAPPING COEFFICIENTs
bh = 0.0029;
c0h = 0.062;
if (lat<0)    %if Southern Hemisphere
   phh  = pi;
   c11h = 0.007;
   c10h = 0.002;
else             %if Northern Hemisphere
     phh  = 0.d0;
     c11h = 0.005;
     c10h = 0.001;
end 
    
ch = c0h + ((cos(doy/365.25d0*2*pi + phh)+1)*c11h/2 + c10h)*(1-cos(lat));

sine    = sin(satELEV);
beta    = bh/( sine + ch  );
gamma   = ah/( sine + beta);
topcon  = (1.d0 + ah/(1.d0 + bh/(1.d0 + ch)));
vmf1h   = topcon/(sine+gamma); 
 
%*********HEIGHT CORRECTION FOR HYDROSTATIC PART [Niell, 1996] 
%***DEFINE COEFFICIENTs
a_ht   = 2.53d-5;
b_ht   = 5.49d-3;
c_ht   = 1.14d-3;
hs_km  = hgt/1000.d0;%convert height to km
beta         = b_ht/( sine + c_ht);
gamma        = a_ht/( sine + beta);
topcon       = (1.d0 + a_ht/(1.d0 + b_ht/(1.d0 + c_ht)));
ht_corr_coef = 1.d0/sine - topcon/(sine + gamma);
ht_corr      = ht_corr_coef * hs_km;
VMFh         = vmf1h + ht_corr; 
 
%(2)*********COMPUTE WET MAPPING FUNCTION(VMFw)
%DEFINE COEFFICIENTs
bw = 0.00146;
cw = 0.04391;
beta   = bw/( sine + cw );
gamma  = aw/( sine + beta);
topcon = (1.d0 + aw/(1.d0 + bw/(1.d0 + cw)));
VMFw   = topcon/(sine+gamma);  
%******************************************END OF vmf1.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%F.*SUBROUTINE TO COMPUTE JULIAN, MODIFIED DAY(JD,MJD) &  DAY OF YEAR(DoY)
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
 
%F_1.*********SUBROUTINE TO COMPUTE DAY OF YEAR(DoY)
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