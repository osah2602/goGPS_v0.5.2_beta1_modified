%**************************************************************************
%DESCRIPTION:                                                              * 
%           "getMETpara_GPT2w" computes surface Meteorological parameters  *
%           such as pressure,temperature,temperature lapse rate,mean ...   *
%           temperature of water vapor,water vapor pressure,water vapour   *
%           decrease factor and geoid undulation for a given latitude,     *
%           height and day of year(doy)                                    *
%           It is based on:                                                *
%           1. a 1 x 1 degree external grid file ('gpt2_1w.grd')           *
%           2. a 5 x 5 degree external grid file ('gpt2_5w.grd') with mean *
%             values as well as sine and cosine amplitudes for the annual  * 
%             and semiannual variation of the coefficients.                *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
%NOTE:  This file goes together with "readGPTgrid.m", which reads the grid *
%       and saves it in cell structures.                                   *
%       +The "readGPTgrid.m" should be ran first to provide the grid input *
%       file (grid)                                                        *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
%USAGE:                                                                    *
%      [P,T,E,TM,lambda,dT,undu,ah,aw]=getMETpara_GPT2w(Time,ReceiverPos,grid,grid_res,Timevar)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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

%2. NOW RUN THE "getMETpara_GPT2w.m" TO COMPUTE MET PARAMETERS             +
%INPUTS(eg):                                                               +
%1.Time        = [2019 1 5]                                                +
%2.ReceiverPos = [6 40 3.7449 -1 34  0.7263  260.2410]                     +
%3.grid        = grid (i.e.extracted from "readGPTgrid.m")                 +
%4.grid_res    = 5                                                         +
%5.Timevar     = 0 i.e. 0 means GPT2w with time variation, 1 means static  +

%==> [P,T,E,TM,lambda,dT,undu,ah,aw]=getMETpara_GPT2w(Time,ReceiverPos,grid,grid_res,Timevar)
%
%i.e.[P,T,E,TM,lambda,dT,undu,ah,aw]=getMETpara_GPT2w([2019 1 5],[6 40 3.7449 -1 34  0.7263  260.2410],grid,5,0)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%**THE FUNCTION CALLS THE FOLLOWING SUBROUTINE(s):                         *
%1.[X,Y,Z] = geo2xyz(latitude,longitude,height,RefEllipsoid);              *
%2.[latRAD,longRAD,h,latDEC,longDEC] = xyz2LLH(X,Y,Z,RefEllipsoid);        *
%3.[a,finv] = Elipsoidpara(RefEllipsoid);                                  *
%5.[P,T,Tm,e,la,undu,dT]=getMETpara(UTCtime,lat,lon,hell,Timevar)          *
%6.[FileName, FilePath] = searchGRID(varargin)                             *
%7.[JD, MJD,DoY]=utc2JulianDay_DoY(Year,Month,Day,Hour,Minute,Seconds)     *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%      Generally, the function accepts Five(5) sets of inputs:              *
%1.           Time: Observation time in[year,Month,Day,Hour,Min,Sec]       *
%2.    ReceiverPos: Receiver position in either Latitude,Longitude &...    *
%                     height or ECEF(XYZ) Coordinates                      *
%3            grid: grid values in cells extracted from the grid file      *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%NOTE: To get the grid values(grid) the "readGPTgrid.m" subroutine has to  *
%      be ran first.
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%4.       grid_res: grid resolution (°) (possible: 1 or 5).indicate 1 for  *
%               1°x1° external grid file and 5 for 5°x5° external grid file*

%5.        Timevar: is an indicator of time variation to either compute Met*
%                   parameters in staic mode or in annual and semiannual   * 
%                   terms. Indicate by specifying 1 or 0.                  *
%              case 1: no time variation but static quantities             *
%              case 0: with time variation (annual and semiannual terms)   *
%NOTE:                                                                     *
%1.   b)HeightS are also given in meters                                   *
%2.   Receiver positions are in the ff formats:                            *
%     a)For Geographic Coordinates (n x m):[Latitude Longitude height(h)]  *
%      eg:DMS:[6 40 21 -1 33 17 187.76]  OR                                *
%          DM:[6 40.35 -1 33.2833 187.76] OR                               *
%         Deg:[6.6725 -1.5547 187.76] Decimal degrees                      *
%*****For Western Longitudes place negative(-) infront of the value        *
% eg:[-1 33 34] or [-001 33 34] or [ 000 -24 16] or [ -1 33.2]OR [-1.3345] *

%*****For Southern Latitudes place negative(-) infront of the nunmber      *
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
%***OUTPUTs:                                                               *
%1.         T : Surface temperature in degree Celcius(°C)                  *
%2.         P : Surface pressure in millibar(mbar)                         *
%3.         E : Surface water vapor pressure in millibar(mbar)             *
%4.        TM : Mean temperature of water vapor in kelvin(k)               *
%5.      beta : temperature `lapse' rate in kelvin per meter[(K/m)]        *
%6.    lambda : water vapour `lapse rate'(dimensionless)]                  * 
%4.        dT :   temperature lapse rate in [degrees per km ]              *
%7.       undu: geoid undulation in meters                                 *
%8.      ah: hydrostatic mapping function coefficient at zero height (VMF1)* 
%9.      aw: wet mapping function coefficient(VMF1)(vector of length nstat)*
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================
%%REFERENCE:                                                               *
%           J. Böhm, G. Möller, M. Schindelegger, G. Pain, R. Weber,...    *
%           Development of an improved blind model for slant delays in the *
%           troposphere (GPT2w),GPS Solutions, 2014,...                    *
%           doi:10.1007/s10291-014-0403-7                                  *

%Original codes by Böhm et al 2014                                         *
%==========================================================================
%Modified by:                                                              +
%            OSAH SAMUEL, MSC GEOMATIC ENGINEERING (PhD STUDENT)           +
%            Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com                +
%            Phone:+233(0)246137410 / +233(0)509438484                     + 
%==========================================================================+
%***************************************************************************

function [P,T,E,TM,lambda,dT,undu,ah,aw]=getMETpara_GPT2w(Time,ReceiverPos,grid,grid_res,Timevar)
                                                      
%**********CHECK INPUTs & REFORMAT INPUT DATA
switch nargin
    
    case {5,4,3,2} %When all inputs are provided        
        
        if (any(nargin == [5,4,3,2]))
            
          if nargin == 5 %for FIVE(5) or all inputs...  
             
             %******Assignments 
             if isempty(grid)%grid data input is empty                
               %Retrieve stored grid 
               grid =getappdata(0,'gpt_grid' ) ;
             end 
                        
             if any([isempty(grid_res),isnan(grid_res)])%grid resolution input is empty                
                grid_res = 5;%set grid resolution to 5
             end
             
             if any([isempty(Timevar) isnan(Timevar)])%Timevar input is empty                
               Timevar = 0;%set time variation incator to 0
             end
             
          elseif nargin == 4 %for ONLY four(4) inputs...
              
                 %******Assignments 
                 
                 if isempty(grid)%grid data input is empty                
                    %Retrieve stored grid 
                    grid =getappdata(0,'gpt_grid' ) ;
                 end 
            
                 if any([isempty(grid_res),isnan(grid_res)])%grid resolution input is empty                
                    grid_res = 5;%set grid resolution to 5
                 end 
                 
                  Timevar = 0;%set time variation incator to 0
             
                  
          elseif nargin ==3 %for ONLY THREE(3) inputs... 
              
                 %******Assignments
                 
                 if isempty(grid)%grid data input is empty                
                    %Retrieve stored grid 
                    grid =getappdata(0,'gpt_grid' ) ;
                 end 
                 
                  grid_res = 5;%set grid resolution to 5
                  Timevar = 0;%set time variation incator to 0   
                  
             
          elseif  nargin ==2 %for ONLY two(2) inputs... 
              
                  %******Assignments                            
                  grid =getappdata(0,'gpt_grid' ) ; %Retrieve stored grid 
                  grid_res = 5;%set grid resolution to 5
                  Timevar = 0;%set time variation incator to 0   
                                      
          end  %if nargin == 5     
        
       %*******************CHECK FOR EMPTY([]) INPUTs
          
        if isempty(Time)
            
           %ISSUE ERROR MESSAGE  
           beep%Give a beep sound 
           errmsg0{1}=sprintf('Observation / Reception Time is not provided i.e it''s Empty ([ ]).\n');
           errmsg0{2}='Please provide Observation / reception Time & Try again.';
           errordlg(errmsg0,'Time Input Error','modal')  
             
           %RETURN EMPTY([]) PARAMETERs
           P =[];T =[];E =[];TM =[];lambda =[];dT =[];undu =[];ah =[];aw =[];
           return
        end    
        
        if isempty(ReceiverPos)
            
            %ISSUE ERROR MESSAGE 
            beep%Give a beep sound 
            errmsg10{1}=sprintf('Reciever / Site Position is not provided i.e it''s Empty ([ ]).\n');
            errmsg10{2}='Please provide Receiver / Site position & Try again.';
            errordlg(errmsg10,'Reciever Position(s) Input Error','modal')  
           
            %RETURN EMPTY([]) PARAMETERs
            P =[];T =[];E =[];TM =[];lambda =[];dT =[];undu =[];ah =[];aw =[];
             
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
                    P =[];T =[];E =[];TM =[];lambda =[];dT =[];undu =[];ah =[];aw =[];
            
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
                     
          %RETURN EMPTY([]) PARAMETERs
          P =[];T =[];E =[];TM =[];lambda =[];dT =[];undu =[];ah =[];aw =[];
           
          return  
          
       end
                         
       if isempty(grid_res)
                      
          %ISSUE ERROR MESSAGE  
          beep%Give a beep sound 
          errmsg20{1}=sprintf('No input for VMF grid resolution, i.e. input is empty ( [ ] ) .\n');
          errmsg20{2}='Please provide VMF grid resolution & Try again.';
          errordlg(errmsg20,'Grid resolution Input Error','modal') 
                     
          %RETURN EMPTY([]) PARAMETERs
          P =[];T =[];E =[];TM =[];lambda =[];dT =[];undu =[];ah =[];aw =[];
           
          return 
                         
       end %//if isempty(grid_resl)
        %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
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
                                    
          end  %if any([~isempty(Time),~isnan(Time)])
          
        end   %if (any(nargin==[3,2]))
          
        %******************************************************************
        %         EXTRACT USER INPUT RECEIVER Coordinates,.....
        %*****************************************************************   
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
          [X,Y,Z] = geo2xyz(lat,lon,h);  %#ok<*ASGLU>
          
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
        
        end  %if exist('Rpos','var')
                                                                                        
    otherwise
             %ISSUE ERROR MESSAGE INPUT IS ONE
              beep%Give a beep sound 
              errmsg{1}=sprintf('Insuficient Data Input / Wrong Data Input format .');
              errmsg{2}='';
              errmsg{3}='Please Check file / Data format & Try Again.';
              errordlg(errmsg,'Input Error','modal')
              
              %Return empty ([]) outputs
              T=[]; P=[]; TM=[]; E=[]; lambda=[]; dT=[];undu=[]; ah=[];aw=[];
              
              return                          
end %switch nargin                                                                                                  
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%-----------------------------------------------------
%********COMPUTE METEOROLOGICAL PARAMETERS USING GPT2w
%------------------------------------------------------
%Call the "getMETpara.m" Function
[P,T,E,TM,lambda,dT,undu,ah,aw] = getMETpara(UTCtime,latR,lonR,h,grid,grid_res,Timevar); 


%A.******SUBROUTINE TO COMPUTE METEOROLOGICAL PARAMETERS 
function [p,T,e,Tm,la,dT,undu,ah,aw] = getMETpara(UTCtime,lat,lon,hell,grid,grid_res,Timevar)

%********COMPUTE METEOROLOGICAL PARAMETERS USING GPT2w MODEL
%1.****NOW COMPUTE MET PARAMETERs 
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
      %Y= diff(X) calculates differences between adjacent elements of X along 
      %the first array dimension whose size does not equal 1: 
      %The elements of Y are the differences between adjacent elements of X.
      %Y = [X(2)-X(1) X(3)-X(2) ... X(m)-X(m-1)].
      %The not (~) converts the differences to logicals(0 difference->true). 
      %Then all sees if all the differences in a column are 0
      %all(bsxfun(@eq,Yr,Yr(1,:))) ====> (This actually has one advantage, 
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
      p=zeros(nrow_time,1);
      [T,dT,Tm,e,ah,aw,la,undu]=deal(p);%create copy of P in dT,Tm,e,ah,aw,la,undu
   
      for i=1:nrow_time
       
          if grid_res == 1 %if grid resolution is 1 x 1 degree
             %Call the "GPT2w_1_x_1.m" Function
             [p(i,1),T(i,1),dT(i,1),Tm(i,1),e(i,1),ah(i,1),aw(i,1),la(i,1),undu(i,1)] = GPT2w_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid,Timevar); 
          
          elseif grid_res == 5 %if grid resolution is 5 x 5 degree
                 %Call the "GPT2w_5_x_5.m" Function
                [p(i,1),T(i,1),dT(i,1),Tm(i,1),e(i,1),ah(i,1),aw(i,1),la(i,1),undu(i,1)] = GPT2w_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid,Timevar);
        
          end 
      end
      
   else 
       %*****INITIALIZE OUTPUTs 
      p=zeros(nrow_pos,nrow_time);
     [T,dT,Tm,e,ah,aw,la,undu]=deal(p);%create copy of P in T,E,TM,lambda,dT,undu
    
     for i=1:nrow_time %LOOP OVER TIME
        for j=1:nrow_pos %LOOP OVER POSITIONS
            
            %****CHECK GRID RESOLUTION
            
            if grid_res == 1 %if grid resolution is 1 x 1 degree
               %Call the the 1 degree grid version; "GPT2w_1_x_1.m"  
               [p(j,i),T(j,i),dT(j,i),Tm(j,i),e(j,i),ah(j,i),aw(j,i),la(j,i),undu(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);
           
            elseif grid_res == 5 %if grid resolution is 5 x 5 degree
                   %Call the the 5 degree grid version; "GPT2w_5_x_5.m"  
                   [p(j,i),T(j,i),dT(j,i),Tm(j,i),e(j,i),ah(j,i),aw(j,i),la(j,i),undu(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);
                   
            end

           
        end
     end
    
   end %//if isequal(Identical,1)
    
else
    %*****INITIALIZE OUTPUTs 
    p=zeros(nrow_pos,nrow_time);
    [T,dT,Tm,e,ah,aw,la,undu]=deal(p);%create copy of P in T,E,TM,lambda,dT,undu
    
    for i=1:nrow_time %LOOP OVER TIME
        for j=1:nrow_pos %LOOP OVER POSITIONS
            
            %****CHECK GRID RESOLUTION
            
            if grid_res == 1 %if grid resolution is 1 x 1 degree
               %Call the the 1 degree grid version; "GPT2w_1_x_1.m"  
               [p(j,i),T(j,i),dT(j,i),Tm(j,i),e(j,i),ah(j,i),aw(j,i),la(j,i),undu(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);
           
            elseif grid_res == 5 %if grid resolution is 5 x 5 degree
                   %Call the the 5 degree grid version; "GPT2w_5_x_5.m"  
                   [p(j,i),T(j,i),dT(j,i),Tm(j,i),e(j,i),ah(j,i),aw(j,i),la(j,i),undu(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);
                   
            end

           
        end
    end
end
%***********************END OF getMETpara.m ***********************    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%A_1.SUBROUTINE TO COMPUTE METEOROLOGICAL PARAMETERS FROM 5x5 GRID MODEL
function [p,T,dT,Tm,e,ah,aw,la,undu] = GPT2w_5_x_5(UTCtime,lat,lon,hell,grid,Timevar)

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
%      [p,T,dT,Tm,e,ah,aw,la,undu] = GPT2w_5(UTCtime,lat,lon,hell,Timevar) *
%INPUTs:                                                                   *
%1.   UTCtime:  UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *                                                                
%2.       lat:  Station geodetic latitude in radians [-pi/2:+pi/2] (n x 1) *
%3.       lon:  Station geodetic longitude in radians [-pi:pi] or [0:2pi]  *
%4.      hell:  ellipsoidal height in [meters] (n x 1)                     *
%5    Timevar:  case 1: no time variation but static quantities            *
%               case 0: with time variation (annual and semiannual terms)  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% OUTPUTs:
%1.       p: pressure in hPa (vector of length nstat)                      *
%2.       T: temperature in degrees Celsius (vector of length nstat)       *
%3.      dT: temperature lapse rate in degrees per km (vector of length nstat)
%4.      Tm: mean temperature of the water vapor in degrees Kelvin         * 
%5.       e: water vapor pressure in hPa (vector of length nstat)          *
%6.      ah: hydrostatic mapping function coefficient at zero height (VMF1)* 
%7.      aw: wet mapping function coefficient (VMF1) (vector of length nstat)
%8.      la: water vapor decrease factor (vector of length nstat)          *
%9.    undu: geoid undulation in m (vector of length nstat)                *
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GPT2w_5 x 5.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%A_2.******SUBROUTINE TO COMPUTE METEOROLOGICAL PARAMETERS FROM 1x1 GRID MODEL
function [p,T,dT,Tm,e,ah,aw,la,undu] = GPT2w_1_x_1(UTCtime,lat,lon,hell,grid,Timevar)

%**************************************************************************
%DESCRIPTION:
%            This subroutine determines pressure, temperature, temperature *
%            lapse rate, mean temperature of the water vapor, water vapor  * 
%            pressure, hydrostatic and wet mapping function coefficients   * 
%            ah and aw, water vapour decrease factor and geoid undulation  * 
%            for specific sites near the Earth surface. It is based on a   *
%            1 x 1 degree external grid file ('gpt2_1w.grd') with mean     *
%            values as well as sine and cosine amplitudes for the annual   * 
%            and semiannual variation of the coefficients.
%USAGE:
%      [p,T,dT,Tm,e,ah,aw,la,undu] = GPT2w_1_x_1(UTCtime,lat,lon,hell,Timevar) 
%INPUTs:                                                                   *
%1.   UTCtime:  UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *                                                                
%2.       lat:  Station geodetic latitude in radians [-pi/2:+pi/2] (n x 1) *
%3.       lon:  Station geodetic longitude in radians [-pi:pi] or [0:2pi]  *
%4.      hell:  ellipsoidal height in [meters] (n x 1)                     *
%5    Timevar:  case 1: no time variation but static quantities            *
%               case 0: with time variation (annual and semiannual terms)  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% OUTPUTs:
%1.       p: pressure in hPa (vector of length nstat)                      *
%2.       T: temperature in degrees Celsius (vector of length nstat)       *
%3.      dT: temperature lapse rate in degrees per km (vector of length nstat)
%4.      Tm: mean temperature of the water vapor in degrees Kelvin         * 
%5.       e: water vapor pressure in hPa (vector of length nstat)          *
%6.      ah: hydrostatic mapping function coefficient at zero height (VMF1)* 
%7.      aw: wet mapping function coefficient (VMF1) (vector of length nstat)
%8.      la: water vapor decrease factor (vector of length nstat)          *
%9.    undu: geoid undulation in m (vector of length nstat)                *
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GPT2w_1_x_1.m  %%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%B.SUBROUTINE TO COMPUTE JULIAN,MODIFIED DAY(JD,MJD) &  DAY OF YEAR(DoY)
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
    JD  = juliandate(Yr, Mn, D);
    MJD = mjuliandate(Yr, Mn, D);
end

%****COMPUTE DAY OF YEAR
%Call the "DayofYear function.m"
 DoY=DayofYear(Year,Month,Day);
%******************************************END OF utc2JulianDay_DoY.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
 
%B_1.*********SUBROUTINE TO COMPUTE DAY OF YEAR(DoY)
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
     jDayJan1  = juliandate(Yr,1,1,0,0,0);%Midnight Morning of January 1st
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
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-