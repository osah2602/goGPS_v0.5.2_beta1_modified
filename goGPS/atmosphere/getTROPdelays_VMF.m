%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            This subroutine determines the hydrostatic and wet mapping ...* 
%           function coefficients ah and aw,as well as the zenith delays   *
%           (ZHD,ZWD)from the gridded VMF1 files, as available from:       *
%           http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/VMFG/ for specific   * 
%           sites near the  Earth surface. It is based on a (2.0 × 2.5)    *                                                
%           degree external grid file include only the two MF coefficients * 
%           as well as the hydrostatic and wet portions of ZPD. They are   * 
%           given for four daily epochs (0h, 6h, 12h, 18h UT) and consist  *
%           of four global (2.0×2.5 deg.) grid files of ah, aw, zh and zw. *                                                     
%           Furthermore, since the ZPD (zh,. zw) correspond to mean grid   *
%           heights, a grid file with these mean ellipsoidal heights       *
%           (orography_ell) is also required.The orography_ell file can be * 
%           downloaded from:                                               *
%           [http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/orography_ell].     *                                                    
%           Example of vmf grid file for January 2018 at(0h,6h,12h,18h UT) *
%           include:[VMFG_20180101.H00,VMFG_20180101.H06,VMFG_20180101.H12 *
%           VMFG_20180101.H18]                                             *                                                                
%USAGE:                                                                    *
%      General:                                                            *
%             [ZTD,ZHD,ZWD,STD,SHD,SWD]=getTROPdelays_VMF(Time,ReceiverPos,*
%                                                          SatPos,varargin)*
%      Others:                                                             *
%             [ZTD,ZHD,ZWD,STD,SHD,SWD]=getTROPdelays_VMF(Time,ReceiverPos,*
%                                                  ElevationAngle,varargin)*

%              varargin(variable argument in)  can include,either:         *

%*************1: a)full path of vmf grid file or just the filename when    * 
%                  working in currennt directory.                          *
%            E.g:1)'C:\Users\SAMUEL\Desktop\TropDelay\VMFG_20180101.H00'   *
%                2) 'VMFG_20180101.H00'   [when in current directory]      *
%                    
%                b)full path of orography_ell file or just the filename    *  
%                  when working in currennt directory.                     *
%            E.g:1)'C:\Users\SAMUEL\Desktop\TropDelay\orography_ell'       *
%                2) 'orography_ell'       [when in current directory]      *

%            In that case inputs become five(5). That is:                  *
%            [ZTD,ZHD,ZWD,STD,SHD,SWD]=getTROPdelays_VMF(Time,ReceiverPos, *
%                                               SatPos,E.g:1:a)1,E.g:1:b)1)*
%         OR:                                                              *
%            [ZTD,ZHD,ZWD,STD,SHD,SWD]=getTROPdelays_VMF(Time,ReceiverPos, *
%                                               SatPos,E.g:1:a)2,E.g:1:b)2)*
%            vmf_grid_file=varargin{1}                                     *
%            orography_ell=varargin{2}                                     *
%         OR:                                                              *
%            inputs become four(4) with the vmf grid file: That is:        *
%            [ZTD,ZHD,ZWD,STD,SHD,SWD]=getTROPdelays_VMF(Time,ReceiverPos, *
%                                                         SatPos,E.g:1:a)1)*
%         or:                                                              *
%            [ZTD,ZHD,ZWD,STD,SHD,SWD]=getTROPdelays_VMF(Time,ReceiverPos, *
%                                                         SatPos,E.g:1:a)2)*
%      Where:                                                              *
%            vmf_grid_file=varargin{1}                                     *
%NOTE:                                                                     *
%     The vmf1 grid provides the mapping function coefficients "a" given   *
%     on a global grid with 2.0 degrees  sampling from north to south and  *
%     2.5 degrees sampling from west to east, and for each parameter there *
%     are four files per day, i.e. at 0, 6, 12, and 18 UT, and they are    *
%     stored in yearly directories.In addition to the "a" coefficients,    * 
%     the hydrostatic and wet zenith delays are provided on the grid, too. *  
%     Their values are in m. These zenith delays correspond to the         * 
%     ellipsoidal heights given in the file orography_ell which can be 
%     downloaded:[http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/orography_ell].*
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%**THE FUNCTION CALLS THE FOLLOWING SUBROUTINE(s):                         *
%1.[X,Y,Z] = geo2xyz(latitude,longitude,height,RefEllipsoid);              *
%2.[latRAD,longRAD,h,latDEC,longDEC] = xyz2LLH(X,Y,Z,RefEllipsoid);        *
%3.[a,finv] = Elipsoidpara(RefEllipsoid);                                  *
%4.[Az_rad,El_rad,Az_deg,El_deg,D]=satAzimuthElevation(UserXYZ,SatXYZ,...  *
%5.[VMFh,VMFw] = VMF(ah,aw,UTCtime,lat,hgt,satEL)                          *  
%6.[ah,aw] = vmfCOE_gpt2w(UTCtime,lat,lon,gridFile,GridPoints,it)          *
%7.[ah,aw,ZHD,ZWD] = Readvmf_grid(UTCtime,lat,lon,hgt,VMF1_grid_file,...   *
%                                                           orography_file)*
%8.[FileName, FilePath]= checkGRID(varargin)                               *
%9.[JD, MJD,DoY]=utc2JulianDay_DoY(Year,Month,Day,Hour,Minute,Seconds)     *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%      Generally, the function accepts FOUR(4) sets of inputs:             *
%1.        Time --> Receiver reception time in[year,Month,Day,Hour,Min,Sec]*
%2.    ReceiverPos--> Receiver position in either Latitude,Longitude &     *
%                     height or ECEF(XYZ) Coordinates                      *
%3.    SatPos  -----> Satellite position(s) in ECEF(XYZ) Coordinates       *
%4.   varargin-----> Which contains a vmf grid file e.g.'VMFG_20180101.H00'*  

%OTHER CONSIDERATIONs:                                                     *

%1.*******Satellites Elevation Angle(satEL) can also be used in place of...*
%         Satellites Position (SatPos).  format should be (nx1)or (n x m)  * 
%     NOTE:
%          Elevation Angles should be decimal degrees(eg:26.981.78.102)    *

%4.*******Where only Ellipsoidal Height is provided as input for receiver  *      
%         position(ReceiverPos), the subroutine returns an error message...* 
         
%5.*******Again, for empty ([]) or not a number (nan) entry/input for ...  *
%         ReceiverPos and Time, the subroutine returns an error message    *
 
%6.*******For no entry / input,error message is given and the program      *
%         Terminates / stop working given empty matrices([]) or zeros as   * 
%         outputs.                                                         *          
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
%                --------------------------------------------------                                      *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-                                                                 *
%OUTPUT:                                                                   *
%1.     ZTD => Zenith Total Tropospheric Delay in meters                   *
%2.     ZHD => Zenith Hydrostaic Tropospheric Delay in meters              *
%3.     ZWD => Zenith Wet Tropospheric Delay  in meters                    *
%4.     STD => Slant Total Tropospheric Delay in meters                    *
%5.     SHD => Slant Hydrostaic Tropospheric Delay in meters               *
%6.     SWD => slant Wet Tropospheric Delay  in meters                     *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================+
% REFERENCE:                                                               *
%           Kouba, J. (2008), Implementation and testing of the gridded    *
%           Vienna Mapping Function 1 (VMF1). J. Geodesy, Vol. 82:193-205, * 
%           DOI: 10.1007/s00190-007-0170-0                                 *

%Original codes by Daniel Landskron (2017/06/28)                           *
%Modified by: Samuel Osah,Msc Geomatic Engineering ,2016                   *    
%             Email: osahsamuel@yahoo.ca                                   *
%             Tel:+233 (0)246137410/+233 (0)509438484                      * 
%==========================================================================
%**************************************************************************+
%**************************************************************************+

function [ZTD,ZHD,ZWD,STD,SHD,SWD] = getTROPdelays_VMF(Time,ReceiverPos,SatPos,varargin)
                                                                                                    
%**********CHECK INPUTs & REFORMAT INPUT DATA
switch nargin
    
    case {5,4} %When all inputs are provided 
          
         %******CHECK FOR EMPTY([]) & Not a Number(nan) Position INPUTs
         flag_empty1=0;%flag indicating ReceiverPos input is not empty([]) 
         flag_empty2=0;%flag indicating SatPos input is not empty([]) 
         flag_nan1=0;%flag indicating ReceiverPos input is not a number(nan)
         flag_nan2=0;%flag indicating SatPos input is not a number(nan)
         
         %ASSIGN FLAGs FOR VARIOUS entry Type
         if isempty(ReceiverPos), flag_empty1=1;
            
         elseif  isnan(ReceiverPos), flag_nan1=1; 
         end 
            
         if isempty(SatPos), flag_empty2=1;
            
         elseif  isnan(SatPos), flag_nan2=1; 
         end 
        
         %*******ISSUE VARIOUS ERROR MESSAGEs
         if flag_empty1==1 & flag_empty2==1
           %ISSUE AN ERROR MESSAGE for [] inputs
           beep %Give a beep sound
           errmsg{1}=sprintf('Wrong Data Inputs. Empty ( [ ] ) inputs for Receiver & Satellite Positions.\n'); 
           errmsg{2}='Please Provide Positions & Try Again.';
           errordlg(errmsg,'Coordinate(s) Input Error','modal') 
           return
           
         elseif (flag_empty1==1 & (flag_empty2==0 | flag_nan2==0))                             
                %ISSUE AN ERROR MESSAGE for [] input
                beep %Give a beep sound
                errmsg{1}=sprintf('Wrong Data Input. Empty ( [ ] ) input for Receiver Position(s).\n'); 
                errmsg{2}='Please Provide Receiver Position(s) & Try Again.';
                errordlg(errmsg,'Coordinate(s) Input Error','modal') 
                return
                
         elseif flag_empty1==0 & flag_empty2==1
              %ISSUE AN ERROR MESSAGE for [] input
              beep %Give a beep sound
              errmsg{1}=sprintf('Wrong Data Input. Empty ( [ ] ) input for Satellite(s) Position.\n'); 
              errmsg{2}='Please Provide Satellite(s) Position & Try Again.';
              errordlg(errmsg,'Coordinate(s) Input Error','modal') 
              return 
              
         elseif flag_nan1==1 & flag_nan2==1
               %ISSUE AN ERROR MESSAGE for nan inputs
               beep %Give a beep sound
               errmsg{1}=sprintf('Wrong Data Inputs. Inputs should be numeric.\n');
               errmsg{2}='Please Provide Numeric Position data & Try Again.';
               errordlg(errmsg,'Coordinate(s) Input Error','modal') 
               return
              
         elseif flag_nan1==1 & flag_nan2==0                                
               %ISSUE AN ERROR MESSAGE for nan input
               beep %Give a beep sound
               errmsg{1}=sprintf('Wrong Data for Receiver Position. Input should be numeric.\n');
               errmsg{2}='Please Provide Numeric Satellite(s) Position data & Try Again.';
               errordlg(errmsg,'Coordinate(s) Input Error','modal')
               return
              
         elseif flag_nan1==0 & flag_nan2==1
                %ISSUE AN ERROR MESSAGE for nan input
                beep %Give a beep sound
                errmsg{1}=sprintf('Wrong Data for Satellite(s) Position. Input should be numeric.\n');
                errmsg{2}='Please Provide Numeric Receiver Position data & Try Again.';
                errordlg(errmsg,'Coordinate(s) Input Error','modal')
                return    
              
         else  % IF INPUTS ARE NOT EMPTY([]), then....
            
             %***ROUND TO THE NEAREST THOUSAND 
             t1=roundn(ReceiverPos,4); 
             t2=roundn(SatPos,3);
                  
             if  size(SatPos,2)==1 %if the number of columns in SatPos is 1                                                      
               if  all(t2(:,1))==0                       
                 satEL_deg=SatPos;%Assigning SatPos to satEL(Satellite Elevation)
                 satEL_rad=satEL_deg.* pi / 180;%Elevation Angles converted to radian
                 
               end   %if all(t2(:,1))==0 
               
             else   
                Satpos=SatPos;%Satellite position in XYZ                   
             end   %\\if size(SatPos,2)==1
              
             if size(ReceiverPos,2)==1 
               
               if all(t1(:,1))==0                   
                                   
                %ISSUE ERROR MESSAGE for HEIGHT INPUT 
                beep%Give a beep sound 
                errmsg{1}=sprintf('Insuficient Input for Receiver / Station Position .\n');               
                errmsg{2}='Please Check file & Try Again.';
                errordlg(errmsg,'Coordinate(s) Input Error','modal')  
                return
               end   %//if (all(t1(:,1))==0 || all(t1(:,1))~=0)                
            else 
                Rpos=ReceiverPos;%Receiver position(XYZ/LAT LONG h) 
                        
             end %//if size(ReceiverPos,2)==1
            
        end %//if flag_empty1==1 & flag_empty2==1
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
          
        end %if exist('Rpos','var')
        
        %*******CONVERT Lat,Long,h IF ANY TO XYZ   
        if exist('lat','var') || exist('lon','var')%Check if lat/lon exist
          
          %Call the 'geo2xyz.m' function 
          [X,Y,Z] = geo2xyz(lat,lon,h); 
          
          %*****CONVERT COORD IN DEGREEs TO RADIAN (deg2rad)
          degrad = pi/180.0; %deg2rad conversion factor
          
          latRAD=lat.*degrad; %latitude in radian
          longRAD=lon.*degrad; %longitude in radian
          
          %SAVE COORDs for use
          setappdata(0,'lat',latRAD)
          setappdata(0,'lon',longRAD)
          setappdata(0,'ht',h)
                    
          %**CONVERT USER POSITION(XYZ) TO lat,long,h
        elseif (exist('X','var')|| exist('Y','var') || exist('Z','var'))
               
              %Call the 'xyz2LLH.m' function
              [latRAD,longRAD,h] = xyz2LLH(X,Y,Z);  
              
              %SAVE COORDs for use
              setappdata(0,'lat',latRAD)
              setappdata(0,'lon',longRAD)
              setappdata(0,'ht',h)
                  
        end %exist('lat','var') || exist('lon','var')
        
          %*******COMPUTE SATELLITE ELEVATION ANGLE IF NOT PROVIDED                   
        if ~exist('satEL_deg','var')%If Satellite Elevation is not given,then 
                                    %Compute Satellite Elevation with the given  
                                    %Receiver & Satellite positions
           
           if (exist('Rpos','var') &&  exist('Satpos','var'))                                  
             %Call the "satAzEl.m" function 
             [satEL_deg,satEL_rad]=satAzEl([X Y Z],Satpos);%compute satellite elevation angle             
           end
         
        end %if ~exist('satEL','var'
                                                                       
    otherwise
             %ISSUE ERROR MESSAGE INPUT IS ONE
              beep%Give a beep sound 
              errmsg{1}=sprintf('Insuficient Data Input / Wrong Data Input format .');
              errmsg{2}='';
              errmsg{3}='Please Check file / Data format & Try Again.';
              errordlg(errmsg,'Input Error','modal')
              
              %Return empty ([]) outputs
              STD=[]; SHD=[]; SWD=[]; ZTD=[]; ZHD=[]; ZWD=[];
              
              return                          
end %switch nargin                                                                                                  
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       
%*****CHECK FILE TYPE IN varargin IS/ARE VMF GRID FILE (S)
%***GET LENGTH OF varargin
lv=length(varargin);

switch lv
    
    case 0
        errmsg{1}=sprintf('No input for VMF Grid file.\n');
        errmsg{2}='Please Provide VMF Grid file & Try Again.';    
        errordlg(errmsg,'Grid file input Error','modal')
        return
        
    case 1 %For single input in varargin
        file=varargin{1};%Get the file or full file path
        
        if ~isempty(file)

        if ischar(file) %if file is a character
          
          %SEARCH FOR GRID FILE IN DIRECTORY
          [fileName, filePath]= searchGRID(file);
                            
          %Compare first 4 characters of strings or matche words('VMFG') 
          %or that optionally end with the extension ".H00" or ".H06"
          %or ".H12" or or ".H18"
          if strncmpi(fileName,'VMFG',4) | any ([regexpi(fileName, '\w*.H00$'),...
                                           regexpi(fileName, '\w*.H06$'), regexpi(fileName, '\w*.H12$'),...
                                           regexpi(fileName, '\w*.H18$')])
             vmf_grid=filePath;
                                                         
          else  %if file is not VMF grid file 
                beep %Give a beep sound
                errmsg{1}=sprintf('The File  " %s "  isn''t VMF Grid file.\n',fileName);                 
                errmsg{2}='Please ensure that file is a VMF Grid file & Try Again.';    
                errordlg(errmsg,'Grid file input Error','modal')                  
                return                                                           
          end 
                                               
        else %if file/input is numeric
             beep %Give a beep sound
             errmsg{1}=sprintf('Wrong Data Input.\n');
             errmsg{2}='Please ensure that file is a VMF Grid file & Try Again.';    
             errordlg(errmsg,'Grid file input Error','modal')                  
             return                                                         
        end %if ischar(file) 
        
        else %if input is empty
             beep %Give a beep sound
             errmsg{1}=sprintf('Wrong Data Input. Empty ( [ ] ) input for VMF Grid file.\n');
             errmsg{2}='Please Provide VMF Grid file & Try Again.';    
             errordlg(errmsg,'Grid file input Error','modal')
             return
        end %if ~isempty(file)
        
    case 2 %if inputs are two
        
        file1=varargin{1};
        file2=varargin{2};
        
        if isnumeric (file1) & isnumeric (file2)
           beep %Give a beep sound
           errmsg{1}=sprintf('Wrong Data Inputs.\n');
           errmsg{2}='Please Provide VMF Grid file & Try Again.';    
           errordlg(errmsg,'Grid file input Error','modal')
           return
           
        elseif ~isnumeric (file1) & isnumeric (file2) 
              %SEARCH FOR GRID FILE IN DIRECTORY
              [fileName, filePath]= searchGRID(file1);
                            
             %Compare first 4 characters of strings or matche words('VMFG') 
             %or that optionally end with the extension ".H00" or ".H06"
             %or ".H12" or or ".H18"
             if strncmpi(fileName,'VMFG',4) | any ([regexpi(fileName, '\w*.H00$'),...
                                           regexpi(fileName, '\w*.H06$'), regexpi(fileName, '\w*.H12$'),...
                                           regexpi(fileName, '\w*.H18$')])
                vmf_grid=filePath;
                                                         
             else  %if file is not VMF grid file 
                 beep %Give a beep sound
                 errmsg{1}=sprintf('The File  " %s "  isn''t VMF Grid file.\n',fileName);                 
                 errmsg{2}='Please ensure that file is a VMF Grid file & Try Again.';    
                 errordlg(errmsg,'Grid file input Error','modal')                  
                 return                                                           
             end
             
        elseif  isnumeric (file1) & ~isnumeric (file2) 
              %SEARCH FOR GRID FILE IN DIRECTORY
              [fileName, filePath]= searchGRID(file2);
                            
             %Compare first 4 characters of strings or matche words('VMFG') 
             %or that optionally end with the extension ".H00" or ".H06"
             %or ".H12" or or ".H18"
             if strncmpi(fileName,'VMFG',4) | any ([regexpi(fileName, '\w*.H00$'),...
                                           regexpi(fileName, '\w*.H06$'), regexpi(fileName, '\w*.H12$'),...
                                           regexpi(fileName, '\w*.H18$')])
                vmf_grid=filePath;
                                                         
             else  %if file is not VMF grid file 
                 beep %Give a beep sound
                 errmsg{1}=sprintf('The File  " %s "  isn''t VMF Grid file.\n',fileName);                 
                 errmsg{2}='Please ensure that file is a VMF Grid file & Try Again.';    
                 errordlg(errmsg,'Grid file input Error','modal')                  
                 return                                                           
             end   
                      
        else %if Grid file and orography_ell file are provided
             %SEARCH FOR GRID & orography_ell FILEs IN DIRECTORY
             [fileName1, filePath1]= searchGRID(file1); 
             [fileName2, filePath2]= searchGRID(file2); 
             
             %*****CHECK FILE 1
             %Compare first 4 characters of strings or matche words('VMFG') 
             %or that optionally end with the extension ".H00" or ".H06"
             %or ".H12" or or ".H18"
               if strncmpi(fileName1,'VMFG',4) | any ([regexpi(fileName1, '\w*.H00$'),...
                                                 regexpi(fileName1, '\w*.H06$'), regexpi(fileName1, '\w*.H12$'),...
                                                 regexpi(fileName1, '\w*.H18$')])
                  vmf_grid=filePath1; 
                 
               elseif strncmpi(fileName1,'orography_ell',13)| strncmpi(fileName1,'orography',9)
                      orography_ell=filePath1;               
               else
                   flag1=1;%flag to indicate file isn't required file
                   errmsg{1}=sprintf('The File  " %s "  isn''t required file.\n',fileName1);                                     
               end
               
               %*****CHECK FILE 2
               %Compare first 4 characters of strings or matche words('VMFG') 
               %or that optionally end with the extension ".H00" or ".H06"
               %or ".H12" or or ".H18"
               if strncmpi(fileName2,'VMFG',4) | any ([regexpi(fileName2, '\w*.H00$'),...
                                                 regexpi(fileName2, '\w*.H06$'), regexpi(fileName2, '\w*.H12$'),...
                                                 regexpi(fileName2, '\w*.H18$')])
                  vmf_grid=filePath2; 
                 
               elseif strncmpi(fileName2,'orography_ell',13)| strncmpi(fileName2,'orography',9)
                      orography_ell=filePath2;               
               else
                   flag2=1;%flag to indicate file isn't required file
                   errmsg{2}=sprintf('The File  " %s "  isn''t required file.\n',fileName2);                                     
               end
               
               if exist('flag1','var') & exist('flag2','var')
                  beep %Give a beep sound                                
                  errmsg{3}='Please ensure that files are required VMF files & Try Again.';    
                  errordlg(errmsg,'Grid file input Error','modal') 
                  return
               end
               
        end % if isnumeric (file1) & isnumeric (file2)
                     
end %switch lv
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%************EXTRACT TROPOHERIC DELAYs FROM VMF GRID FILE
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   
if exist('vmf_grid','var') & exist('orography_ell','var')                
      %COMPUTE ZENITH DELAYs & MAPPING coefficients from vmf grid  model
      %Call the "Readvmf_grid.m" function     
      [ZTD,ZHD,ZWD,ah,aw]=Readvmf_grid(UTCtime,latRAD,longRAD,h,vmf_grid, orography_ell);
      
elseif  exist('vmf_grid','var') & ~exist('orography_ell','var')  
        %COMPUTE ZENITH DELAYs & MAPPING coefficients from vmf grid  model
        %Call the "Readvmf_grid.m" function      
        [ZTD,ZHD,ZWD,ah,aw]=Readvmf_grid(UTCtime,latRAD,longRAD,h,vmf_grid);                                                         
end 

%*****************COMPUTE SLANT TROPOSPHERIC DELAYs
if any( [exist('satEL_deg','var') exist('satEL_rad','var')]) 
%(1)*****COMPUTE VIENNA MAPPING FUNCTIONs(VMF)
[MFh,MFw] = VMF(ah,aw,UTCtime,latRAD,h,satEL_rad);

%***************NOW THE SLANT DELAYs
[STD, SHD, SWD] =SlantTropDelay(ZHD,ZWD,MFh,MFw);

else
    %Create zero matrix
     STD=zeros(size(ZHD,1),size(ZHD,2));
     SHD=zeros(size(ZHD,1),size(ZHD,2));
     SWD=zeros(size(ZHD,1),size(ZHD,2));
 end

%***********************END OF getTROPdelays_VMF.m ***********************    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%A.SUBROUTINE TO COMPUTE SLANT TROPOSPHERIC DELAY (SATELLITE IN VIEW)
function [STD, SHD, SWD] =SlantTropDelay(ZHD,ZWD,MFh,MFw)                                                 
                                             
%DESCRIPTION:
%            "SlantTropDelay" Computes Slant Tropospheric Delay from       * 
%             Zenith Hydrostaic and Wet Tropospheric Delays Extracted from *
%             Vienna  Mapping function 1(VMF1) grid file using the Vienna  * 
%             Mapping function 1(VMF1).                                    *
%******SYNTAX:                                                             *
%             [STD, SHD, SWD] =SlantTropDelay(ZHD,ZWD,MFh,MFw)             *
%******INPUT:
%            ZHD = Zenith Hydrostatic Delay                                *
%            ZWD = Zenith Wet Delay                                        *
%            MFh = Hydrostatic Mapping Function                            *
%            MFw = Wet Mapping Function                                    *
%NOTE:                                                                     *
%      MFh or MFw is in the format (nxm);                                  *
%Where:                                                                    *
%      n = Number of rows representing each satellite                      *
%      m = Number of columns representing each station/receiver            *
%    
%****OUTPUT:
%           SHD = Slant Hydrostatic Delay                                  *
%           SWD = Slant Wet Delay                                          *
%           STD = Slant Total Delay                                        *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%*********COMPUTE SLANT TROPOSPHERIC DELAY (SATELLITE IN VIEW)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%(1)*****INITIALIZING OUTPUT VARIABLEs
SHD=zeros(size(MFh,1),size(MFh,2));%Assign zeros of nxm to SHD
[STD,SWD]=deal(SHD);%copy the contents of SHD to all the requested outputs
                                                  
%(2)*****************LOOP OVER ZENITH DELAYs
%(2.1)FIND NUMBER OF ROWS & COLUMNs IN MFh
[nrow,ncol]=size(MFh);

for i=1:ncol %Loop over the Number of rows(corresponds to Receiver positions)
    
   for j=1:nrow %Loop over MF VALUEs(corresponds to sat elevations)
    
       %COMPUTE Slant Hydrostatic Delay(SHD)
       SHD(j,i)=ZHD(i,1).*MFh(j,i);

       %COMPUTE Slant Wet Delay(ZWD)
       SWD(j,i)=ZWD(i,1).*MFw(j,i);

       %****COMPUTE TOTAL SLANT TROPOSPHERIC DELAY(STD)
       STD(j,i)=SHD(j,i)+ SWD(j,i);            
                  
   end %for j=1:nrow
          
end    %for i=1:ncol

%***********************END OF SlantTropDelay.m ***********************    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%B.SUBROUTINE TO READ VMF GRID FILE & EXTRACT COEFFICIENTs(ah,aw)
function [ZTD,ZHD,ZWD,ah,aw] = Readvmf_grid(UTCtime,lat,lon,hgt,VMF1_grid_file,...
                                                            orography_file)                                                     
%**************************************************************************
%DESCRIPTION:
%           This subroutine determines the hydrostatic and wet mapping ... * 
%           function coefficients ah and aw,as well as the zenith delays   *
%           (ZHD,ZWD)from the gridded VMF1 files, as available from:       *
%           http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/VMFG/ for specific   * 
%           sites near the  Earth surface. It is based on a (2.0 × 2.5)    *                                                
%           degree external grid file include only the two MF coefficients * 
%           as well as the hydrostatic and wet portions of ZPD. They are   * 
%           given for four daily epochs (0h, 6h, 12h, 18h UT) and consist  *
%           of four global (2.0×2.5 deg.) grid files of ah, aw, zh and zw. *                                                     
%           Furthermore, since the ZPD (zh,. zw) correspond to mean grid   * 
%           heights, a grid file with these mean ellipsoidal heights       *
%           (orography_ell) is also required.The orography_ell file can be * 
%           downloaded from:                                               *
%           [http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/orography_ell].     *                                                    
%           Example of vmf grid file for January 2018 at(0h,6h,12h,18h UT) *
%           include:[VMFG_20180101.H00,VMFG_20180101.H06,VMFG_20180101.H12 *
%           VMFG_20180101.H18]                                             *
%           ***************************************************************
%           ***************************************************************
%           On the temporal scale, the values from the two surrounding NWM * 
%           epochs are linearly interpolated to the respective mjd.        *
%           In the horizontal, a bilinear interpolation is done for the    * 
%           mapping function coefficients as well as for the zenith delays.*
%           In the vertical,on the one hand the height correction by       *
%           Niell(1996) is applied in order to "lift" the hydrostatic      *
%           mapping function from zero height to hell; On the other hand,  *
%           to "lift" the zenith delays from the respective heights of the *
%           grid points (orography_ell) to that of the desired location.   *
%           specific formulae as suggested by Kouba (2008) are applied.    *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%USAGE:
%      [ZTD,ZHD,ZWD,ah,aw] = Readvmf_grid(UTCtime,lat,lon,hell,VMF1_grid_file...
%                                                          ,orography_file)
%INPUTs:
%1.     UTCtime:.........UTC time in [Year,Month,Day,Hour,Minute,Seconds]  *
%2.     lat: ............station ellipsoidal latitude in [radians]         *
%3.     lon: ............station ellipsoidal longitude in [radians]        *
%4.     hgt: ............station ellipsoidal height in [meters]            *
%5.    VMF1_grid_file:...VMF grid file eg:'VMFG_20180101.H00'              *
%6.    orography_file:...ellipsoidal orography. eg:'orography_ell'         *

%OUTPUTs:                                                                  *
%        ZHD ............... Zenith Hydrostatic Delay, valid at hell       *
%        ZWD ............... Zenith Wet Delay, valid at hell               *
%        ZTD ................Zenith Total Delay, valid at hell             *
%        ah ............... hydrostatic mapping coefficient, valid at hell *
%        aw ............... wet mapping coefficient, valid at hell         *

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% REFERENCE:                                                               *
%           Kouba, J. (2008), Implementation and testing of the gridded    *
%           Vienna Mapping Function 1 (VMF1). J. Geodesy, Vol. 82:193-205, * 
%           DOI: 10.1007/s00190-007-0170-0                                 *

%Original codes by Daniel Landskron (2017/06/28)                           *
%Modified by: Samuel Osah,Msc Geomatic Engineering ,2016                   *    
%             Email: osahsamuel@yahoo.ca                                   *
%             Tel:+233 (0)246137410/+233 (0)509438484                      *        
%==========================================================================+
%**************************************************************************+
%**************************************************************************+

 switch nargin
     case 6
         orogfile=orography_file;
     case  5
         orogfile=[];
 end 
   
%****GET VARIOUS INPUT PARAMETERs
%1.******UTCtime
Yr=UTCtime(:,1);%get Hour
Mn=UTCtime(:,3);%get Month
Day=UTCtime(:,3);%get Day
H=UTCtime(:,4);%get Hour
MIN=UTCtime(:,5);%get Minute
SECs=UTCtime(:,6);%get Seconds

vmf_grid=VMF1_grid_file;

%INITIALIZE OUTPUTs
ZHD=zeros(size(lat,1),1);
[ZWD,ZTD,ah,aw]=deal(ZHD);%copy content of ZHD to ZWD,ah,aw

%lOOP OVER STATIONs
for l=1:length(lat)
    
   %***CONVERT lat and lon TO degrees    
   latDeg = lat(l).*180/pi; %[degrees]
   lonDeg = lon(l).*180/pi; %[degrees]                                                     
   
  for j=1:size(UTCtime,1)%lOOP OVER TIME       
    %********COMPUTE MODIFIED JULIAN DATE(MJD)
    %Call the "utc2JulianDay_DoY.m" function
    [~, MJD ,~]=utc2JulianDay_DoY(Yr(j),Mn(j),Day(j),H(j),MIN(j),SECs(j));
  end 

%****FIND THE TWO SURROUNDING EPOCHs
if mod(MJD,0.25)==0
    MJD_all = MJD;
else
    MJD_int = floor(MJD*4)/4 : 0.25 : ceil(MJD*4)/4;
    MJD_all = [MJD MJD_int];
end

%*******ONLY +VE LONGITUDE IN [degrees] is REQUIRED
if lonDeg < 0
    lon = lon + 2*pi;
    lonDeg = (lonDeg + 360);
end

%****READ OROGRAPHY_ell FILE
if ~isempty(orogfile)
  flag_orog=0; %flag to indicate  orogfile is not empty 
  fileID = fopen(orogfile);
  orography_ell_temp = textscan(fileID,'%f%f%f%f%f%f%f%f%f%f','HeaderLines',1,'CollectOutput',1);
  fclose(fileID);
  orography_ell = reshape(orography_ell_temp{1}',13650,1);
  orography_ell(isnan(orography_ell)) = []; % to get rid of all NaN values
  orography_ell(145:145:length(orography_ell))=[];% delete every 145th value, as these coincide (lon=0=360)
else
    flag_orog=1;%flag to indicate  orogfile does not exist
end    
    
%(3)************* FIND INDICES OF FOUR(4) SURROUNDING GRID POINTs

%*****FIND COORDs (LAT,LON) OF THE SURROUNDING GRID POINTs
lat_int_deg = [floor(latDeg/2)*2 ceil(latDeg/2)*2];
lon_int_deg = [floor(lonDeg/2.5)*2.5 ceil(lonDeg/2.5)*2.5];

%******NOW, FIND THE INDICES OF THESE SURROUNDING GRID POINTS
num_lon = 144;
ind_lat_int_deg = -(lat_int_deg-90)/2; %auxiliary variable
ind_lon_int_deg = lon_int_deg/2.5+1;   %auxiliary variable

index(1) = ind_lat_int_deg(1)*num_lon+ind_lon_int_deg(1);
index(2) = ind_lat_int_deg(1)*num_lon+ind_lon_int_deg(2);
index(3) = ind_lat_int_deg(2)*num_lon+ind_lon_int_deg(1);
index(4) = ind_lat_int_deg(2)*num_lon+ind_lon_int_deg(2);

%(4)******* READ THE DATA AND PERFORM A LINEAR TIME INTERPOLATION 
            %FROM THE SURROUNDING TWO EPOCHS READ IN WITH TEXTSCAN, BUT 
            %ONLY UP TO MAXIMUM INDEX, EVERYTHING BEFORE WILL BE TREATED AS 
            %HEADERLINES
            
%******OPEN THE VMFG FILE            
fileID = fopen(vmf_grid);%GET FILE ID
                        %Only read data up to the maximum index in order to save time
VMF1_data_all = textscan(fileID,'%f%f%f%f%f%f',max(index),'CommentStyle','!','CollectOutput',1);   
                      
                        %Reduce to the indices of the surrounding grid points
VMF1_data = cellfun(@(c) c(index,:),VMF1_data_all,'UniformOutput',false);

%CLOSE THE OPENED FILE
fclose(fileID);

%******INITIALIZE GRID OUTPUT
%NOTE: The vmf grid file has 6 columns i.e.[(lat lon ah aw zhd zwd)]
VMF1_data_int_h0 = zeros(4,7);%Grid values @ grid height (7TH column for grid pressure values)
VMF1_data_int_h1 = zeros(4,9);%Grid values reduced to site/station height
                              %(7TH column for site pressure values)
%DO THE LINEAR TIME INTERPOLATION FOR EACH ARGUMENT; THE RESULTS ARE THE 
%VMF1 VALUES FOR THE SURROUNDING GRID POINTS AT THE TIME OF THE MEASUREMENT
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
iv_ind = 1:4;
if length(MJD_all)==1 %if the observation epoch coincides with an NWM epoch
    VMF1_data_int_h0(iv_ind,1:6) = VMF1_data{1}(iv_ind,1:6);
    
else %otherwise,perform the linear interpolation
    iv_line = 1:6;
    VMF1_data_int_h0(iv_ind,iv_line) = VMF1_data{1}(iv_ind,iv_line) + (VMF1_data{2}(iv_ind,iv_line)-VMF1_data{1}(iv_ind,iv_line))*(mjd-mjd_int(1))/(mjd_int(2)-mjd_int(1));   % the appendix 'h0' means that the values are valid at zero height
end

%ASSIGNING THE FIRST FOUR(4) COLUMNs(They are equal)
VMF1_data_int_h1(:,1:4) = VMF1_data_int_h0(:,1:4);

%(5) BRING MFH, MFW, ZHD AND ZWD OF THE SURROUNDING GRID POINTS TO THE... 
%                  RESPECTIVE HEIGHT OF THE LOCATION

if isequal(flag_orog,0)
    
%(a)*********ZENITH HYDROSTATIC DELAY(ZHD)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%    (1)Convert the hydrostatic zenith delay at grid height to the 
%    respective pressure value using station latitude and height
%NOTE:
%    to be exact, the latitudes of the respective grid points would have to 
%    be used instead of the latitude of the station (lat). However,the loss
%    of accuracy is only in the sub-micrometer range.
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%USE DAVIS ET AL(1985) FORMULAR TO COMPUTE PRESSURE
VMF1_data_int_h0(iv_ind,7) = (VMF1_data_int_h0(iv_ind,5)/0.0022768) .* (1-0.00266*cos(2*lat(1))-0.28*10^-6*orography_ell(index));   

%(2)Transfer the pressure from each grid height to site height using
%   Berg(1948) pressure formular
VMF1_data_int_h1(iv_ind,7) = VMF1_data_int_h0(iv_ind,7).*(1-0.0000226.*(hgt(1)-orography_ell(index))).^5.225;   

%(3)NOW Obtain ZHD at station/site height
VMF1_data_int_h1(iv_ind,5) = 0.0022768*VMF1_data_int_h1(iv_ind,7) / (1-0.00266*cos(2*lat(1))-0.28*10^-6*hgt(1));   % (3) convert the lifted pressure to zhd again (as proposed by Kouba, 2008)

%(b)*********ZENITH WET DELAY(ZWD)
%Use a simple exponential decay approximation function
VMF1_data_int_h1(iv_ind,6) = VMF1_data_int_h0(iv_ind,6) .* exp(-(hgt(1)-orography_ell(index))/2000);

else
    %(a)*********ZENITH HYDROSTATIC DELAY(ZHD)
    VMF1_data_int_h1(iv_ind,5) =VMF1_data_int_h0(iv_ind,5);
    
    %(b)*********ZENITH WET DELAY(ZWD)
    VMF1_data_int_h1(iv_ind,6) =VMF1_data_int_h0(iv_ind,6);
   
end %if isequal(flag_orog,0)

%(c)*********VMF COEFFICIENTs [ah, aw]
hcoe = VMF1_data_int_h0(iv_ind,3);%Hydrostatic coefficient
wcoe = VMF1_data_int_h0(iv_ind,4);%Wet coefficient

VMF1_data_int_h1(iv_ind,8) = hcoe;
VMF1_data_int_h1(iv_ind,9) = wcoe;

%(6)******* PERFORM THE BILINEAR INTERPOLATION

if length(unique(index)) == 1   % if the point is directly on a grid point
    
    zhd = VMF1_data_int_h1(1,5);
    zwd = VMF1_data_int_h1(1,6);
     hcoe = VMF1_data_int_h1(1,8);
     wcoe = VMF1_data_int_h1(1,9);

else
    %BILINEAR INTERPOLATION (INTERPRETED AS TWO 1D LINEAR INTERPOLATIONS 
    %FOR LAT AND LON, BUT PROGRAMMED WITHOUT SUBFUNCTIONS) 
    
    %(a)*******LINEAR INTERPOLATION FOR LONGITUDE
    if ~isequal(VMF1_data_int_h1(1,2), VMF1_data_int_h1(2,2))%if longitude must be interpolated (that is, the point does not have a longitude on the interval [0:2.5:357.5])
        zhd_lon1 = VMF1_data_int_h1(1,5) + (VMF1_data_int_h1(2,5)-VMF1_data_int_h1(1,5))*(lonDeg-VMF1_data_int_h1(1,2))/(VMF1_data_int_h1(2,2)-VMF1_data_int_h1(1,2));
        zhd_lon2 = VMF1_data_int_h1(3,5) + (VMF1_data_int_h1(4,5)-VMF1_data_int_h1(3,5))*(lonDeg-VMF1_data_int_h1(3,2))/(VMF1_data_int_h1(4,2)-VMF1_data_int_h1(3,2));
        zwd_lon1 = VMF1_data_int_h1(1,6) + (VMF1_data_int_h1(2,6)-VMF1_data_int_h1(1,6))*(lonDeg-VMF1_data_int_h1(1,2))/(VMF1_data_int_h1(2,2)-VMF1_data_int_h1(1,2));
        zwd_lon2 = VMF1_data_int_h1(3,6) + (VMF1_data_int_h1(4,6)-VMF1_data_int_h1(3,6))*(lonDeg-VMF1_data_int_h1(3,2))/(VMF1_data_int_h1(4,2)-VMF1_data_int_h1(3,2));
         ah_lon1 = VMF1_data_int_h1(1,8) + (VMF1_data_int_h1(2,8)-VMF1_data_int_h1(1,8))*(lonDeg-VMF1_data_int_h1(1,2))/(VMF1_data_int_h1(2,2)-VMF1_data_int_h1(1,2));
         ah_lon2 = VMF1_data_int_h1(3,8) + (VMF1_data_int_h1(4,8)-VMF1_data_int_h1(3,8))*(lonDeg-VMF1_data_int_h1(3,2))/(VMF1_data_int_h1(4,2)-VMF1_data_int_h1(3,2));
         aw_lon1 = VMF1_data_int_h1(1,9) + (VMF1_data_int_h1(2,9)-VMF1_data_int_h1(1,9))*(lonDeg-VMF1_data_int_h1(1,2))/(VMF1_data_int_h1(2,2)-VMF1_data_int_h1(1,2));
         aw_lon2 = VMF1_data_int_h1(3,9) + (VMF1_data_int_h1(4,9)-VMF1_data_int_h1(3,9))*(lonDeg-VMF1_data_int_h1(3,2))/(VMF1_data_int_h1(4,2)-VMF1_data_int_h1(3,2));

    else %if the station coincides with the longitude of the grid
        zhd_lon1 = VMF1_data_int_h1(1,5);
        zhd_lon2 = VMF1_data_int_h1(3,5);
        zwd_lon1 = VMF1_data_int_h1(1,6);
        zwd_lon2 = VMF1_data_int_h1(3,7);
        ah_lon1 = VMF1_data_int_h1(1,8);
        ah_lon2 = VMF1_data_int_h1(3,8);
        aw_lon1 = VMF1_data_int_h1(1,9);
        aw_lon2 = VMF1_data_int_h1(3,9);
    end
    
    %*****LINEAR INTERPOLATION FOR LATITUDE
    if ~isequal(VMF1_data_int_h1(1,1), VMF1_data_int_h1(3,1))
        zhd = zhd_lon1 + (zhd_lon2-zhd_lon1)*(latDeg-VMF1_data_int_h1(1,1))/(VMF1_data_int_h1(3,1)-VMF1_data_int_h1(1,1));
        zwd = zwd_lon1 + (zwd_lon2-zwd_lon1)*(latDeg-VMF1_data_int_h1(1,1))/(VMF1_data_int_h1(3,1)-VMF1_data_int_h1(1,1));
        hcoe = ah_lon1 + (ah_lon2-ah_lon1)*(latDeg-VMF1_data_int_h1(1,1))/(VMF1_data_int_h1(3,1)-VMF1_data_int_h1(1,1));
        wcoe = aw_lon1 + (aw_lon2-aw_lon1)*(latDeg-VMF1_data_int_h1(1,1))/(VMF1_data_int_h1(3,1)-VMF1_data_int_h1(1,1));
    else %if the station coincides with the latitude of the grid
        zhd = zhd_lon1;%Zenith Hydrostatic Delay
        zwd = zwd_lon1;%Zenith Wet Delay
        hcoe = ah_lon1;%Hydrostatic Coefficient
        wcoe = aw_lon1;%Wet Coefficient
    end
      
end 
%ASSIGN VALUEs FOR ALL STATIONs
ZHD(l,1)=zhd;%Zenith Hydrostatic Delay
ZWD(l,1)=zwd;%Zenith Wet Delay
ah(l,1)=hcoe;%Hydrostatic Coefficient
aw(l,1)=wcoe;%Wet Coefficient

%COMPUTE ZENITH TOTAL DELAY
ZTD(l,1)=ZHD(l,1) + ZWD(l,1);

end %for l=1:length(lat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF Readvmf_grid.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-     
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%C.SUBROUTINE TO COMPUTE VIENNA MAPPING FUNCTION(VMF)
function [VMFh,VMFw] = VMF(ah,aw,UTCtime,lat,hgt,satEL)

%**************************************************************************
%***DESCRIPTION:
%              This subroutine determines the Tropospheric Mapping         *
%              Function using the Vienna Mapping Functions 1(VMF1) for     * 
%              Specific site given receiver/station position in [lat hgt], *
%              Satellites Elevation Angle(el),and receiver reception time  *
%              in the utc time format                                      *                                       
%USAGE:                                                                    *
%      General:                                                            *
%              [VMFh,VMFw] = VMF(ah,aw,UTCtime,lat,hgt,satEL)              *       

%****INPUTs:                                                               *
%1.         ah: hydrostatic coefficient a (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)
%2.         aw: wet coefficient a         (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)                           *
%3.    UTCtime: UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *
%4.        lat: station ellipsoidal latitude in [radians]                  *
%5.        hgt: station ellipsoidal height in [meters]                     *
%6.      satEL: Satellites Elevation Angle in [radians]                    *
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
Yr=UTCtime(:,1);%get Hour
Mn=UTCtime(:,3);%get Month
Day=UTCtime(:,3);%get Day
%2.UTCtime
H=UTCtime(:,4);%get Hour
MIN=UTCtime(:,5);%get Minute
SECs=UTCtime(:,6);%get Seconds

%******LOOP OVER STATIONs & SATELLITEs ELEVATION

%****DEFINE pi
pi = 3.14159265359d0;

%1.*****INITIALIZING OUTPUT VARIABLEs
satel=zeros(size(satEL,1),size(satEL,2));%Assign zeros of nxm to satel
[VMFh,VMFw]=deal(satel);%copy the contents of satel to MFh,MFw

%FIND NUMBER OF ROWS & COLUMNs IN satEL
[nrow,ncol]=size(satEL);

for i=1:ncol %Loop over Station
    
   for j=1:nrow %Loop over Satellites 
              
      %1.COMPUTE ZENITH ANGLE(zd) FROM ELEVATION ANGLEs       
      zd =(pi/2)-satEL(j,i); % zd in [radian]
                                                          
    for k=1:size(UTCtime,1)%lOOP OVER TIME       
       %********COMPUTE MODIFIED JULIAN DATE(MJD)
       %Call the "utc2JulianDay_DoY.m" function
       [~, dmjd,~]=utc2JulianDay_DoY(Yr(k),Mn(k),Day(k),H(k),MIN(k),SECs(k));
    
       %*********COMPUTE DAY OF YEAR(doy)
       %NOTE:Reference day is 28 January. This is taken from Niell (1996) 
       %to be consistent     
       doy = dmjd  - 44239.d0 + 1 - 28;    
    end 
    
    %***********COMPUTE HYDROSTATIC MAPPING FUNCTION(VMFh)
    %*****DEFINE MAPPING COEFFICIENTs
    bh = 0.0029;
    c0h = 0.062;
    if (lat(i)<0)    %if Southern Hemisphere
      phh  = pi;
      c11h = 0.007;
      c10h = 0.002;
    else            %if Northern Hemisphere
        phh  = 0.d0;
        c11h = 0.005;
        c10h = 0.001;
    end
    
    ch = c0h + ((cos(doy/365.25d0*2*pi + phh)+1)*c11h/2 + c10h)*(1-cos(lat(i)));

    sine   = sin(pi/2.d0 - zd);
    beta   = bh/( sine + ch  );
    gamma  = ah(i)/( sine + beta);
    topcon = (1.d0 + ah(i)/(1.d0 + bh/(1.d0 + ch)));
    vmf1h   = topcon/(sine+gamma); 
 
    %*********HEIGHT CORRECTION FOR HYDROSTATIC PART [Niell, 1996] 
    %***DEFINE COEFFICIENTs
    a_ht = 2.53d-5;
    b_ht = 5.49d-3;
    c_ht = 1.14d-3;
    hs_km  = hgt(i)/1000.d0;%convert height to km
    beta         = b_ht/( sine + c_ht);
    gamma        = a_ht/( sine + beta);
    topcon       = (1.d0 + a_ht/(1.d0 + b_ht/(1.d0 + c_ht)));
    ht_corr_coef = 1.d0/sine - topcon/(sine + gamma);
    ht_corr      = ht_corr_coef * hs_km;
    VMFh(j,i)        = vmf1h + ht_corr; 
 
    %(2)*********COMPUTE WET MAPPING FUNCTION(VMFw)
    %DEFINE COEFFICIENTs
    bw = 0.00146;
    cw = 0.04391;
    beta   = bw/( sine + cw );
    gamma  = aw(i)/( sine + beta);
    topcon = (1.d0 + aw(i)/(1.d0 + bw/(1.d0 + cw)));
    VMFw(j,i)   = topcon/(sine+gamma);

   end %for j=1:nrow
   
end %for k=1:nrow
%******************************************END OF VMF.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%D.SUBROUTINE TO COMPUTE JULIAN,MODIFIED DAY(JD,MJD) &  DAY OF YEAR(DoY)
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
 
%D_1.*********SUBROUTINE TO COMPUTE DAY OF YEAR(DoY)
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
     jDay = juliandate(Yr,Mn,D,0,0,0);%Midnight Morning of this Day
     jDayJan1   = juliandate(Yr,1,1,0,0,0);%Midnight Morning of January 1st
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

%E.******SUBROUTINE TO CHECK IF GRID FILE IS IN CURRENT DIRECTORY
function [FileName, FilePath]= searchGRID(varargin)

file=varargin{1};%get input file

if ischar(file)%Check if file is character
  
  %GET LIST OF FILE(S) FROM DIRECTORY
  listFile = dir(file);%List of file(s) 
  
  if isempty(listFile) %if file is not in Current Directory 
     %**CHECK IF MAIN FOLDER HAS FILES & SUBFOLDERS
     folder=cd;%GET MAIN FOLDER
     listFiles = dir(folder);%Check list of file(s)
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
               
       %***OPEN ALL SUB-FOLDERS & GET VMF grid FILES
                
       p=1;%Loop index       
                
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
                              
          end %if Sizesublist(1,1)~=2
          
          Filelists{p,1}= fileLists; %File Lists
          Filepaths{p,1}=Filepath; %File Paths
           
          p=p+1;%Update index
                    
         %COMBINE ALL FILES FROM VARIOUS SUB-FOLDERS  
         filelists = vertcat(Filelists{:}) ;   
         filepaths=vertcat(Filepaths{:}) ; 
          
         %COMPARE INPUT FILE WITH LIST OF FILEs IN SUB-DIRECTORIEs
         for k =1:length(filelists)%Loop over list of files
                        
            if strcmpi(file,filelists(k))%Compare files
              FileName=filelists(k);
              FilePath=filepaths(k);
              
              %CONVERT CELL ARRAY TO CHARACTER
              if iscell(FilePath)
                FilePath=char(FilePath);
              end      
              if iscell(FileName)
                FileName=char(FileName);
              end  
            
            end 
                        
         end %for k =1:length(filelists)
         
         if ~exist('FileName','var') | ~exist('FilePath','var')
            FileName=[];
            FilePath=[];
         end 
         
       end  %for iDir = find(validIndex)
                    
     end %if ~isempty(fileList)
     
  else
      FileName=file;
      FilePath=file;
      
      %CONVERT CELL ARRAY TO CHARACTER
      if iscell(FilePath)
         FilePath=char(FilePath);
      end        
      if iscell(FileName)
         FileName=char(FileName);
      end    
          
  end %if isempty(listFile)
  
end %if ischar(file)
%******************************************END OF searcGRID.m 
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
  
%The function Calls the following Subroutine(s):                             *
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