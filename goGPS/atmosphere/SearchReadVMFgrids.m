%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            This subroutine reads and extract Vienna Mapping Function(VMF)* 
%            Zenith Tropospheric Delays,ZPDs(i.e.ZHD,ZWD and ZTD) and ...  *
%            mapping coefficients(ah,aw) from both gridded VMF1& VMF3 files*
%            as available from: 
%            http://vmf.geo.tuwien.ac.at/trop_products/GRID/ for specific  * 
%            sites near the  Earth surface.The VMF1 grid file is based on  *                                                
%            a 2.0° × 2.5° external grid file while the VMF3 grid file     *
%            is based on either a 1° × 1° OR 5° × 5° grid.Both grids(VMF1 &*
%            VMF3) are given for four daily epochs (0h, 6h, 12h, 18h UT).  *
%            The VMF grid file consist two MF coefficients (ah,aw)as well  * 
%            as the Hydrostatic and Wet portions of ZPD(zh,zw).            *
%            The Zenith Total Delay(ZTD) is the sum of ZHD & ZWD.          *
%            i.e. ZTD = ZHD + ZWD                                          *      
%--------------------------------------------------------------------------
%*******COLUMN DECRIPTION FOR THE VMF GRIDs:                               *
%The VMF grid files have six(6) columns[lat lon ah aw zhd zwd],where:      *
%COLUMN                                                                    *
%      (1)   = latitude (°)                                                *
%      (2)   = longitude (°)                                               *
%      (3)   = hydrostatic "a" coefficient                                 *
%      (4)   = wet "a" coefficient                                         *
%      (5)   = hydrostatic zenith delay (m)                                *
%      (6)   = wet zenith delay (m)                                        *
%--------------------------------------------------------------------------*
%--------------------------------------------------------------------------*
%INPUT:                                                                    *
%       The function accepts FIVE(5) sets of inputs:                       *
%1.     UTCtime : Receiver reception time in[year,Month,Day,Hour,Min,Sec]  *
%       OR         [year,Month,Day]                                        *
%2.     lat     : ellipsoidal latitude (radians)                           *
%3.     lon     : ellipsoidal longitude (radians)                          *
%4.     h       : ellipsoidal height (m)                                   *
%5.VMFgrids     : Structure array of VMF grid files which can be :         *
%--------------------------------------------------------------------------*
%A.    h0files  : Zero(0) Hour VMFG File Name(s) eg:'VMFG_20180101.H00'    * 
%B.    h6files  : Six(6) Hour VMFG File  Name(s) eg:'VMFG_20180101.H06'    *  
%C.    h12files : Twelve(12) Hour VMFG File Name(s) eg:'VMFG_20180101.H12' *  
%D.    h18files : Eighteen(18) Hour VMFG File Name(s)eg:'VMFG_20180101.H18'*
%E.    H0files  : Zero(0) Hour VMFG File Path(s)                           * 
%                 eg:'C:\Users\...data\VMFG_20180101.H00'                  * 
%F.    H6files  : Six(6) Hour VMFG File  Path(s)                           *  
%G.    H12files : Twelve(12) Hour VMFG File Path(s)                        *  
%H.    H18files : Eighteen(18) Hour VMFG File Path(s)                      * 
%I.    Orography: Orography file                                           *
%--------------------------------------------------------------------------+
%--------------------------------------------------------------------------+
%*******NOTE:                                                              +
%            Situation where VMFgrids are not provided,
%            The "TropModel_VMF" subroutine searches Vienna...             +
%            Mapping Function Grid (VMF)files in a given Directory/Folder  +
%            as indicated in the sub-routine. Situation where VMFG files   + 
%            are not found in the default directory/folder,it searches     +   
%            recursively through all Sub-Directories/Folders of the...     +  
%            given Directory. VMFG files are extracted by looping ...      + 
%            through all the listed files in the provided folder or ...    +
%            sub-folders.Finally,if VMFG files are still not found in the..+ 
%            default directory/folder and its sub-folders,the search for   +
%            VMFG files is extended to the current directory.              +
%--------------------------------------------------------------------------
%OTHER CONSIDERATIONs:                                                     *
%1.   Receiver positions are in the ff formats for SINGLE INPUT:           *
%     a)For Geographic Coordinates (n x m):[Latitude Longitude height(h)]  *
%      eg:DMS:[6 40 21 -1 33 17 187.76]  OR                                *
%          DM:[6 40.35 -1 33.2833 187.76] OR                               *
%         Deg:[6.6725 -1.5547 187.76] Decimal degrees                      *

%For Western Longitudes place negative(-) infront of the value             *
% eg:[-1 33 34] or [-001 33 34] or [ 000 -24 16] or [ -1 33.2]OR [-1.3345] *

%For Southern Latitudes place negative(-) infront of the nunmber           *
%  eg:[-8 23 14] or [-008 23 14] or [ 000 -23 26] or [ -8 23.6]OR [-8.2315]*

%     b)For ECEF(XYZ) Coordinates (3 x n) OR (n x 3)matrix :               *
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
%--------------------------------------------------------------------------+
%OUTPUTs                                                                   *
%4.      ZHD : Zenith Hydrostatic Delay (m),valid at h                     *
%5.      ZWD : Zenith Wet Delay (m), valid at h                            *
%6.      ZTD : Zenith Total Delay (m), valid at h                          *
%4.       ah : Hydrostatic mapping function coefficient valid at h         *
%5.       aw : Wet mapping function coefficient valid at h                 *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%Original codes by Daniel Landskron (2017/06/28)                           *
%Modified by:                                                              +
%            OSAH SAMUEL, MSC GEOMATIC ENGINEERING (PhD STUDENT)           +
%            Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com                +
%            Phone:+233(0)246137410 / +233(0)509438484                     +         
%==========================================================================+
%**************************************************************************+
%**************************************************************************+

function [ZHD,ZWD,ZTD,ah,aw,VMF_model] = SearchReadVMFgrids(Time,Lat,Lon,h,VMFgrids)

%**********CHECK & REFORMAT INPUTs DATA
switch nargin
    
   case {5,4,3,2} %Various inputs format 
       
        if (any(nargin==[5,4,3,2])) 
            
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
           
            if (any(nargin==[5,4])) %IF 5 or 4 inputs are provided
                
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
               
               if nargin == 4   %IF 4 INPUTS ARE PROVIDED
                  
                  %CHECK IF h POSITION IS NOT VMF GRID VALUES(VMFgrids)
                  %OR VMFgrids ASSIGNED TO h POSITION----------------------
                  
                  if isstruct(h)%Determine whether h input is structure array
                   
                     h        = zeros(size(Lon,1),1);%Assigning zeros to  heights(h) 
                     VMFgrids = h;%ASSUMING h POSITION TO BE VMFgrids
                     
                  elseif ~isstruct(h)%Determine whether h input is not structure array
                      
                         VMFgrids = [];%Assigning empty matrix ([]) to VMFgrids
                         
                  elseif isempty(h) %Determine whether h array is empty
                         
                          h        = zeros(size(Lon,1),1);%Assigning zeros to  heights(h)
                          VMFgrids = [];%Assigning empty matrix ([]) to VMFgrids
                          
                  end %//if isstruct(h)
                  
               end %//if nargin == 4
               
            elseif (any(nargin==[3,2])) %IF THREE INPUTS ARE PROVIDED
                     
                    if nargin == 3 %IF 3 INPUTS ARE PROVIDED
                            
                       if isstruct(Lon)%Determine whether Lon input is structure array
                          
                           Rpos     = Lat; %ASSIGNING Lat TO Rpos(RECEIVER/SITE POSITON) 
                           VMFgrids = Lon;%ASSUMING Lon POSITION TO BE VMFgrids
                                    
                       else 
                           
                            VMFgrids = []; %Assigning empty matrix ([]) to VMFgrids
                                 
                            if ~isempty(Lon)
                                     
                               if any(ncol1==[3,2,1])
                                     
                                  h = zeros(size(Lon,1),1);%Assigning zeros to  heights(h)
                                       
                                  %*****CHECK LONGITUDE INPUT 
                                  ncol_lon = size(Lon,2); %Get number of longitude entry 
   
                                  switch ncol_lon
                                      case 3     %if input is in DMS,Convert to degrees
                                            lon = dms2degrees(Lon);%Lat in degrees 
                                            
                                      case 2   %if input is in DM,Convert to degrees    
                                            lon = dm2degrees(Lon);%Lat in degrees
                   
                                      otherwise     
                                               lon = Lon;
                                  end  
                                        
                                  %****CHECK LATITUDE INPUT 
                                  ncol_lat= size(Lat,2); %Get number of latitude entry  
   
                                  switch ncol_lat
                                
                                      case 3     %if input is in DMS,Convert to degrees
                                            lat = dms2degrees(Lat);%Lat in degrees 
                                            
                                      case 2   %if input is in DM,Convert to degrees    
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
                                 
                       end %//if isstruct(Lon)
                                       
                     elseif nargin == 2
                             
                            Rpos     = Lat; %ASSIGNING Lat TO Rpos(RECEIVER/SITE POSITON) 
                            VMFgrids = []; %Assigning empty matrix ([]) to VMFgrids
                                 
                    end %//if nargin == 3    
                      
            end %if (any(nargin==[5,4]))
            
        end %//if (any(nargin==[5,4,3,2]))
          
        %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
        
        else
            latD = lat;%Latitude in Decimal degrees
            lonD = lon;%Longitude in Decimal degrees
            
        end   %if exist('Rpos','var') 
        
    otherwise   
             %ISSUE ERROR MESSAGE for 1 INPUT
              beep%Give a beep sound 
              errmsg{1}=sprintf('Insuficient Data Input / Wrong Data Input format .');
              errmsg{2}='';
              errmsg{3}='Please Check inputs / Data format & Try Again.';
              errordlg(errmsg,'Coordinate(s) Input Error','modal')  
              
              %Return empty ([]) outputs
              STD = []; SHD = []; SWD = []; ZTD = []; ZHD = []; ZWD = [];
              ah  = []; aw  = []; VMF_model = [];
              
              return                          
        
end  %switch nargin
%===============================END OF INPUT REFORMATTING
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%CHECK IF VMFgrids IS EMPTY([]) & IMPORT GRIDS USING "SearchVMFgrids.m"
%THIS IS USER FAILS TO PROVIDE VMFgrids IN THE INPUT ARGUMENT

if isempty(VMFgrids)%Determine whether array is empty
    
   %Call the "SearchVMFgrids.m" fxn
   [VMF_grid_found,VMFgrids] = SearchVMFgrids();
   
end

%CHECK IF VMFgrids IS EMPTY([]) & IMPORT GRIDS USING "SearchVMFgrids.m"
%THIS IS USER FAILS TO PROVIDE VMFgrids IN THE INPUT ARGUMENT

if any([all([~isfield(VMFgrids,'H00'),~isfield(VMFgrids,'H06'),~isfield(VMFgrids,'H12'),~isfield(VMFgrids,'H18')]),...
        isempty(VMFgrids)])%Determine whether array is empt
   
   %Call the "SearchVMFgrids.m" fxn
   [VMF_grid_found,VMFgrids] = SearchVMFgrids();
   
end

if any([any([isfield(VMFgrids,'H00'),isfield(VMFgrids,'H06'),isfield(VMFgrids,'H12'),isfield(VMFgrids,'H18')]),...
        ~isempty(VMFgrids)])%Determine whether array is empt
    
%RETRIEVE VMF GRID & OROGRAPHY FILES FROM STRUCTURE ARRAY
%*******00 HOUR
H00     = VMFgrids.H00;%struct of 00 Hour file name & path
h0files = H00.name;%00 hour file name
H0files = H00.path;%00 Hour file path

%*******06 HOUR
H06     = VMFgrids.H06;%struct of 06 Hour file name & path
h6files = H06.name;%06 Hour file name
H6files = H06.path;%06 Hour file path

%*******12 HOUR
H12      = VMFgrids.H12;%struct of 12 Hour file name & path
h12files = H12.name;%12 Hour file name
H12files = H12.path;%12 Hour file path

%*******18 HOUR
H18      = VMFgrids.H18;%struct of 18 Hour file name & path
h18files = H18.name;%18 Hour file name
H18files = H18.path;%18 Hour file path

%*******OROGRAPHY FILE
Oro      = VMFgrids.Orography;%struct of Orography file name & path
orofileN = Oro.name;%Orography file name
orofileP = Oro.path;%Orography file path

else 
    H0files = []; H6files = []; H12files = []; H18files = [];
    
end 

%***********CHECK IF ALL GRID FILES ARE EMPTY([]) 
if all([isempty(H0files),isempty(H6files),isempty(H12files),isempty(H18files)])
   
   %RETURN EMPTY([]) MATRICES AS OUTPUT
   ZHD = []; ZWD = []; ZTD = [];  ah = []; aw = []; VMF_model = [];
        
   return
   
end 
                                                                   
%%%%%%%%%%%%%%%%%%%%%%END OF FILE IMPORT
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*

%*****************CONCANTENATE LAT LON & h ALL TOGETHER 
%1******RECEIVER / STATION POSTION
Rpos  = [latD lonD h]; 

%CHECK THE NUMBER OF STATIONS
nrow_Rpos = size(Rpos,1);  
 
%1.GET [YEAR(Yr) MONTH(Mn) DAY]
Yr  = UTCtime(:,1);%get Year
Mn  = UTCtime(:,2);%get Month
Day = UTCtime(:,3);%get Day
   
%ASSIGNMENT
Date = UTCtime;%UTC TIME

%**************LOOK FOR VMF FILE TYPE(VMF1 OR VMF3)

%RETRIEVE STORED VMF GRID TYPE & RESOLUTION('VMF1(2°x2.5°)','VMF3(1°x1°)','VMF3(5°x5°)')
VMFgrid_type=getappdata(0,'VMFgrid_type');%VMF GRID TYPE(VMF1 or VMF3)
VMFgrid_res = getappdata(0,'VMFgrid_res');%VMF GRID RESOLUTION especially for VMF3(1°x1° OR 5°x5°)
      
%********DERIVE RELATED VMF1 OR VMF3 FILENAME(S)

epoch = [0 6 12 18];%VMF Observation epochs
 
%INITIALIZE OUTPUT
h00 = cell(size(Date,1),1);%create cell array
[h06,h12,h18,h00_txt,h06_txt,h12_txt,h18_txt] = deal(h00);%Make copy of h00 

for q = 1 : size(Date,1)
           
     yr  = Yr(q);
     mn  = Mn(q);
     day = Day(q);
            
     if any([strcmpi(VMFgrid_type,'VMF1'),strfind(VMFgrid_type,'VMF1'),any(VMFgrid_res ==[ 2,2.5])]) 
        %DERIVED VMFG FILENAME USING USER INPUT DATE
        h00{q,1}     = ['VMFG_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(1)))];
        h06 {q,1}    = ['VMFG_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(2)))] ;
        h12{q,1}     = ['VMFG_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(3)))] ;
        h18{q,1}     = ['VMFG_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(4)))];

        %DERIVED VMFG FILENAME WITH '.txt' EXTENSION USING USER INPUT DATE
        h00_txt{q} = ['VMFG_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(1))) '.txt']; 
        h06_txt{q} = ['VMFG_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(2))) '.txt']; 
        h12_txt{q} = ['VMFG_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(3))) '.txt']; 
        h18_txt{q} = ['VMFG_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(4))) '.txt'];
            
     elseif any([strcmpi(VMFgrid_type,'VMF3'),strfind(VMFgrid_type,'VMF3'),any(VMFgrid_res ==[ 1,5])]) 
            %DERIVED VMF3 FILENAME USING USER INPUT DATE
            h00{q,1}     = ['VMF3_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(1)))];
            h06 {q,1}    = ['VMF3_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(2)))] ;
            h12{q,1}     = ['VMF3_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(3)))] ;
            h18{q,1}     = ['VMF3_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(4)))];

            %DERIVED VMF3 FILENAME WITH '.txt' EXTENSION USING USER INPUT DATE
            h00_txt{q} = ['VMF3_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(1))) '.txt']; 
            h06_txt{q} = ['VMF3_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(2))) '.txt']; 
            h12_txt{q} = ['VMF3_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(3))) '.txt']; 
            h18_txt{q} = ['VMF3_' num2str(yr) sprintf('%02s',num2str(mn)) sprintf('%02s',num2str(day)) '.H' sprintf('%02s',num2str(epoch(4))) '.txt'];
                
     end  %//if isequal(VMF1file,1)    
         
end  %//for q = 1 : size(Date,1) 
         
%********SORT THE DERIVED VMF grid FILES & REMOVE DUPLICATE FILES
        
%***********0 HOUR FILES(0h)
[~,I_01,~]=unique(h00);%GET IDENTICAL 0 HOUR FILES INDEX
[~,i_01,~]=unique(h00_txt);
           
%***GET UNIQUE 0h FILES
h00=h00(I_01);% sort in Ascending order
h00_txt=h00_txt(i_01);% sort in Ascending order
                     
%***********6 HOUR FILES(6h)
[~,I_61,~]=unique(h06);%GET IDENTICAL 6 HOUR FILES INDEX
[~,i_61,~]=unique(h06_txt);
           
%***GET UNIQUE 6h FILES
h06=h06(I_61);% sort in Ascending order
h06_txt=h06_txt(i_61);% sort in Ascending order
            
%***********12 HOUR FILES(12h)
[~,I_121,~]=unique(h12);%GET IDENTICAL 12 HOUR FILES INDEX
[~,i_121,~]=unique(h12_txt);
           
%***GET UNIQUE 12h FILES
h12=h12(I_121);% sort in Ascending order
h12_txt=h12_txt(i_121);% sort in Ascending order
                   
%***********18 HOUR FILES(18h)      
[~,I_181,~]=unique(h18);%GET IDENTICAL 18 HOUR FILES INDEX
[~,i_181,~]=unique(h18_txt);
           
%***GET UNIQUE 18h FILES
h18=h18(I_181);% sort in Ascending order
h18_txt=h18_txt(i_181);% sort in Ascending order
             
%COMPARE CORRESPONDING ROWS IN IMPORTED VMF FILEs & DERIVED VMF FILEs &...
%RETURN DATA COMMON TO BOTH IMPORTED & DERIVED VMF FILES  
%--------------------------------------------------------------------------
%NOTE:THE "upper" function is used to Convert VMF string file to 
%     uppercase in case imported VMF file is in lower case.E.G.:
%    'VMFG_20180101.h00' or VMF3_20180101.h00 instead of
%     VMFG_20180101.H00 or VMF3_20180101.H00
%--------------------------------------------------------------------------

%*******EXTRACT CORRESPONDING ROWS OF DATA IN EACH FILE
if ~isempty(h0files)
   [H0fileN,i_h00,i_h0files] = intersect(h00,upper(h0files));%H0 filenames (e.g.VMFG_20180401.H00)
   H0Files = H0files(i_h0files);%H0 filenames & path(e.g.C:\Users\VMFG_20180401.H00)
           
   %******VMF FILES WITH .txt EXTENSION
   [H0fileN_t,i_h00_t,i_h0files_t] = intersect(h00_txt,upper(h0files));%H0 filenames (e.g.VMFG_20180401.H00)
   H0Files_t = H0files(i_h0files_t);%H0 filenames & path(e.g.C:\Users\VMFG_20180401.H00)
     
end 
        
if ~isempty(h6files)
    
    [H6fileN,i_h06,i_h6files] = intersect(h06,upper(h6files));%H0 filenames (e.g.VMFG_20180401.H06)
     H6Files = H6files(i_h6files);%H0 filenames & path(e.g.C:\Users\VMFG_20180401.H06)
           
   %******VMF FILES WITH .txt EXTENSION
   [H6fileN_t,i_h06_t,i_h6files_t] = intersect(h06_txt,upper(h6files));%H0 filenames (e.g.VMFG_20180401.H06)
   H6Files_t = H6files(i_h6files_t);%H0 filenames & path(e.g.C:\Users\VMFG_20180401.H06)
   
end 
        
if ~isempty(h12files)
    
   [H12fileN,i_h12,i_h12files] = intersect(h12,upper(h12files));%H0 filenames (e.g.VMFG_20180401.H12)
    H12Files = H12files(i_h12files);%H0 filenames & path(e.g.C:\Users\VMFG_20180401.H12)
           
   %******VMF FILES WITH .txt EXTENSION
   [H12fileN_t,i_h12_t,i_h12files_t] = intersect(h12_txt,upper(h12files));%H0 filenames (e.g.VMFG_20180401.H12)
   H12Files_t = H12files(i_h12files_t);%H0 filenames & path(e.g.C:\Users\VMFG_20180401.H12)

end 
        
if ~isempty(h18files)
    
   [H18fileN,i_h18,i_h18files] = intersect(h18,upper(h18files));%H0 filenames (e.g.VMFG_20180401.H18)
   H18Files = H18files(i_h18files);%H0 filenames & path(e.g.C:\Users\VMFG_20180401.H18)
           
   %******VMF FILES WITH .txt EXTENSION
   [H18fileN_t,i_h18_t,i_h18files_t] = intersect(h18_txt,upper(h18files));%H0 filenames (e.g.VMFG_20180401.H18)
   H18Files_t = H18files(i_h18files_t);%H0 filenames & path(e.g.C:\Users\VMFG_20180401.H18)

end 
              
%*********COMBINING VARIOUS VMF FILES EXTENSIONS
%1.********0 HOUR FILEs
if ~isempty(H0fileN)
           
   if all([~isempty(H0fileN),~isempty(H0fileN_t)])
      H0filesN = [H0fileN;H0fileN_t];%H0 filenames(e.g.VMFG_20180401.H00) 
      H0filesP = [H0Files;H0Files_t];%H0 filenames & path(e.g.C:\Users\VMFG_20180401.H00)
      index_t  = [i_h00;i_h00_t]; %GET INDEX TO EXTRACT CORRESPONDING DATE
               
   elseif all([~isempty(H0fileN),isempty(H0fileN_t)]) 
          H0filesN = H0fileN;%H0 filenames(e.g.VMFG_20180401.H00) 
          H0filesP = H0Files;%H0 filenames & path(e.g.C:\Users\VMFG_20180401.H00)
          index_t  = i_h00; %GET INDEX TO EXTRACT CORRESPONDING DATE 
          
   elseif  all([isempty(H0fileN),~isempty(H0fileN_t)]) 
           H0filesN = H0fileN_t;%H0 filenames(e.g.VMFG_20180401.H00.txt) 
           H0filesP = H0Files_t;%H0 filenames & path(e.g.C:\Users\VMFG_20180401.H00.txt)
           index_t  = i_h00_t; %GET INDEX TO EXTRACT CORRESPONDING DATE
   
   end %//if ~isempty(H0fileN) & ~isempty(H0fileN_t) 
   
   size_vmf = size(H0filesN,1);
   
else 
    H0filesN = [];
         
end %//if ~isempty(h0files)
         
%2.********6 HOUR FILEs
if ~isempty(H6fileN)
           
   if all([~isempty(H6fileN),~isempty(H6fileN_t)])
      H6filesN = [H6fileN;H6fileN_t];%H6 filenames(e.g.VMFG_20180401.H06) 
      H6filesP = [H6Files;H6Files_t];%H6 filenames & path(e.g.C:\Users\VMFG_20180401.H06)
      index_t  = [i_h06;i_h06_t]; %GET INDEX TO EXTRACT CORRESPONDING DATE 
      
   elseif all([~isempty(H6fileN),isempty(H6fileN_t)])
          H6filesN = H6fileN;%H6 filenames(e.g.VMFG_20180401.H06) 
          H6filesP = H6Files;%H6 filenames & path(e.g.C:\Users\VMFG_20180401.H06)
          index_t  = i_h06; %GET INDEX TO EXTRACT CORRESPONDING DATE 
          
   elseif all([isempty(H6fileN),~isempty(H6fileN_t)])
          H6filesN = H6fileN_t;%H6 filenames(e.g.VMFG_20180401.H06.txt) 
          H6filesP = H6Files_t;%H6 filenames & path(e.g.C:\Users\VMFG_20180401.H06.txt)
          index_t  = i_h06_t; %GET INDEX TO EXTRACT CORRESPONDING DATE 
   
   end %//if ~isempty(H6fileN) & ~isempty(H6fileN_t)   
 
  size_vmf = size(H6filesN,1);
  
else
     H6filesN = [];

  
end %//if ~isempty(h6files)
       
%3.********12 HOUR FILEs
if ~isempty(H12fileN)
           
  if all([~isempty(H12fileN),~isempty(H12fileN_t)])
     H12filesN = [H12fileN;H12fileN_t];%H12 filenames(e.g.VMFG_20180401.H12) 
     H12filesP = [H12Files;H12Files_t];%H12 filenames & path(e.g.C:\Users\VMFG_20180401.H12)
     index_t   = [i_h12;i_h12_t]; %GET INDEX TO EXTRACT CORRESPONDING DATE 
     
  elseif all([~isempty(H12fileN),isempty(H12fileN_t)])
         H12filesN = H12fileN;%H12 filenames(e.g.VMFG_20180401.H12) 
         H12filesP = H12Files;%H12 filenames & path(e.g.C:\Users\VMFG_20180401.H12)
         index_t   = i_h12; %GET INDEX TO EXTRACT CORRESPONDING DATE 
         
  elseif all([isempty(H12fileN),~isempty(H12fileN_t)])
         H12filesN = H12fileN_t;%H12 filenames(e.g.VMFG_20180401.H12.txt) 
         H12filesP = H12Files_t;%H12 filenames & path(e.g.C:\Users\VMFG_20180401.H12.txt)
         index_t   = i_h12_t; %GET INDEX TO EXTRACT CORRESPONDING DATE 
  
  end %if ~isempty(H12fileN) & ~isempty(H12fileN_t)   
   
 size_vmf = size(H12filesN,1);
 
else 
    H12filesN = [];
 
end %//if ~isempty(h12files)
       
%4.********18 HOUR FILEs
if ~isempty(H18fileN)
           
  if all([~isempty(H18fileN),~isempty(H18fileN_t)])
     H18filesN = [H18fileN;H18fileN_t];%H18 filenames(e.g.VMFG_20180401.H18) 
     H18filesP = [H18Files;H18Files_t];%H18 filenames & path(e.g.C:\Users\VMFG_20180401.H18)
     index_t   = [i_h18;i_h18_t]; %GET INDEX TO EXTRACT CORRESPONDING DATE  
     
  elseif all([~isempty(H18fileN),isempty(H18fileN_t)])
         H18filesN = H18fileN;%H18 filenames(e.g.VMFG_20180401.H18) 
         H18filesP = H18Files;%H18 filenames & path(e.g.C:\Users\VMFG_20180401.H18)
         index_t   = i_h18; %GET INDEX TO EXTRACT CORRESPONDING DATE 
                
  elseif all([isempty(H18fileN),~isempty(H18fileN_t)])
         H18filesN = H18fileN_t;%H18 filenames(e.g.VMFG_20180401.H18.txt) 
         H18filesP = H18Files_t;%H18 filenames & path(e.g.C:\Users\VMFG_20180401.H18.txt)
        index_t   = i_h18_t; %GET INDEX TO EXTRACT CORRESPONDING DATE 
  
  end %//if ~isempty(H18fileN) & ~isempty(H18fileN_t)  
 
 size_vmf = size(H18filesN,1);
 
else 
     H18filesN = [];
 
end %//if ~isempty(h6files)

if all([isempty(H0filesN),isempty(H6filesN),isempty(H12filesN),isempty(H18filesN)])
   
    %RETURN EMPTY([]) MATRICES AS OUTPUT
    ZHD = []; ZWD = []; ZTD = []; 
    ah = []; aw = []; VMF_model = [];
   
   beep %Give a beep sound
   errmsg900{1}=sprintf('VMF grid (s) Epoch doesn''t Match Rinex Observation Time.\n');
   errmsg900{2}=sprintf('Please make sure grid (s) Epoch and Rinex observation Time Matches & Try Again.\n');
   errordlg(errmsg900,'VMF grid file Error','modal')
   
   return
   
end

%***EXTRACT CORRESPONDING DATE / TIME
%SORT TIME INDEX IN ASCENDING ORDING
[~,I_time]=unique(index_t );%GET IDENTICAL TIME INDEX

%***GET UNIQUE TIME INDEX
index_t=index_t(I_time);% sort in Ascending order

%GET DATE
date = Date(index_t,:) ;

% %GET NROWS IN date
% nrow_time = size(date,1);

%*****GET ALL ALL CORRESPONDING ROWS IN Date if # ROWS > 1

%EXTRACT ONLY DATE(Year Month Day)
Date1 = Date(:,1:3);%original input date from user
date1 = date(:,1:3);%Extracted correspong date with VMF grid epoch 

%EXTRACT ONLY TIME(Hour Minute Seconds)
UTCtime1 = Date(:,4:6);

%CHECK & EXTRACT SIMILARITY IN date1 & Date1
try
   Date1(date1(:,:)==Date1(:,:));
   UTCtime1(date1(:,:)==Date1(:,:));%get associated/corresponding  time
catch
     date1 = repmat(date1,size(Date,1),1);
     Date1(date1(:,:)==Date1(:,:));
     UTCtime1(date1(:,:)==Date1(:,:));%get associated/corresponding  time
end

%CONCATENATE Date1 & UTCtime1
date = [Date1  UTCtime1];

%GET SIZE OF USER INPUT POSITION & EXTRAXTED TIME  
nrow_time = size(date,1);
nrow_pos = nrow_Rpos;%Assignment
       
%INCLUDE OROGRAPHY FILE
if ~isempty(orofileN)
   oro = orofileP{1};%CHANGE FROM CELL TO CHARACTER 
else 
    oro = [];
end

if isequal(nrow_time,nrow_pos)%IF USER POSITION & OBSERVATION TIME ARE EQUAL[SAME SIZE]
      
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
          
     %LOOP OVER STATIONS
     for i = 1 : nrow_Rpos 
    
        %********EXTRACTING TROPOSPHERIC DELAYS
               
        %COMBINE VARIOUS EPOCHs FROM IMPORTED VMF grid FILEs
        Y  = date(i,1); %GET YEAR
        M  = date(i,2); %GET MONTH
        D  = date(i,3); %GET DAY
        H  = date(i,4);%get Hour
      MIN  = date(i,5);%get Minute
      SECs = date(i,6);%get Seconds
           
        if D > size_vmf 
           
           try
              D = index_t(i);
           catch
                D = index_t;
                
           end
             
        end %if D > size_vmf 
           
        %*********CONVERT TIME IN [Hour Min Sec] TO HOURs               
        UT = H + (MIN / 60) + (SECs / 3600);% Time in Hours
          
        %********COMBINE AVAILABLE VMF GRID EPOCHS(00 06 12 18) FILEs FOR EACH AVAILABLE DAY
        %IF ALL FILES ARE AVAILABLE[00H 06H 12H 18H]
        if all([~isempty(h0files),~isempty(h6files),~isempty(h12files),~isempty(h18files)])
            
           if size(Time,2) > 3 
                  
              %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
              %------------------------------------------------------------------
              %NOTE: 
              %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
              %    I.E 00Hr 06Hr 12Hr & 18Hr
              %------------------------------------------------------------------
              if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                 VMFGfileN = H0filesN(D);
                 VMFGfileP = H0filesP(D);
            
              elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                     VMFGfileN = H6filesN(D);
                     VMFGfileP = H6filesP(D);
                
              elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                     VMFGfileN = H12filesN(D);
                     VMFGfileP = H12filesP(D);
                      
              elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                     VMFGfileN = H18filesN(D);
                     VMFGfileP = H18filesP(D);
                       
              end  
               
           elseif size(Time,2) == 3
                     
                  VMFGfileN = [H0filesN(D);H6filesN(D);H12filesN(D);H18filesN(D)];%#ok<*NASGU> %VMFG FILE WITHOUT PATH OR DIRECTORY
                  VMFGfileP = [H0filesP(D);H6filesP(D);H12filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY
           
           end  
                %IF THREE FILES ARE AVAILABLE[00H 06H 12H ]
        elseif all([~isempty(h0files),~isempty(h6files),~isempty(h12files),isempty(h18files)])
                  
               if size(Time,2) > 3 
                  
                  %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                  %------------------------------------------------------------------
                  %NOTE: 
                  %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                  %    I.E 00Hr 06Hr 12Hr & 18Hr
                  %------------------------------------------------------------------
                  if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                
                  elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                
                  elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                      
                  elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th Hour epoch. System will use the 12th (12) hour grid file.\nOtherwise you can provide the 18th(12) hour grid files & try again\n\n') 
                         
                  end   
               
               elseif size(Time,2) == 3
                     
                       VMFGfileN = [H0filesN(D);H6filesN(D);H12filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                       VMFGfileP = [H0filesP(D);H6filesP(D);H12filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY
    
               end 
                 
             %IF THREE FILES ARE AVAILABLE[00H 06H 18H ]
        elseif  all([~isempty(h0files),~isempty(h6files),isempty(h12files),~isempty(h18files)])
                  
                if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th Hour epoch. System will use the 6th (06) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n') 
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                   end   
               
                elseif  size(Time,2) == 3
                   
                        VMFGfileN = [H0filesN(D);H6filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H0filesP(D);H6filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY
           
                end 
                  
                %IF THREE FILES ARE AVAILABLE[00H 12H 18H ]      
        elseif  all([~isempty(h0files),isempty(h6files),~isempty(h12files),~isempty(h18files)])
                  
                if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H0filesN(D);
                        VMFGfileP = H0filesP(D);
                        
                        beep
                        fprintf('\n\nNo VMF grid file at the 6th Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')   
                        
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                   end   
               
                elseif  size(Time,2) == 3
               
                        VMFGfileN = [H0filesN(D);H12filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H0filesP(D);H12filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY
            
                end 
                 
                %IF THREE FILES ARE AVAILABLE[06H 12H 18H ]       
        elseif  all([isempty(h0files),~isempty(h6files),~isempty(h12files),~isempty(h18files)])
                   
                 if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H6filesN(D);
                     VMFGfileP = H6filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                   end   
               
                 elseif size(Time,2) == 3
                     
                        VMFGfileN = [H6filesN(D);H12filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H6filesP(D);H12filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY   
               
                 end 
                  
                %IF TWO FILES ARE AVAILABLE[00H 06H ]       
        elseif   all([~isempty(h0files),~isempty(h6files),isempty(h12files),isempty(h18files)])
                   
                 if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                         
                   end   
               
                 elseif size(Time,2) == 3
                      
                        VMFGfileN = [H0filesN(D);H6filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H0filesP(D);H6filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY   
           
                 end 
                  
                 %IF TWO FILES ARE AVAILABLE[00H 12H ]
        elseif   all([~isempty(h0files),isempty(h6files),~isempty(h12files),isempty(h18files)])
                 
                 if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                         VMFGfileN = H0filesN(D);
                         VMFGfileP = H0filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                         
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                         
                   end   
               
                 elseif size(Time,2) == 3
                      
                        VMFGfileN = [H0filesN(D);H12filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H0filesP(D);H12filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY   
             
                 end 
                  
                 %IF TWO FILES ARE AVAILABLE[00H 18H ] 
        elseif   all([~isempty(h0files),isempty(h6files),isempty(h12files),~isempty(h18files)])
                  
                 if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                         VMFGfileN = H0filesN(D);
                         VMFGfileP = H0filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                                       
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the 18th(18) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                           
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                   end  
               
                 elseif size(Time,2) == 3
                      
                        VMFGfileN = [H0filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H0filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY 
                 end 
                  
                 %IF TWO FILES ARE AVAILABLE[06H 12H ]
        elseif   all([isempty(h0files),~isempty(h6files),~isempty(h12files),isempty(h18files)])
                  
                 if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H6filesN(D);
                     VMFGfileP = H6filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                          beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                         
                   end  
               
                 elseif size(Time,2) == 3
                      
                        VMFGfileN = [H6filesN(D);H12filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H6filesP(D);H12filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY 
                
                 end 
                  
                 %IF TWO FILES ARE AVAILABLE[06H 18H ] 
        elseif  all([isempty(h0files),~isempty(h6files),isempty(h12files),~isempty(h18files)])
                  
                if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H6filesN(D);
                     VMFGfileP = H6filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the 18th(18) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                   end  
               
                elseif size(Time,2) == 3
                      
                       VMFGfileN = [H6filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                       VMFGfileP = [H6filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY
                  
                end 
                 
                 %IF TWO FILES ARE AVAILABLE[12H 18H ] 
        elseif   all([isempty(h0files),isempty(h6files),~isempty(h12files),~isempty(h18files)])

                 if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H12filesN(D);
                     VMFGfileP = H12filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H12filesN(D);
                        VMFGfileP = H12filesP(D);
                        
                         beep
                         fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                   end  
               
                 elseif size(Time,2) == 3
                      
                        VMFGfileN = [H12filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H12filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY 
                  
                 end 
                  
                 %IF ONLY ONE FILE IS AVAILABLE[00H] 
        elseif  all([~isempty(h0files),isempty(h6files),isempty(h12files),isempty(h18files)])
                  
                if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H0filesN(D);
                        VMFGfileP = H0filesP(D);
                        
                         beep
                         fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H0filesN(D);
                         VMFGfileP = H0filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H0filesN(D);
                         VMFGfileP = H0filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                           
                   end  
               
                elseif  size(Time,2) == 3 
                      
                        VMFGfileN = H0filesN(D);%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = H0filesP(D);%VMFG FILE WITH PATH OR DIRECTORY   
                end 
                  
                 %IF ONLY ONE FILE IS AVAILABLE[06H] 
        elseif  all([isempty(h0files),~isempty(h6files),isempty(h12files),isempty(h18files)])
                  
                if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H6filesN(D);
                     VMFGfileP = H6filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                               
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                                                       
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                           
                   end  
               
                elseif  size(Time,2) == 3 
                      
                        VMFGfileN = H6filesN(D);%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = H6filesP(D);%VMFG FILE WITH PATH OR DIRECTORY  
                         
                end 
                
                 %IF ONLY ONE FILE IS AVAILABLE[12H] 
        elseif  all([isempty(h0files),isempty(h6files),~isempty(h12files),isempty(h18files)])
                  
                if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H12filesN(D);
                     VMFGfileP = H12filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                               
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H12filesN(D);
                        VMFGfileP = H12filesP(D);
                        
                        beep
                        fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                                                             
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                           
                   end  
               
                elseif  size(Time,2) == 3
                      
                        VMFGfileN = H12filesN(D);%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = H12filesP(D);%VMFG FILE WITH PATH OR DIRECTORY        
                  
                end 
                  
                 %IF ONLY ONE FILE IS AVAILABLE[18H] 
        elseif  all([isempty(h0files),isempty(h6files),isempty(h12files),~isempty(h18files)])
                  
                if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H18filesN(D);
                     VMFGfileP = H18filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 18th(18) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                               
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H18filesN(D);
                        VMFGfileP = H18filesP(D);
                        
                        beep
                        fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the 18th(18) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                                                             
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the 18th(18) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                           
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                          
                   end  
               
                elseif  size(Time,2) == 3
                      
                        VMFGfileN = H18filesN(D);%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = H18filesP(D);%VMFG FILE WITH PATH OR DIRECTORY  
                end 
                  
        end  %//if all([~isempty(h0files),~isempty(h6files),~isempty(h12files),~isempty(h18files)])
                  
        %LOOP OVER VMFG FILESV
        for k = 1:size(VMFGfileP,1)
               
            VMFgridP = VMFGfileP{k,1};%CHANGE FROM CELL ARRAY TO CHARACTER
              
            %EXTRACT ZENITH DELAYS & CREATE AN ARRAY OF II(POSITION EPOCHs) ROWS,k(VMF EPOCHs) COLUMNs    
            %*****CHECK VMF FILE TYPE(VMF1/VMF3)getTROPdelays
            if any([strcmpi(VMFgrid_type,'VMF1'),strfind(VMFgrid_type,'VMF1'),any(VMFgrid_res ==[ 2,2.5])]) 
                   
               VMF_model = 'VMF1'; %flag to indicate VMF file type
                  
               [ZTD(i,k),ZHD(i,k),ZWD(i,k),ah(i,k),aw(i,k)] = readVMF1grid(date(i,1:3),lat(i),lon(i),h(i),VMFgridP,oro);
                                                                  
            elseif any([strcmpi(VMFgrid_type,'VMF3'),strfind(VMFgrid_type,'VMF3'),any(VMFgrid_res ==[ 1,5])]) 
                   
                   VMF_model = 'VMF3';%flag to indicate VMF file type
                      
                   [ZTD(i,k),ZHD(i,k),ZWD(i,k),ah(i,k),aw(i,k)] = readVMF3grid(date(i,1:3),lat(i),lon(i),h(i),VMFgridP,VMFgrid_res,oro);
                          
            end %//if any([strcmpi(VMFgrid_type,'VMF1'),strfind(VMFgrid_type,'VMF1'),any(VMFgrid_res ==[ 2,2.5])]) 
                                     

        end  %//for k = 1:size(VMFGfile,1)       
         
     end   %//for z = 1 : nrow_Rpos

%********CHECK USER DID NOT PROVIDE TIME ELEMENTS
%--------------------------------------------------------------------------
%NOTE:
%     THE TIME INPUT IS SUPPOSE TO BE IN SIX(6) COLUMNS[Yr Mn D Hr Min Sec]
%     HOWEVER, WHEN PROVIDES ONLY THREE(3)[Yr Mn D],THE COLUMNS IN THE ZTDs
%     OUTPUTS WILL CONTAIN THE VMF GRID EPOCHs[00 06 12 18].IN THIS CASE, WE
%     FIND THE MEAN DELAYS OF ALL THE EPOCHS TO OBTAIN THE TOTAL DELAYS FOR
%     THE DAY. I.E ZHD =(DELAY AT 00 + DELAY AT 06 +DELAY AT 12 + DELAY AT 18)/4
%--------------------------------------------------------------------------
if size(Time,2) == 3
    
   %FIND MEAN & REFORMAT COLUMNS INTO DAYS(TIME/DATE)
   zhd_mean  = (round(mean(ZHD,2),5))';
   zwd_mean  = (round(mean(ZWD,2),5))';
   ztd_mean  = (round(mean(ZTD,2),5))';
    
   %ASSIGNMENT
   ZHD = zhd_mean;
   ZWD = zwd_mean;
   ZTD = ztd_mean; 
    
end
     
%--------------------------------------------------------------------------   
  else %IF TIMES ARE NOT SAME/IDENTICAL
%--------------------------------------------------------------------------
      %LOOP OVER STATIONS
      for i = 1 : nrow_Rpos 
    
       %********EXTRACTING TROPOSPHERIC DELAYS
               
       %COMBINE VARIOUS EPOCHs FROM IMPORTED VMF grid FILEs
       for j = 1 : size(date,1)
           Y  = date(j,1); %GET YEAR
           M  = date(j,2); %GET MONTH
           D  = date(j,3); %GET DAY
           H  = date(j,4);%get Hour
         MIN  = date(j,5);%get Minute
         SECs = date(j,6);%get Seconds
           
           if D > size_vmf 
              
              try
                 D = index_t(j);
              catch
                  D = index_t;
              end
             
           end %if D > size_vmf 
           
           %*********CONVERT TIME IN [Hour Min Sec] TO HOURs               
           UT = H + (MIN / 60) + (SECs / 3600);% Time in Hours
          
           %********COMBINE AVAILABLE VMF GRID EPOCHS(00 06 12 18) FILEs FOR EACH AVAILABLE DAY
           %IF ALL FILES ARE AVAILABLE[00H 06H 12H 18H]
           if all([~isempty(h0files),~isempty(h6files),~isempty(h12files),~isempty(h18files)])
            
              if size(Time,2) > 3 
                  
                 %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                 %------------------------------------------------------------------
                 %NOTE: 
                 %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                 %    I.E 00Hr 06Hr 12Hr & 18Hr
                 %------------------------------------------------------------------
                 if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                    VMFGfileN = H0filesN(D);
                    VMFGfileP = H0filesP(D);
            
                 elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                       VMFGfileN = H6filesN(D);
                       VMFGfileP = H6filesP(D);
                
                  
                 elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                       VMFGfileN = H12filesN(D);
                       VMFGfileP = H12filesP(D);
                      
                 elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                       VMFGfileN = H18filesN(D);
                       VMFGfileP = H18filesP(D);
                       
                 end 
               
              elseif size(Time,2) == 3
                     
                    VMFGfileN = [H0filesN(D);H6filesN(D);H12filesN(D);H18filesN(D)];%#ok<*NASGU> %VMFG FILE WITHOUT PATH OR DIRECTORY
                    VMFGfileP = [H0filesP(D);H6filesP(D);H12filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY
           
              end 
                %IF THREE FILES ARE AVAILABLE[00H 06H 12H ]
           elseif all([~isempty(h0files),~isempty(h6files),~isempty(h12files),isempty(h18files)])
                  
                 if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                      
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th Hour epoch. System will use the 12th (12) hour grid file.\nOtherwise you can provide the 18th(12) hour grid files & try again\n\n') 
                         
                   end   
               
                 elseif size(Time,2) == 3
                     
                       VMFGfileN = [H0filesN(D);H6filesN(D);H12filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                       VMFGfileP = [H0filesP(D);H6filesP(D);H12filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY
    
                 end
                 
             %IF THREE FILES ARE AVAILABLE[00H 06H 18H ]
           elseif all([~isempty(h0files),~isempty(h6files),isempty(h12files),~isempty(h18files)])
                  
                 if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th Hour epoch. System will use the 6th (06) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n') 
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                   end  
               
                  elseif size(Time,2) == 3
                   
                        VMFGfileN = [H0filesN(D);H6filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H0filesP(D);H6filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY
           
                  end
                  
                %IF THREE FILES ARE AVAILABLE[00H 12H 18H ]      
           elseif all([~isempty(h0files),isempty(h6files),~isempty(h12files),~isempty(h18files)])
                  
                 if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H0filesN(D);
                        VMFGfileP = H0filesP(D);
                        
                        beep
                        fprintf('\n\nNo VMF grid file at the 6th Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')   
                        
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                   end  
               
                 elseif size(Time,2) == 3
               
                        VMFGfileN = [H0filesN(D);H12filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H0filesP(D);H12filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY
            
                 end
                 
                %IF THREE FILES ARE AVAILABLE[06H 12H 18H ]       
           elseif all([isempty(h0files),~isempty(h6files),~isempty(h12files),~isempty(h18files)])
                   
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H6filesN(D);
                     VMFGfileP = H6filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                   end  
               
                  elseif size(Time,2) == 3
                     
                        VMFGfileN = [H6filesN(D);H12filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H6filesP(D);H12filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY   
               
                  end
                  
                %IF TWO FILES ARE AVAILABLE[00H 06H ]       
           elseif all([~isempty(h0files),~isempty(h6files),isempty(h12files),isempty(h18files)])
                   
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                         
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = [H0filesN(D);H6filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H0filesP(D);H6filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY   
           
                  end
                  
                 %IF TWO FILES ARE AVAILABLE[00H 12H ]
           elseif all([~isempty(h0files),isempty(h6files),~isempty(h12files),isempty(h18files)])
                 
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                         VMFGfileN = H0filesN(D);
                         VMFGfileP = H0filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                         
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                         
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = [H0filesN(D);H12filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H0filesP(D);H12filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY   
             
                  end
                  
                 %IF TWO FILES ARE AVAILABLE[00H 18H ] 
           elseif all([~isempty(h0files),isempty(h6files),isempty(h12files),~isempty(h18files)])
                  
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                         VMFGfileN = H0filesN(D);
                         VMFGfileP = H0filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                                       
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the 18th(18) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                           
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = [H0filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H0filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY 
                  end
                  
                 %IF TWO FILES ARE AVAILABLE[06H 12H ]
           elseif all([isempty(h0files),~isempty(h6files),~isempty(h12files),isempty(h18files)])
                  
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H6filesN(D);
                     VMFGfileP = H6filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                          beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                         
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = [H6filesN(D);H12filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H6filesP(D);H12filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY 
                
                  end
                  
                 %IF TWO FILES ARE AVAILABLE[06H 18H ] 
           elseif all([isempty(h0files),~isempty(h6files),isempty(h12files),~isempty(h18files)])
                  
                 if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H6filesN(D);
                     VMFGfileP = H6filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the 18th(18) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = [H6filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                       VMFGfileP = [H6filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY
                  
                 end
                 
                 %IF TWO FILES ARE AVAILABLE[12H 18H ] 
           elseif all([isempty(h0files),isempty(h6files),~isempty(h12files),~isempty(h18files)])

                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H12filesN(D);
                     VMFGfileP = H12filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H12filesN(D);
                        VMFGfileP = H12filesP(D);
                        
                         beep
                         fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = [H12filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H12filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY 
                  
                  end
                  
                 %IF ONLY ONE FILE IS AVAILABLE[00H] 
           elseif all([~isempty(h0files),isempty(h6files),isempty(h12files),isempty(h18files)])
                  
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H0filesN(D);
                        VMFGfileP = H0filesP(D);
                        
                         beep
                         fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H0filesN(D);
                         VMFGfileP = H0filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H0filesN(D);
                         VMFGfileP = H0filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                           
                   end  
               
                  elseif size(Time,2) == 3 
                      
                         VMFGfileN = H0filesN(D);%VMFG FILE WITHOUT PATH OR DIRECTORY
                         VMFGfileP = H0filesP(D);%VMFG FILE WITH PATH OR DIRECTORY   
                  end
                  
                 %IF ONLY ONE FILE IS AVAILABLE[06H] 
           elseif all([isempty(h0files),~isempty(h6files),isempty(h12files),isempty(h18files)])
                  
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H6filesN(D);
                     VMFGfileP = H6filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                               
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                                                       
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                           
                   end  
               
                  elseif size(Time,2) == 3 
                      
                         VMFGfileN = H6filesN(D);%VMFG FILE WITHOUT PATH OR DIRECTORY
                         VMFGfileP = H6filesP(D);%VMFG FILE WITH PATH OR DIRECTORY  
                         
                  end
                
                 %IF ONLY ONE FILE IS AVAILABLE[12H] 
           elseif all([isempty(h0files),isempty(h6files),~isempty(h12files),isempty(h18files)])
                  
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H12filesN(D);
                     VMFGfileP = H12filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                               
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H12filesN(D);
                        VMFGfileP = H12filesP(D);
                        
                        beep
                        fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                                                             
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                           
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = H12filesN(D);%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = H12filesP(D);%VMFG FILE WITH PATH OR DIRECTORY        
                  
                  end
                  
                 %IF ONLY ONE FILE IS AVAILABLE[18H] 
           elseif all([isempty(h0files),isempty(h6files),isempty(h12files),~isempty(h18files)])
                  
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H18filesN(D);
                     VMFGfileP = H18filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 18th(18) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                               
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H18filesN(D);
                        VMFGfileP = H18filesP(D);
                        
                        beep
                        fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the 18th(18) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                                                             
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the 18th(18) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                           
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                          
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = H18filesN(D);%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = H18filesP(D);%VMFG FILE WITH PATH OR DIRECTORY  
                  end
                  
           end %//if all([~isempty(h0files),~isempty(h6files),~isempty(h12files),~isempty(h18files)])
             
           %LOOP OVER VMFG FILES
           for k = 1:size(VMFGfileP,1)
               
               VMFgridP = VMFGfileP{k,1};%CHANGE FROM CELL ARRAY TO CHARACTER
              
               %EXTRACT ZENITH DELAYS & CREATE AN ARRAY OF II(POSITION EPOCHs) ROWS,k(VMF EPOCHs) COLUMNs    
               %*****CHECK VMF FILE TYPE(VMF1/VMF3)getTROPdelays
               if any([strcmpi(VMFgrid_type,'VMF1'),strfind(VMFgrid_type,'VMF1'),any(VMFgrid_res ==[ 2,2.5])]) 
                   
                  VMF_model = 'VMF1'; %flag to indicate VMF file type
                  
                  [ZTDvmf(j,k),ZHDvmf(j,k),ZWDvmf(j,k),AH(j,k),AW(j,k)] = readVMF1grid(date(j,1:3),lat(i),lon(i),h(i),VMFgridP,oro);
                                                                  
               elseif any([strcmpi(VMFgrid_type,'VMF3'),strfind(VMFgrid_type,'VMF3'),any(VMFgrid_res ==[ 1,5])])
                   
                      VMF_model = 'VMF3';%flag to indicate VMF file type
                      
                     [ZTDvmf(j,k),ZHDvmf(j,k),ZWDvmf(j,k),AH(j,k),AW(j,k)] = readVMF3grid(date(j,1:3),lat(i),lon(i),h(i),VMFgridP,VMFgrid_res,oro); 
                            
               end %//if any([strcmpi(VMFgrid_type,'VMF1'),strfind(VMFgrid_type,'VMF1'),any(VMFgrid_res ==[ 2,2.5])]) 
                   
              %**************CREATE NUMBER OF COLUMNs IN OUTPUT FILEs
              %CREATE A 3D CELL ARRAY FOR THE OUTPUT
              %Create a cell array of No.of observational period/days/date/
              %time(j) as rows,No.of grid epoch[00,06,12,18 UT] files as...
              %columns(k) and the No. of User stations(i) as number of arrays.
              %i.e.Create a j x k x i cell array matrices.
              %------------------------------------------------------------ 
              zhd{j,k,i} = ZHDvmf(j,k);
              zwd{j,k,i} = ZWDvmf(j,k);
              ztd{j,k,i} = ZTDvmf(j,k);
              Ah{j,k,i}  = AH(j,k);
              Aw{j,k,i}  = AW(j,k);

           end  %//for k = 1:size(VMFGfile,1)
           
           
       end    %//for r = 1 : size(date,1)         
         
      end  %//for z = 1 : nrow_Rpos
           
%**********************END OF VMF GRIDDED ZENITH DELAYS EXTRACTION 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%         -----------------------------------------------------------------
%*********EXTRACT ALL STATION DATA & COMBINE AS ONE DATA SET
%         -----------------------------------------------------------------

%LOOP OVER SITE POSITIONs
for q = 1 : nrow_Rpos  
   
    %EXTRACT ALL STATION DATA & COMBINE AS ONE DATA SET
    %ZENITH DELAYS
    zhd_epoch = cell2mat(zhd(:,:,q)); 
    zwd_epoch = cell2mat(zwd(:,:,q)); 
    ztd_epoch = cell2mat(ztd(:,:,q));
    
    %VMF COEFFICIENTS
    ah_epoch = cell2mat(Ah(:,:,q));
    aw_epoch = cell2mat(Aw(:,:,q));
    
    ZPD_epochs = [zhd_epoch  zwd_epoch ztd_epoch ah_epoch aw_epoch];
    
    %*****FIND MEAN OF ALL AVAILABLE EPOCHS[00 06 12 18] FOR DAILY ZPDs
    %I.E. FIND MEAN OF EACH ROW IN ZTD, ZHD & ZWD (MEAN OF VMF EPOCHs ...
    %                                     TO OBTAIN TOTAL DELAY FOR A DAY)
    %NOTE:
    %     THE NUMBER ROWS IN EACH MATRIX REPRESENTS THE NUMBERS OF
    %     DAYS(TIME/DATE) & NUMBER OF COLUMNS ALSO REPRESENTS THE NO. OF
    %     VMF EPOCHS[i.e. 00 06 12 18]
    %----------------------------------------------------------------------
    %*****************ZENITH DELAYS
    %FIND MEAN & REFORMAT COLUMNS INTO DAYS(TIME/DATE)
    zhd_mean  = (round(mean(cell2mat(zhd(:,:,q)),2),5))';
    zwd_mean  = (round(mean(cell2mat(zwd(:,:,q)),2),5))';
    ztd_mean  = (round(mean(cell2mat(ztd(:,:,q)),2),5))';
    
    %ASSIGNMENT
    ZHD(q,:) = zhd_mean;
    ZWD(q,:) = zwd_mean; %#ok<*AGROW>
    ZTD(q,:) = ztd_mean;
    
    %***************VMF COEFFICIENTS
    %FIND MEAN & REFORMAT COLUMNS INTO DAYS(TIME/DATE)
    ah_mean  = (mean(cell2mat(Ah(:,:,q)),2))';
    aw_mean  = (mean(cell2mat(Aw(:,:,q)),2))'; 
    
    %ASSIGNMENT
    ah(q,:) = ah_mean;
    aw(q,:) = aw_mean;
       
end %//for q = 1 : nrow_Rpos  

ZPD_mfc = [ZHD ZWD ZTD ah aw];

  end  %//if isequal(Identical,1)

%--------------------------------------------------------------------------    
else %IF USER POSITION & OBSERVATION TIME ARE UNEQUAL
%--------------------------------------------------------------------------
    %LOOP OVER STATIONS
      for i = 1 : nrow_Rpos 
    
       %********EXTRACTING TROPOSPHERIC DELAYS
               
       %COMBINE VARIOUS EPOCHs FROM IMPORTED VMF grid FILEs
       for j = 1 : size(date,1)
           Y  = date(j,1); %GET YEAR
           M  = date(j,2); %GET MONTH
           D  = date(j,3); %GET DAY
           H  = date(j,4);%get Hour
         MIN  = date(j,5);%get Minute
         SECs = date(j,6);%get Seconds
           
           if D > size_vmf 
              
              try
                 D = index_t(j);
              catch
                  D = index_t;
              end
             
           end %if D > size_vmf 
           
           %*********CONVERT TIME IN [Hour Min Sec] TO HOURs               
           UT = H + (MIN / 60) + (SECs / 3600);% Time in Hours
          
           %********COMBINE AVAILABLE VMF GRID EPOCHS(00 06 12 18) FILEs FOR EACH AVAILABLE DAY
           %IF ALL FILES ARE AVAILABLE[00H 06H 12H 18H]
           if all([~isempty(h0files),~isempty(h6files),~isempty(h12files),~isempty(h18files)])
            
              if size(Time,2) > 3 
                  
                 %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                 %------------------------------------------------------------------
                 %NOTE: 
                 %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                 %    I.E 00Hr 06Hr 12Hr & 18Hr
                 %------------------------------------------------------------------
                 if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                    VMFGfileN = H0filesN(D);
                    VMFGfileP = H0filesP(D);
            
                 elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                       VMFGfileN = H6filesN(D);
                       VMFGfileP = H6filesP(D);
                
                  
                 elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                       VMFGfileN = H12filesN(D);
                       VMFGfileP = H12filesP(D);
                      
                 elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                       VMFGfileN = H18filesN(D);
                       VMFGfileP = H18filesP(D);
                       
                 end 
               
              elseif size(Time,2) == 3
                     
                    VMFGfileN = [H0filesN(D);H6filesN(D);H12filesN(D);H18filesN(D)];%#ok<*NASGU> %VMFG FILE WITHOUT PATH OR DIRECTORY
                    VMFGfileP = [H0filesP(D);H6filesP(D);H12filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY
           
              end 
                %IF THREE FILES ARE AVAILABLE[00H 06H 12H ]
           elseif all([~isempty(h0files),~isempty(h6files),~isempty(h12files),isempty(h18files)])
                  
                 if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                      
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th Hour epoch. System will use the 12th (12) hour grid file.\nOtherwise you can provide the 18th(12) hour grid files & try again\n\n') 
                         
                   end   
               
                 elseif size(Time,2) == 3
                     
                       VMFGfileN = [H0filesN(D);H6filesN(D);H12filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                       VMFGfileP = [H0filesP(D);H6filesP(D);H12filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY
    
                 end
                 
             %IF THREE FILES ARE AVAILABLE[00H 06H 18H ]
           elseif all([~isempty(h0files),~isempty(h6files),isempty(h12files),~isempty(h18files)])
                  
                 if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th Hour epoch. System will use the 6th (06) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n') 
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                   end  
               
                  elseif size(Time,2) == 3
                   
                        VMFGfileN = [H0filesN(D);H6filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H0filesP(D);H6filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY
           
                  end
                  
                %IF THREE FILES ARE AVAILABLE[00H 12H 18H ]      
           elseif all([~isempty(h0files),isempty(h6files),~isempty(h12files),~isempty(h18files)])
                  
                 if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H0filesN(D);
                        VMFGfileP = H0filesP(D);
                        
                        beep
                        fprintf('\n\nNo VMF grid file at the 6th Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')   
                        
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                   end  
               
                 elseif size(Time,2) == 3
               
                        VMFGfileN = [H0filesN(D);H12filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H0filesP(D);H12filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY
            
                 end
                 
                %IF THREE FILES ARE AVAILABLE[06H 12H 18H ]       
           elseif all([isempty(h0files),~isempty(h6files),~isempty(h12files),~isempty(h18files)])
                   
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H6filesN(D);
                     VMFGfileP = H6filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                   end  
               
                  elseif size(Time,2) == 3
                     
                        VMFGfileN = [H6filesN(D);H12filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H6filesP(D);H12filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY   
               
                  end
                  
                %IF TWO FILES ARE AVAILABLE[00H 06H ]       
           elseif all([~isempty(h0files),~isempty(h6files),isempty(h12files),isempty(h18files)])
                   
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                         
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = [H0filesN(D);H6filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H0filesP(D);H6filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY   
           
                  end
                  
                 %IF TWO FILES ARE AVAILABLE[00H 12H ]
           elseif all([~isempty(h0files),isempty(h6files),~isempty(h12files),isempty(h18files)])
                 
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                         VMFGfileN = H0filesN(D);
                         VMFGfileP = H0filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                         
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                         
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = [H0filesN(D);H12filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H0filesP(D);H12filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY   
             
                  end
                  
                 %IF TWO FILES ARE AVAILABLE[00H 18H ] 
           elseif all([~isempty(h0files),isempty(h6files),isempty(h12files),~isempty(h18files)])
                  
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                         VMFGfileN = H0filesN(D);
                         VMFGfileP = H0filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                                       
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the 18th(18) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                           
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = [H0filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H0filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY 
                  end
                  
                 %IF TWO FILES ARE AVAILABLE[06H 12H ]
           elseif all([isempty(h0files),~isempty(h6files),~isempty(h12files),isempty(h18files)])
                  
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H6filesN(D);
                     VMFGfileP = H6filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                          beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                         
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = [H6filesN(D);H12filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H6filesP(D);H12filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY 
                
                  end
                  
                 %IF TWO FILES ARE AVAILABLE[06H 18H ] 
           elseif all([isempty(h0files),~isempty(h6files),isempty(h12files),~isempty(h18files)])
                  
                 if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H6filesN(D);
                     VMFGfileP = H6filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the 18th(18) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = [H6filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                       VMFGfileP = [H6filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY
                  
                 end
                 
                 %IF TWO FILES ARE AVAILABLE[12H 18H ] 
           elseif all([isempty(h0files),isempty(h6files),~isempty(h12files),~isempty(h18files)])

                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H12filesN(D);
                     VMFGfileP = H12filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H12filesN(D);
                        VMFGfileP = H12filesP(D);
                        
                         beep
                         fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = [H12filesN(D);H18filesN(D)];%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = [H12filesP(D);H18filesP(D)];%VMFG FILE WITH PATH OR DIRECTORY 
                  
                  end
                  
                 %IF ONLY ONE FILE IS AVAILABLE[00H] 
           elseif all([~isempty(h0files),isempty(h6files),isempty(h12files),isempty(h18files)])
                  
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H0filesN(D);
                     VMFGfileP = H0filesP(D);
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H0filesN(D);
                        VMFGfileP = H0filesP(D);
                        
                         beep
                         fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                                      
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H0filesN(D);
                         VMFGfileP = H0filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H0filesN(D);
                         VMFGfileP = H0filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the zero(00) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                           
                   end  
               
                  elseif size(Time,2) == 3 
                      
                         VMFGfileN = H0filesN(D);%VMFG FILE WITHOUT PATH OR DIRECTORY
                         VMFGfileP = H0filesP(D);%VMFG FILE WITH PATH OR DIRECTORY   
                  end
                  
                 %IF ONLY ONE FILE IS AVAILABLE[06H] 
           elseif all([isempty(h0files),~isempty(h6files),isempty(h12files),isempty(h18files)])
                  
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H6filesN(D);
                     VMFGfileP = H6filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                               
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H6filesN(D);
                        VMFGfileP = H6filesP(D);
                                                       
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H6filesN(D);
                         VMFGfileP = H6filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the 6th(06) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                           
                   end  
               
                  elseif size(Time,2) == 3 
                      
                         VMFGfileN = H6filesN(D);%VMFG FILE WITHOUT PATH OR DIRECTORY
                         VMFGfileP = H6filesP(D);%VMFG FILE WITH PATH OR DIRECTORY  
                         
                  end
                
                 %IF ONLY ONE FILE IS AVAILABLE[12H] 
           elseif all([isempty(h0files),isempty(h6files),~isempty(h12files),isempty(h18files)])
                  
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H12filesN(D);
                     VMFGfileP = H12filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                               
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H12filesN(D);
                        VMFGfileP = H12filesP(D);
                        
                        beep
                        fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                                                             
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H12filesN(D);
                         VMFGfileP = H12filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 18th(18) Hour epoch. System will use the 12th(12) hour grid file.\nOtherwise you can provide the 18th(18) hour grid files & try again\n\n')     
                           
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = H12filesN(D);%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = H12filesP(D);%VMFG FILE WITH PATH OR DIRECTORY        
                  
                  end
                  
                 %IF ONLY ONE FILE IS AVAILABLE[18H] 
           elseif all([isempty(h0files),isempty(h6files),isempty(h12files),~isempty(h18files)])
                  
                  if size(Time,2) > 3 
                  
                   %********FIND THE CORRECT FILES BASED ON THE OBSERVATION TIME(UTC)
                   %------------------------------------------------------------------
                   %NOTE: 
                   %    VMF GRID FILES ARE GIVEN AT FOUR EPOCHS WITH 6 HOUR RESOLUTION
                   %    I.E 00Hr 06Hr 12Hr & 18Hr
                   %------------------------------------------------------------------
                   if (UT >= 0 & UT < 6)%IF UTC TIME RANGE IS BETWEEN(0-6)=>00H FILES

                     VMFGfileN = H18filesN(D);
                     VMFGfileP = H18filesP(D);
                     
                     beep
                     fprintf('\n\nNo VMF grid file at the zero(00) Hour epoch. System will use the 18th(18) hour grid file.\nOtherwise you can provide the zero(00) hour grid files & try again\n\n')     
                               
                     
                   elseif (UT >= 6 & UT < 12)%IF UTC TIME RANGE IS BETWEEN(6-12)=>06H FILES  
                   
                        VMFGfileN = H18filesN(D);
                        VMFGfileP = H18filesP(D);
                        
                        beep
                        fprintf('\n\nNo VMF grid file at the 6th(06) Hour epoch. System will use the 18th(18) hour grid file.\nOtherwise you can provide the 6th(06) hour grid files & try again\n\n')     
                                                             
                   elseif (UT >= 12 & UT < 18) %IF UTC TIME RANGE IS BETWEEN(12-18) =>12H FILES  
                      
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                         
                         beep
                         fprintf('\n\nNo VMF grid file at the 12th(12) Hour epoch. System will use the 18th(18) hour grid file.\nOtherwise you can provide the 12th(12) hour grid files & try again\n\n')     
                           
                   elseif (UT >= 18 & UT <= 24) %IF UTC TIME RANGE IS BETWEEN(18-24) =>18H FILES 
                     
                         VMFGfileN = H18filesN(D);
                         VMFGfileP = H18filesP(D);
                          
                   end  
               
                  elseif size(Time,2) == 3
                      
                        VMFGfileN = H18filesN(D);%VMFG FILE WITHOUT PATH OR DIRECTORY
                        VMFGfileP = H18filesP(D);%VMFG FILE WITH PATH OR DIRECTORY  
                  end
                  
           end %//if all([~isempty(h0files),~isempty(h6files),~isempty(h12files),~isempty(h18files)])
                
           %LOOP OVER VMFG FILES
           for k = 1:size(VMFGfileP,1)
               
               VMFgridP = VMFGfileP{k,1};%CHANGE FROM CELL ARRAY TO CHARACTER
              
               %EXTRACT ZENITH DELAYS & CREATE AN ARRAY OF II(POSITION EPOCHs) ROWS,k(VMF EPOCHs) COLUMNs    
               %*****CHECK VMF FILE TYPE(VMF1/VMF3)getTROPdelays
               if any([strcmpi(VMFgrid_type,'VMF1'),strfind(VMFgrid_type,'VMF1'),any(VMFgrid_res ==[ 2,2.5])]) 
                   
                  VMF_model = 'VMF1'; %flag to indicate VMF file type
                  
                  [ZTDvmf(j,k),ZHDvmf(j,k),ZWDvmf(j,k),AH(j,k),AW(j,k)] = readVMF1grid(date(j,1:3),lat(i),lon(i),h(i),VMFgridP,oro);
                                                                  
               elseif any([strcmpi(VMFgrid_type,'VMF3'),strfind(VMFgrid_type,'VMF3'),any(VMFgrid_res ==[ 1,5])])
                   
                      VMF_model = 'VMF3';%flag to indicate VMF file type
                      
                     [ZTDvmf(j,k),ZHDvmf(j,k),ZWDvmf(j,k),AH(j,k),AW(j,k)] = readVMF3grid(date(j,1:3),lat(i),lon(i),h(i),VMFgridP,VMFgrid_res,oro); 
                            
               end %//if any([strcmpi(VMFgrid_type,'VMF1'),strfind(VMFgrid_type,'VMF1'),any(VMFgrid_res ==[ 2,2.5])]) 
                   
             
              %**************CREATE NUMBER OF COLUMNs IN OUTPUT FILEs
              %CREATE A 3D CELL ARRAY FOR THE OUTPUT
              %Create a cell array of No.of observational period/days/date/
              %time(j) as rows,No.of grid epoch[00,06,12,18 UT] files as...
              %columns(k) and the No. of User stations(i) as number of arrays.
              %i.e.Create a j x k x i cell array matrices.
              %------------------------------------------------------------
              zhd{j,k,i} = ZHDvmf(j,k);
              zwd{j,k,i} = ZWDvmf(j,k);
              ztd{j,k,i} = ZTDvmf(j,k);
              Ah{j,k,i}  = AH(j,k);
              Aw{j,k,i}  = AW(j,k);

           end  %//for k = 1:size(VMFGfile,1)
           
           
       end    %//for r = 1 : size(date,1)         
         
      end  %//for z = 1 : nrow_Rpos
           
%**********************END OF VMF GRIDDED ZENITH DELAYS EXTRACTION 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%         -----------------------------------------------------------------
%*********EXTRACT ALL STATION DATA & COMBINE AS ONE DATA SET
%         -----------------------------------------------------------------

%LOOP OVER SITE POSITIONs
for q = 1 : nrow_Rpos  
   
    %ZENITH DELAYS
    zhd_epoch = cell2mat(zhd(:,:,q)); 
    zwd_epoch = cell2mat(zwd(:,:,q)); 
    ztd_epoch = cell2mat(ztd(:,:,q));
    
    %VMF COEFFICIENTS
    ah_epoch = cell2mat(Ah(:,:,q));
    aw_epoch = cell2mat(Aw(:,:,q));
    
    ZPD_epochs = [zhd_epoch  zwd_epoch ztd_epoch ah_epoch aw_epoch];
    
    %*****FIND MEAN OF ALL AVAILABLE EPOCHS[00 06 12 18] FOR DAILY ZPDs
    %I.E. FIND MEAN OF EACH ROW IN ZTD, ZHD & ZWD (MEAN OF VMF EPOCHs ...
    %                                     TO OBTAIN TOTAL DELAY FOR A DAY)
    %NOTE:
    %     THE NUMBER ROWS IN EACH MATRIX REPRESENTS THE NUMBERS OF
    %     DAYS(TIME/DATE) & NUMBER OF COLUMNS ALSO REPRESENTS THE NO. OF
    %     VMF EPOCHS[i.e. 00 06 12 18]
    %----------------------------------------------------------------------
    %*****************ZENITH DELAYS
    %FIND MEAN & REFORMAT COLUMNS INTO DAYS(TIME/DATE)
    zhd_mean  = (round(mean(cell2mat(zhd(:,:,q)),2),5))';
    zwd_mean  = (round(mean(cell2mat(zwd(:,:,q)),2),5))';
    ztd_mean  = (round(mean(cell2mat(ztd(:,:,q)),2),5))';
    
    %ASSIGNMENT
    ZHD(q,:) = zhd_mean;
    ZWD(q,:) = zwd_mean; %#ok<*AGROW>
    ZTD(q,:) = ztd_mean;
    
    %***************VMF COEFFICIENTS
    %FIND MEAN & REFORMAT COLUMNS INTO DAYS(TIME/DATE)
    ah_mean  = (mean(cell2mat(Ah(:,:,q)),2))';
    aw_mean  = (mean(cell2mat(Aw(:,:,q)),2))'; 
    
    %ASSIGNMENT
    ah(q,:) = ah_mean;
    aw(q,:) = aw_mean;
       
end %//for q = 1 : nrow_Rpos  

ZPD_mfc = [ZHD ZWD ZTD ah aw]; 

end %//if isequal(nrow_time,nrow_pos)

%%%%%%%%%%%%%%%%%%%%%%%%END OF ZPD EXTRACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%********SUB-ROUTINE TO READ  AND EXTRACT ZENITH TROPOSPHERIC DELAYS

%A.SUBROUTINE TO READ VMF1 GRID FILE & EXTRACT COEFFICIENTs(ZHD,ZWD,ZTD,ah,aw)
function [ZTD,ZHD,ZWD,ah,aw] = readVMF1grid(UTCtime,Lat,Lon,hgt,VMF1_grid,...
                                                            orography_file)                                                     
%**************************************************************************
%DESCRIPTION:
%***********This subroutine determines the Hydrostatic and Wet Mapping ... * 
%           Function(MF)Coefficients ah and aw,as well as the Zenith Delays*
%           (ZHD,ZWD)from the gridded VMF1 files, as available from:       *
%           http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/VMFG/ NOW replace with* 
%           http://vmf.geo.tuwien.ac.at/trop_products/GRID/2.5x2/VMF1/ for * 
%           specific sites near the  Earth surface. The VMF1 grid file is  *                                                
%           based on a 2.0 × 2.5 degree external grid file given for four  * 
%           daily epochs (0h, 6h, 12h, 18h UT).consist two MF coefficients * 
%           (ah,aw)as well as the Hydrostatic and Wet portionsof ZPD(zh,zw)* 

%***********Furthermore, since the ZPDs(zh,zw) correspond to mean grid     * 
%           heights, a grid file with these mean ellipsoidal heights       *
%           (orography_ell) is also required.The orography_ell file can be * 
%           downloaded from:                                               *
%           [http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/orography_ell] NOW  *  
%           replace with http://vmf.geo.tuwien.ac.at/station_coord_files/  *
%           Example of vmf grid file for January 2018 at(0h,6h,12h,18h UT) *
%           include:[VMFG_20180101.H00,VMFG_20180101.H06,VMFG_20180101.H12 *
%           VMFG_20180101.H18]                                             *
%           ***************************************************************
%           ***************************************************************
%********** On the temporal scale, the values from the two surrounding NWM * 
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
%      [ZTD,ZHD,ZWD,ah,aw] = Readvmf_grid(UTCtime,lat,lon,hell,VMF1_grid...
%                                                          ,orography_file)
%INPUTs:
%1.     UTCtime:.........UTC time in [Year,Month,Day,Hour,Minute,Seconds]  *
%2.     Lat: ............station ellipsoidal latitude in [degrees]         *
%3.     Lon: ............station ellipsoidal longitude in [degrees]        *
%4.     hgt: ............station ellipsoidal height in [meters]            *
%5.     VMF1grid_file:...VMF grid file eg:'VMFG_20180101.H00'              *
%6.     orography_file:...ellipsoidal orography. eg:'orography_ell'        *

%OUTPUTs:                                                                  *
%        ZHD ............... Zenith Hydrostatic Delay, valid at hell       *
%        ZWD ............... Zenith Wet Delay, valid at hell               *
%        ZTD ................Zenith Total Delay, valid at hell             *
%        ah ............... hydrostatic mapping coefficient, valid at hell *
%        aw ............... wet mapping coefficient, valid at hell         *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% REFERENCE:                                                               *
%o          Kouba, J. (2008), Implementation and testing of the gridded    *
%           Vienna Mapping Function 1 (VMF1). J. Geodesy, Vol. 82:193-205, * 
%           DOI: 10.1007/s00190-007-0170-0                                 *
%o          Niell, A.E. (1996), Global mapping functions for the atmosphere 
%           delay at 310 radio wavelengths. J. Geophys. Res.,101,3227-3246 *

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
   
%****GET VARIOUS TIME INPUTs 
%1.******DATE[YEAR MONTH DAY]
Yr  = UTCtime(:,1);%get Yr
Mn  = UTCtime(:,2);%get Month
Day = UTCtime(:,3);%get Day

%2.******TIME[HOUR MINUTE SECONDs]
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

%********ASSIGNMENTs
%1.STATION COORDs
latDeg = Lat; %assign Lat in degrees to latDeg
lonDeg = Lon; %assign Lon in degrees to lonDeg

%2.VMF grid FILE
VMFgrid=VMF1_grid;

%***CONVERT Lat and Lon TO RADIAN 
lat = (Lat./180).*pi; %[radian]
lon = (Lon./180).*pi ;%[radian]                                                     
 
%********COMPUTE MODIFIED JULIAN DATE(MJD)
%Call the "utc2JulianDay_DoY.m" function
[~, MJD ,~]=utc2JulianDay_DoY(Yr,Mn,Day,H,MIN,SECs);
   

%****FIND THE TWO SURROUNDING EPOCHs
if mod(MJD,0.25)==0
    MJD_all = MJD;
else
    MJD_int = floor(MJD*4)/4 : 0.25 : ceil(MJD*4)/4;
    MJD_all = [MJD MJD_int];
end

MJD_all(H==24)=MJD_all(H==24)+1;

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
fileID = fopen(VMFgrid);%GET FILE ID
                        %Only read data up to the maximum index in order to save time
VMF1_data_all = textscan(fileID,'%f%f%f%f%f%f',max(index),'CommentStyle','!','CollectOutput',1);  
                      
                        %Reduce to the indices of the surrounding grid points
vec =VMF1_data_all{1} ;                       
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
    VMF1_data_int_h0(iv_ind,iv_line) = VMF1_data{1}(iv_ind,iv_line) + (VMF1_data{2}(iv_ind,iv_line)-VMF1_data{1}(iv_ind,iv_line))*(MJD-MJD_int(1))/(MJD_int(2)-MJD_int(1));   % the appendix 'h0' means that the values are valid at zero height
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
VMF1_data_int_h0(iv_ind,7) = (VMF1_data_int_h0(iv_ind,5)/0.0022768) .* (1-0.00266*cos(2*lat)-0.28*10^-6*orography_ell(index));% (1) convert the hydrostatic zenith delay at grid height to the respective pressure value   

%(2)Transfer the pressure from each grid height to site height using
%   Berg(1948) pressure formular
VMF1_data_int_h1(iv_ind,7) = VMF1_data_int_h0(iv_ind,7).*(1-0.0000226.*(hgt-orography_ell(index))).^5.225;% (2) lift the pressure each from grid height to site height   

%(3)NOW Obtain ZHD at station/site height
VMF1_data_int_h1(iv_ind,5) = 0.0022768*VMF1_data_int_h1(iv_ind,7) / (1-0.00266*cos(2*lat)-0.28*10^-6*hgt);   % (3) convert the lifted pressure to zhd again (as proposed by Kouba, 2008)

%(b)*********ZENITH WET DELAY(ZWD)
%Use a simple exponential decay approximation function
VMF1_data_int_h1(iv_ind,6) = VMF1_data_int_h0(iv_ind,6) .* exp(-(hgt-orography_ell(index))/2000);

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
      
end %//if length(unique(index)) == 1 

%*******ASSIGNMENT
ZHD=zhd;%Zenith Hydrostatic Delay
ZWD=zwd;%Zenith Wet Delay
ah=hcoe;%Hydrostatic Coefficient
aw=wcoe;%Wet Coefficient

%COMPUTE ZENITH TOTAL DELAY
ZTD=ZHD + ZWD;

%%%%%%%%%%%%%%%%%%%%%%%%%END OF readVMF1grid.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-     
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%B.SUBROUTINE TO READ VMF3 GRID FILE & EXTRACT COEFFICIENTs(zhd,zwd,ztd,ah,aw)
function [ZTD,ZHD,ZWD,ah,aw] = readVMF3grid(UTCtime,Lat,Lon,hgt,VMF3_grid,...
                                                            gridRESOLUTION,orography_file)                                                     
%**************************************************************************
%DESCRIPTION:
%***********This subroutine determines the Hydrostatic and Wet Mapping ... * 
%           Function(MF)Coefficients ah and aw,as well as the Zenith Delays*
%           (ZHD,ZWD)from the gridded VMF3 files, as available from:       *
%           http://vmf.geo.tuwien.ac.at/trop_products/GRID/ for specific   * 
%           sites near the  Earth surface. The VMF3 grid file is based on  *                                                
%           a 1° × 1° and 5° × 5° external grid files given for four daily* 
%           epochs (0h, 6h, 12h, 18h UT).Consist two MF coefficients(ah,aw)* 
%           as well as the Hydrostatic and Wet portions of ZPD(zh,zw)      * 

%***********Furthermore, since the ZPDs(zh,zw) correspond to mean grid     * 
%           heights, a grid file with these mean ellipsoidal heights       *
%           (orography_ell_1x1 or orography_ell_5x5 for 1° and 5° grid ... *
%           resolutions respectively) is also required.The orography_ell   * 
%           file can as well be downloaded from:                           * 
%           http://vmf.geo.tuwien.ac.at/station_coord_files/               *
%           Example of vmf3 grid file for January 2018 at(0h,6h,12h,18h UT)*
%           include:[VMF3_20180101.H00,VMF3_20180101.H06,VMF3_20180101.H12 *
%           VMF3_20180101.H18]                                             *
%           ***************************************************************
%           ***************************************************************
%********** On the temporal scale, the values from the two surrounding NWM * 
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
%      [ZTD,ZHD,ZWD,ah,aw] = Readvmf_grid(UTCtime,lat,lon,hell,VMF1_grid...
%                                             gridRESOLUTION,orography_file)
%INPUTs:
%1.     UTCtime:.........UTC time in [Year,Month,Day,Hour,Minute,Seconds]  *
%2.     Lat: ............station ellipsoidal latitude in [degrees]         *
%3.     Lon: ............station ellipsoidal longitude in [degrees]        *
%4.     hgt: ............station ellipsoidal height in [meters]            *
%5.     VMF1_grid:.......VMF grid file eg:'VMFG_20180101.H00'              *
%6.     gridRESOLUTION....Grid resolution(°) (possible: 1 or 5)            *
%7.     orography_file:..ellipsoidal orography. eg:'orography_ell'         *

%OUTPUTs:                                                                  *
%        ZHD ............... Zenith Hydrostatic Delay, valid at hell       *
%        ZWD ............... Zenith Wet Delay, valid at hell               *
%        ZTD ................Zenith Total Delay, valid at hell             *
%        ah ............... hydrostatic mapping coefficient, valid at hell *
%        aw ............... wet mapping coefficient, valid at hell         *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% REFERENCE:                                                               *
%o          Reference for conversion of zenith delays:
%           Kouba, J. (2008), Implementation and testing of the gridded    *
%           Vienna Mapping Function 1 (VMF1). J. Geodesy, Vol. 82:193-205, * 
%           DOI: 10.1007/s00190-007-0170-0                                 *

%o          Reference for conversion of mapping functions:
%           Niell, A.E. (1996),Global mapping functions for the atmosphere * 
%           delay at 310 radio wavelengths. J. Geophys. Res.,101,3227-3246 *

%o          Reference for VMF3:
%           Landskron, D. & Böhm, J. J Geod (2017). https://doi.org/10.1007
%                                                        /s00190-017-1066-2
%Original codes by Daniel Landskron (2017/06/28)                           *
%Modified by: Samuel Osah,Msc Geomatic Engineering ,2016                   *    
%             Email: osahsamuel@yahoo.ca                                   *
%             Tel:+233 (0)246137410/+233 (0)509438484                      *        
%==========================================================================+
%**************************************************************************+
%**************************************************************************+

%*************************STEPS TO READ VMF1 GRID FILE
%(1)CHECK IF OROGRAPHY FILE IS PROVIDED
 switch nargin
     case 7
         orogfile=orography_file;
     case  6
         orogfile=[];
 end 
   
%(2)****GET VARIOUS TIME INPUTs 

%(2.1)******DATE[YEAR MONTH DAY]
Yr  = UTCtime(:,1);%get Hour
Mn  = UTCtime(:,2);%get Month
Day = UTCtime(:,3);%get Day

%(2.2)******TIME[HOUR MINUTE SECONDs]
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

%(3)**ASSIGNMENT
%(3.1)STATION COORDs
latDeg = Lat; %assign Lat in degrees to latDeg
lonDeg = Lon; %assign Lon in degrees to lonDeg

%(3.2)VMF1 grid FILE
vmf_grid=VMF3_grid;

%(3.3)GRID RESOLUTION
grid_res = gridRESOLUTION;

%(4)***CONVERT Lat and Lon TO RADIAN 
lat = (Lat./180).*pi; %[radian]
lon = (Lon./180).*pi ;%[radian] 

%(5)********COMPUTE MODIFIED JULIAN DATE(MJD)
%Call the "utc2JulianDay_DoY.m" function
[~, MJD ,~]=utc2JulianDay_DoY(Yr,Mn,Day,H,MIN,SECs);

%(6)****FIND THE TWO SURROUNDING EPOCHs
if mod(MJD,0.25)==0
    MJD_all = MJD;
else
    MJD_int = floor(MJD*4)/4 : 0.25 : ceil(MJD*4)/4;
    MJD_all = [MJD MJD_int];
end

MJD_all(H==24)=MJD_all(H==24)+1;

%epoch = (MJD_all-floor(MJD_all))*24;

%(7)*******ONLY +VE LONGITUDE IN [degrees] is REQUIRED
if lonDeg < 0
    lon = lon + 2*pi;
    lonDeg = (lonDeg + 360);
end

%(8)****READ OROGRAPHY_ell FILE
if ~isempty(orogfile)
    
  flag_orog=0; %flag to indicate  orogfile is not empty
  
  try  
     fileID = fopen(orogfile);%Open Orography file
     
     orography_ell = textscan(fileID,'%f');
     orography_ell = cell2mat(orography_ell);
     
     fclose(fileID);%Close fileID
     
  catch
       orography_ell = importOROGRAPHY(orogfile);
  end %//try
  
else
     flag_orog=1;%flag to indicate  orogfile does not exist
end     
    
%(9)************* FIND INDICES OF FOUR(4) SURROUNDING GRID POINTs

%(9.1) FIND COORDINATES (LAT,LON) OF THE SURROUNDING GRID POINTS
lat_all = 90-grid_res/2 : -grid_res : -90;
lon_all = 0+grid_res/2 : grid_res : 360;

%9.2)FIND THE TWO(20 CLOSEST LATITUDES 
lat_temp = latDeg-lat_all;
[~,ind_lat_int(1)] = min(abs(lat_temp));
ind_lat_int(2) = ind_lat_int(1)-sign(lat_temp(ind_lat_int(1)));

%9.2)FIND THE TWO(20 CLOSEST LONGITUDES 
lon_temp = lonDeg-lon_all;
[~,ind_lon_int(1)] = min(abs(lon_temp));
ind_lon_int(2) = ind_lon_int(1)+sign(lon_temp(ind_lon_int(1)));

%(9.3)CORRECT THE INDICES OUT OF RANGE 
for i_ind = 1:2
    if ind_lat_int(i_ind)>length(lat_all); ind_lat_int(i_ind) = length(lat_all);                    end
    if ind_lat_int(i_ind)<1;               ind_lat_int(i_ind) = 1;                                  end
    if ind_lon_int(i_ind)>length(lon_all); ind_lon_int(i_ind) = ind_lon_int(i_ind)-length(lon_all); end
    if ind_lon_int(i_ind)<1;               ind_lon_int(i_ind) = ind_lon_int(i_ind)+length(lon_all); end
end

%(9.4)DEFINE INDICES
index(1) = (ind_lat_int(1)-1)*length(lon_all)+ind_lon_int(1);
index(2) = (ind_lat_int(1)-1)*length(lon_all)+ind_lon_int(2);
index(3) = (ind_lat_int(2)-1)*length(lon_all)+ind_lon_int(1);
index(4) = (ind_lat_int(2)-1)*length(lon_all)+ind_lon_int(2);

%(10) READ THE CORRECT DATA AND PERFORM A LINEAR TIME INTERPOLATION FROM 
%     THE SURROUNDING TWO EPOCHS AND READ IN WITH TEXTSCAN, BUT ONLY UP TO
%     MAXIMUM INDEX, EVERYTHING BEFORE WILL BE TREATED AS HEADERLINES

%******OPEN THE VMFG FILE
FID = fopen(vmf_grid);%GET FILE ID
                        %Only read data up to the maximum index in order to save time
VMF3_data_all = textscan(FID,'%f%f%f%f%f%f',max(index),'CommentStyle','!','CollectOutput',1);  
                      
                        %Reduce to the indices of the surrounding grid points                       
VMF3_data = cellfun(@(c) c(index,:),VMF3_data_all,'UniformOutput',false);

%CLOSE THE OPENED FILE
fclose(FID);

%******INITIALIZE GRID OUTPUT
%NOTE: The vmf grid file has 6 columns i.e.[(lat lon ah aw zhd zwd)]
VMF3_data_int_h0 = zeros(4,7);%Grid values @ grid height (7TH column for grid pressure values)
VMF3_data_int_h1 = zeros(4,9);%Grid values reduced to site/station height
                              %(7TH column for site pressure values)
%DO THE LINEAR TIME INTERPOLATION FOR EACH ARGUMENT; THE RESULTS ARE THE 
%VMF1 VALUES FOR THE SURROUNDING GRID POINTS AT THE TIME OF THE MEASUREMENT
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
iv_ind = 1:4;
if length(MJD_all)==1 %if the observation epoch coincides with an NWM epoch
    VMF3_data_int_h0(iv_ind,1:6) = VMF3_data{1}(iv_ind,1:6);
    
else %otherwise,perform the linear interpolation
    iv_line = 1:6;
    VMF3_data_int_h0(iv_ind,iv_line) = VMF3_data{1}(iv_ind,iv_line) + (VMF3_data{2}(iv_ind,iv_line)-VMF3_data{1}(iv_ind,iv_line))*(MJD-MJD_int(1))/(MJD_int(2)-MJD_int(1));   % the appendix 'h0' means that the values are valid at zero height
end

%ASSIGNING THE FIRST FOUR(4) COLUMNs(They are equal)
VMF3_data_int_h1(:,1:4) = VMF3_data_int_h0(:,1:4);


%(10.1) BRING MFH, MFW, ZHD AND ZWD OF THE SURROUNDING GRID POINTS TO THE... 
%                  RESPECTIVE HEIGHT OF THE LOCATION
%--------------------------------------------------------------------------
if isequal(flag_orog,0) %if orography file is not empty
    
%(A)*********ZENITH HYDROSTATIC DELAY(ZHD)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%(1)Convert the hydrostatic zenith delay at grid height to the respective
%   pressure value using station latitude and height
%NOTE:
%    to be exact, the latitudes of the respective grid points would have to 
%    be used instead of the latitude of the station (lat). However,the loss
%    of accuracy is only in the sub-micrometer range.
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%USE DAVIS ET AL(1985) FORMULAR TO COMPUTE PRESSURE
VMF3_data_int_h0(iv_ind,7) = (VMF3_data_int_h0(iv_ind,5)/0.0022768) .* (1-0.00266*cos(2*lat)-0.28*10^-6*orography_ell(index));% (1) convert the hydrostatic zenith delay at grid height to the respective pressure value   

%(2)Transfer the pressure from each grid height to site height using
%   Berg(1948) pressure formular
VMF3_data_int_h1(iv_ind,7) = VMF3_data_int_h0(iv_ind,7).*(1-0.0000226.*(hgt-orography_ell(index))).^5.225;% (2) lift the pressure each from grid height to site height   

%(3)NOW Obtain ZHD at station/site height
VMF3_data_int_h1(iv_ind,5) = 0.0022768*VMF3_data_int_h1(iv_ind,7) / (1-0.00266*cos(2*lat)-0.28*10^-6*hgt);   % (3) convert the lifted pressure to zhd again (as proposed by Kouba, 2008)

%(B)*********ZENITH WET DELAY(ZWD)
%Use a simple exponential decay approximation function
VMF3_data_int_h1(iv_ind,6) = VMF3_data_int_h0(iv_ind,6) .* exp(-(hgt-orography_ell(index))/2000);

else %if orography file is empty
    %(a)*********ZENITH HYDROSTATIC DELAY(ZHD)
    VMF3_data_int_h1(iv_ind,5) =VMF3_data_int_h0(iv_ind,5);
    
    %(b)*********ZENITH WET DELAY(ZWD)
    VMF3_data_int_h1(iv_ind,6) =VMF3_data_int_h0(iv_ind,6);
   
end %if isequal(flag_orog,0)

%(C)*********VMF COEFFICIENTs [ah, aw]
hcoe = VMF3_data_int_h0(iv_ind,3);%Hydrostatic coefficient
wcoe = VMF3_data_int_h0(iv_ind,4);%Wet coefficient

VMF3_data_int_h1(iv_ind,8) = hcoe;
VMF3_data_int_h1(iv_ind,9) = wcoe;

%(10.2)******* PERFORM THE BILINEAR INTERPOLATION

if length(unique(index)) == 1   % if the point is directly on a grid point
    
    zhd  = VMF3_data_int_h1(1,5);
    zwd  = VMF3_data_int_h1(1,6);
    hcoe = VMF3_data_int_h1(1,8);
    wcoe = VMF3_data_int_h1(1,9);

else
    %BILINEAR INTERPOLATION (INTERPRETED AS TWO 1D LINEAR INTERPOLATIONS 
    %FOR LAT AND LON, BUT PROGRAMMED WITHOUT SUBFUNCTIONS) 
    
    %(a)*******LINEAR INTERPOLATION FOR LONGITUDE
    if ~isequal(VMF3_data_int_h1(1,2), VMF3_data_int_h1(2,2))%if longitude must be interpolated (that is, the point does not have a longitude on the interval [0:2.5:357.5])
        zhd_lon1 = VMF3_data_int_h1(1,5) + (VMF3_data_int_h1(2,5)-VMF3_data_int_h1(1,5))*(lonDeg-VMF3_data_int_h1(1,2))/(VMF3_data_int_h1(2,2)-VMF3_data_int_h1(1,2));
        zhd_lon2 = VMF3_data_int_h1(3,5) + (VMF3_data_int_h1(4,5)-VMF3_data_int_h1(3,5))*(lonDeg-VMF3_data_int_h1(3,2))/(VMF3_data_int_h1(4,2)-VMF3_data_int_h1(3,2));
        zwd_lon1 = VMF3_data_int_h1(1,6) + (VMF3_data_int_h1(2,6)-VMF3_data_int_h1(1,6))*(lonDeg-VMF3_data_int_h1(1,2))/(VMF3_data_int_h1(2,2)-VMF3_data_int_h1(1,2));
        zwd_lon2 = VMF3_data_int_h1(3,6) + (VMF3_data_int_h1(4,6)-VMF3_data_int_h1(3,6))*(lonDeg-VMF3_data_int_h1(3,2))/(VMF3_data_int_h1(4,2)-VMF3_data_int_h1(3,2));
         ah_lon1 = VMF3_data_int_h1(1,8) + (VMF3_data_int_h1(2,8)-VMF3_data_int_h1(1,8))*(lonDeg-VMF3_data_int_h1(1,2))/(VMF3_data_int_h1(2,2)-VMF3_data_int_h1(1,2));
         ah_lon2 = VMF3_data_int_h1(3,8) + (VMF3_data_int_h1(4,8)-VMF3_data_int_h1(3,8))*(lonDeg-VMF3_data_int_h1(3,2))/(VMF3_data_int_h1(4,2)-VMF3_data_int_h1(3,2));
         aw_lon1 = VMF3_data_int_h1(1,9) + (VMF3_data_int_h1(2,9)-VMF3_data_int_h1(1,9))*(lonDeg-VMF3_data_int_h1(1,2))/(VMF3_data_int_h1(2,2)-VMF3_data_int_h1(1,2));
         aw_lon2 = VMF3_data_int_h1(3,9) + (VMF3_data_int_h1(4,9)-VMF3_data_int_h1(3,9))*(lonDeg-VMF3_data_int_h1(3,2))/(VMF3_data_int_h1(4,2)-VMF3_data_int_h1(3,2));

    else %if the station coincides with the longitude of the grid
        zhd_lon1 = VMF3_data_int_h1(1,5);
        zhd_lon2 = VMF3_data_int_h1(3,5);
        zwd_lon1 = VMF3_data_int_h1(1,6);
        zwd_lon2 = VMF3_data_int_h1(3,7);
        ah_lon1 = VMF3_data_int_h1(1,8);
        ah_lon2 = VMF3_data_int_h1(3,8);
        aw_lon1 = VMF3_data_int_h1(1,9);
        aw_lon2 = VMF3_data_int_h1(3,9);
        
    end
    
    %*****LINEAR INTERPOLATION FOR LATITUDE
    if ~isequal(VMF3_data_int_h1(1,1), VMF3_data_int_h1(3,1))
        zhd = zhd_lon1 + (zhd_lon2-zhd_lon1)*(latDeg-VMF3_data_int_h1(1,1))/(VMF3_data_int_h1(3,1)-VMF3_data_int_h1(1,1));
        zwd = zwd_lon1 + (zwd_lon2-zwd_lon1)*(latDeg-VMF3_data_int_h1(1,1))/(VMF3_data_int_h1(3,1)-VMF3_data_int_h1(1,1));
        hcoe = ah_lon1 + (ah_lon2-ah_lon1)*(latDeg-VMF3_data_int_h1(1,1))/(VMF3_data_int_h1(3,1)-VMF3_data_int_h1(1,1));
        wcoe = aw_lon1 + (aw_lon2-aw_lon1)*(latDeg-VMF3_data_int_h1(1,1))/(VMF3_data_int_h1(3,1)-VMF3_data_int_h1(1,1));
    else %if the station coincides with the latitude of the grid
        zhd = zhd_lon1;%Zenith Hydrostatic Delay
        zwd = zwd_lon1;%Zenith Wet Delay
        hcoe = ah_lon1;%Hydrostatic Coefficient
        wcoe = aw_lon1;%Wet Coefficient
    end
      
end 
%ASSIGN VALUEs FOR ALL STATIONs
ZHD = zhd;%Zenith Hydrostatic Delay
ZWD = zwd;%Zenith Wet Delay
ah  = hcoe;%Hydrostatic Coefficient
aw  = wcoe;%Wet Coefficient

%COMPUTE ZENITH TOTAL DELAY
ZTD=ZHD + ZWD;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF readVMF3grid.m  %%%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-     
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%C.SUBROUTINE TO IMPORT VMF3 OROGRAPHY FILE
function oroData = importOROGRAPHY(filename, startRow, endRow)
%**************************************************************************
%**************************************************************************
%DESCRIPTION:
%            "importOROGRAPHY" Import numeric data from a orography_ell    *
%            text file as a matrix.                                        *
%USAGE:                                                                    *
%       oroData = importOROGRAPHY(filename, startRow, endRow) Reads data   *
%       from rows STARTROW through ENDROW of text file FILENAME.           * 
%    OR                                                                    *
%       oroData = importOROGRAPHY(filename)without specifying startRow and * 
%                 endRow
%Example:
%        oroData = importOROGRAPHY(('orography_ell_5x5', 1, 2592);
%        oroData = importOROGRAPHY(('orography_ell_5x5');
%
% Auto-generated by MATLAB on 2019/06/16 18:00:03

% Initialize variables.
delimiter = ' ';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%Read columns of data as text:
%For more information, see the TEXTSCAN documentation.
formatSpec = '%s%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

% Close the text file.
fclose(fileID);

% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

% Converts text in the input cell array to numbers. Replaced non-numeric
% text with NaN.
rawData = dataArray{1};
for row=1:size(rawData, 1)
    % Create a regular expression to detect and remove non-numeric prefixes and
    % suffixes.
    regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
    try
        result = regexp(rawData(row), regexstr, 'names');
        numbers = result.numbers;
        
        % Detected commas in non-thousand locations.
        invalidThousandsSeparator = false;
        if numbers.contains(',')
            thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
            if isempty(regexp(numbers, thousandsRegExp, 'once'))
                numbers = NaN;
                invalidThousandsSeparator = true;
            end
        end
        % Convert numeric text to numbers.
        if ~invalidThousandsSeparator
            numbers = textscan(char(strrep(numbers, ',', '')), '%f');
            numericData(row, 1) = numbers{1};
            raw{row, 1} = numbers{1};
        end
    catch
        raw{row, 1} = rawData{row};
    end
end


% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

% Create output variable
oroData = cell2mat(raw);
%%%%%%%%%%%%%%%%%%%%%%%%%%END OF importOROGRAPHY.m  %%%%%%%%%%%%%%%%%%%%%%
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

%D_2.*********SUBROUTINE TO CONVERT TWO DIGITs TO YEAR 4 DIGITS
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


%1.******SUBROUTINE TO CONVERT XYZ COORDs TO GEOGRAPHIC COORDs
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

%1.1.******SUBROUTINE TO GET ELLIPSOID PARAMETERS
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
