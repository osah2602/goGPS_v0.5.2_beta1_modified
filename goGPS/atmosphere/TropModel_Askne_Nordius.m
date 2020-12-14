%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            "TropModel_Askne_Nordius" is a subroutine that Computes the   *
%             wet,dry and Total Tropospheric Delays Using  Askne and ...   *
%             Nordius Tropospheric delay correction Model                  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%USAGE:                                                                    *
%      General:                                                            *
%             [STD,SHD,SWD,ZTD,ZHD,ZWD]=TropModel_Askne_Nordius(ReceiverPos,SatPos,Pres,...
%                                             WaterVaporPres,meanTemp,WaterVaporLapseRate)
%      Others:                                                             *
%I.         [STD,SHD,SWD,ZTD,ZHD,ZWD]=_Askne_Nordius(ReceiverPos,SatPos,...*
%                                                         [P  e Tm lambda])*
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%**THE FUNCTION CALLS THE FOLLOWING SUBROUTINE(s):                         *
%1.[X,Y,Z] = geo2xyz(latitude,longitude,height,RefEllipsoid);              *
%2.[latRAD,longRAD,h,latDEC,longDEC] = xyz2LLH(X,Y,Z,RefEllipsoid);        *
%3.[a,finv] = Elipsoidpara(RefEllipsoid);                                  *
%4.[Az_rad,El_rad,Az_deg,El_deg,D]=satAzimuthElevation(UserXYZ,SatXYZ,...  *
%                                                            RefEllipsoid);*
%5.[ZTD, ZHD, ZWD] = ZenithTropDelay(T,P, es);                             *
%6.[STD, SHD, SWD] = SlantTropDelay(T,P,es,satELEV);                       *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%      Generally, the function accepts five(5) sets of inputs:             *
%A.    POSITIONs:                                                          *
%      =-=-=-=-=                                                           *
%1.    ReceiverPos : Receiver position in either Latitude,Longitude &      *
%                     height or ECEF(XYZ) Coordinates                      *
%2.    SatPos      : Satellite position(s) in ECEF(XYZ) Coordinates        *

%B.    METEOROLOGICAL PARAMETER:                                           *
%      =-=-=-=-=-=-=-=-=-=-=-=-=                                           *
%3.    Pressure(P)         : Atmospheric Pressure in millibars(mbar /hPa)  *
%4.    WaterVaporPres(e)   : Water Vapor Partial Pressure in [mbar]        *
%5.    meanTemp            : Mean atmospheric temperature(Tm) in Kelvin(K) *
%6.    WaterVaporLapseRate : water vapor lapse rate(lambda)
%OTHER CONSIDERATIONs:                                                     *
%*******Any missing Meteorological parameter, standard Meteorological value* 
%       is/are computed using the UNB3m model                              *

%*******The function also accept heights and elevation angle as inputs for *
%       Receiver position and Satellite position(s) respectively.          *
%       i.e:                                                               *
%       ReceiverPos--> height in (m)                                       *
%       SatPos  -----> elevation angle in decimal degrees                  *

%*******Where user provides empty([]) or Not a Number(nan) as entry for ...*
%       Receiver position and elevation angle as entry for Satellite...    * 
%       position(s), the function directly uses the elevation angle(s) to  *
%       compute the Mapping Function(s) for the slant delays. However,...  * 
%       where Satellite position(s) are provided instead, the subroutine   *  
%       computes only the Zenith delays(ZTD,ZHD,ZWD).                      *     
%       That is:                                                           *
%               ReceiverPos =[] or nan  and                                *
%            if                                                            *
%               SatPos  = elevation angle in decimal degrees ,             *
                                                                        
%               The ff are computed:[STD,SHD,SWD,ZTD,ZHD,ZWD]              *
%       However:                                                           *
%               ReceiverPos = [] or nan,  and                              *
%             if                                                           *
%               SatPos  = ECEF(XYZ) Coordinate                             * 

%               The ff are computed:[ZTD,ZHD,ZWD] and [STD,SHD,SWD,=0 ]    *

%*******Again,Where user provides empty([]) or Not a Number(nan) as entry  *
%       for BOTH Receiver and Satellite positions, the subroutine          *  
%       computes only the Zenith delays(ZTD,ZHD,ZWD)                       * 
%       That is:                                                           *
%               ReceiverPos =[] or nan  and ,                              *
%            if                                                            *
%               SatPos  = [] or nan ,                                      *
%               The ff are computed:[ZTD,ZHD,ZWD] and [STD,SHD,SWD,=0 ]    *

%*******Inputs can formated as:                                            *
%1.     ReceiverPos = ECEF(XYZ) Coordinate  or latitude,longitude & height *
%       SatPos  = ECEF(XYZ) Coordinate                                     *

%2.     ReceiverPos = ECEF(XYZ) Coordinate  or latitude,longitude & height *
%       SatPos  = elevation angle in decimal degrees                       *

%3.     ReceiverPos = height in (m)                                        *
%       SatPos  = elevation angle in decimal degrees                       *

%4.     ReceiverPos = empty ( [ ] ) or nan                                 *
%       SatPos  = elevation angle in decimal degrees                       *

%       Any entry or input other than what is specified above (i.e. 1-4),  *
%       the program terminates by giving an error message alert !          *

%*******For three(3)inputs,the function checks the number of columns in    *
%       Temperature and if the number of columns are three(3),then all     *
%       meteolorogical parameters are available or provided.In that case,  *  
%       Temperature(Temp) will be in column 1, Pressure in column 2 &      *
%       RelHumidity in column 3. i.e.[Temperature  Pressure  RelHumidity]. *

%*******For two(2)inputs,the function,examines the columns in the first two*
%       input arguements( ReceiverPos and SatPos),that is if nargin==2,and * 
%       if the columns are single(or nx1) decimal values each in ReceiverPos    
%       and SatPos,then the function assumes height  and  elevation angle  * 
%       respectively as inputs.                                            *
%       i.e. where:    ht = ReceiverPos i.e. Ellipsoidal Height            *
%                   satEL = SatPos i.e. satellite elevation (satEL)        *
%       otherwise:, entry will be seen as Receiver & satellite positions   *

%*******For five (5) inputs with single(or nx1) decimal values each in     *
%       ReceiverPos and SatPos, the function  still assumes height  and    *
%       elevation angle respectively as inputs for ReceiverPos and SatPos  *

%*******For only one(1) input, the program examines the  ReceiverPos input *
%       and extract the user positions.It also employs standard atmospheric* 
%       parameters [1013.25(mbar) 291.15(K) 50(%)] in the computation.
%       Here again, only the Zenith delays are computed i.e.[ZTD,ZHD,ZWD]

%*******For no input, error message is given and the program               *
%       Terminates / stop working given empty matrices([]) or zeros as     * 
%       outputs                                                            *                                                        

%REMARKS:                                                                  *
%        All these considerations are in away, other input formats one can *
%        consider when using this subroutine.                              *            
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
%REFERENCE:                                                                *
%          Askne and Nordius, Estimation of tropospheric delay for         *
%          microwaves from surface weather data, Radio Science, Vol 22(3): * 
%          379-386, 1987.                                                  *
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%**************************************************************************+
%**************************************************************************+

function [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_Askne_Nordius(ReceiverPos,SatPos,Pres,...
                                                             WaterVaporPres,meanTemp,WaterVaporLapseRate)
%**********CHECK INPUTs & REFORMAT INPUT DATA
switch nargin
    
    case {6,5,4,3,2,1} %When all inputs are provided        
        
        if (any(nargin==[6,5,4,3,2,1]))
            
            if nargin == 6
               %******Assignments   
               P = Pres;%Atmospheric Pressure in millibars(mbar)/Hecto pascal(hpa)
              es = WaterVaporPres;%partial water vapor pressure(e)  
              Tm = meanTemp; % mean temperature in Kelvin
          lambda = WaterVaporLapseRate; %water vapor lapse rate (see definition in Askne and Nordius 1987)
          
            elseif nargin == 5
                   %******Assignments   
                   P = Pres;%Atmospheric Pressure in millibars(mbar)/Hecto pascal(hpa)
                  es = WaterVaporPres;%partial water vapor pressure(e)  
                  Tm = meanTemp; % mean temperature in Kelvin
              lambda = []; %Assign empty matrix([]) to water vapor lapse rate 
               
            elseif nargin == 4 %for four(4) or all inputs...              
                   %******Assignments   
                   P = Pres;%Atmospheric Pressure in millibars(mbar)/Hecto pascal(hpa)
                  es = WaterVaporPres;%partial water vapor pressure(e) 
                  Tm = []; %Assign empty matrix([]) to mean temperature 
              lambda = []; %Assign empty matrix([]) to water vapor lapse rate 
                          
          elseif nargin==3 %For three(3) inputs...                              
                %CHECK THE NUMBER OF COLUMNs IN Temperature
                Tcol=size(Pres,2); 
                
                switch Tcol 
                    
                    case 1 %if num of col is 1
                         P = Pres;% Atmospheric Pressure 
                        es = [];%Assign empty matrix([]) to Relative Humidity 
                        Tm = []; %Assign empty matrix([]) to mean temperature 
                    lambda = []; %Assign empty matrix([]) to water vapor lapse rate 
                        
                    case 2 %if num of col is 2
                         P = Pres(:,1);%Assign 1st Column of Pres to Atmospheric Pressure              
                        es = Pres(:,2);%Assign 2nd Column of Pres to partial water vapor pressure(e)        
                        Tm = []; %Assign empty matrix([]) to mean temperature 
                    lambda = []; %Assign empty matrix([]) to water vapor lapse rate 
                        
                    case 3 %if num of col is 3
                         P = Pres(:,1);%Assign 1st Column of Pres to Atmospheric Pressure              
                        es = Pres(:,2);%Assign 2nd Column of Pres to partial water vapor pressure(e)                   
                        Tm = Pres(:,3);%Assign 3rd Column of Pres to mean temperature
                    lambda = []; %Assign empty matrix([]) to water vapor lapse rate   
                    
                    case 4 %if num of col is 4
                         P = Pres(:,1);%Assign 1st Column of Pres to Atmospheric Pressure              
                        es = Pres(:,2);%Assign 2nd Column of Pres to partial water vapor pressure(e)                   
                        Tm = Pres(:,3);%Assign 3rd Column of Pres to mean temperature
                    lambda = Pres(:,4);%Assign 4th Column of Pres to water vapor lapse rate  
 
                    case 5 %if num of col is 5
                         T = Pres(:,1);%Assign 1st Column of Pres in Degree Celcius                  
                         P = Pres(:,2);%Assign 2nd Column of Pres to Atmospheric Pressure                   
                        es = Pres(:,3);%Assign 3rd Column of Pres to partial water vapor pressure(e) 
                        Tm = Pres(:,4);%Assign 4th Column of Pres to mean temperature 
                    lambda = Pres(:,5);%Assign 5th Column of Pres to water vapor lapse rate   
                   
                    case 0
                         T = [];%Assign empty matrix([]) to Temperature in Degree Celcius
                         P = [];%Assign empty matrix([]) to Atmospheric Pressure 
                        es = [];%Assign empty matrix([]) to partial water vapor pressure(e)
                        Tm = []; %Assign empty matrix([]) to mean temperature 
                    lambda = []; %Assign empty matrix([]) to water vapor lapse rate
                          
  
                end %switch Tcol  
                
          elseif  nargin==2 %For two(2) inputs...              
                %*****Assignments
                T = [];%Assign empty matrix([]) to Temperature in Degree Celcius
                P = [];%Assign empty matrix([]) to Atmospheric Pressure 
               es = [];%Assign empty matrix([]) to partial water vapor pressure(e)
               Tm = []; %Assign empty matrix([]) to mean temperature 
           lambda = []; %Assign empty matrix([]) to water vapor lapse rate 
                
          elseif  nargin==1 %for single input            
                %*****Assignments                
                T = [];%Assign empty matrix([]) to Temperature in Degree Celcius
                P = [];%Assign empty matrix([]) to Atmospheric Pressure 
               es = [];%Assign empty matrix([]) to partial water vapor pressure(e)
               Tm = []; %Assign empty matrix([]) to mean temperature 
           lambda = []; %Assign empty matrix([]) to water vapor lapse rate           
           SatPos = [];%Assign empty matrix([]) to SatPos
                                               
            end   %if nargin == 7
          
        end  %if (any(nargin==[7,6,5,4,3,2,1]))              
        %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-             
       
         %*************IF INPUTS ARE NOT EMPTY([])     
        %1.RECEIVER POSITION
        if ~isempty(ReceiverPos) | ~isnan(ReceiverPos)
            
           %***ROUND TO THE NEAREST THOUSAND 
           t1=roundn(ReceiverPos,4) ;
           
           if size(ReceiverPos,2)==1                
             if all(t1(:,1))==0                    
               h=ReceiverPos;%Assigning ReceiverPos to h(Elipsoidal Height)
             end 
             
           else
               Rpos=ReceiverPos;%Receiver position(XYZ/LAT LONG h) 
                        
           end %//if size(ReceiverPos,2)==1
           
        end %\\if any([flag_empty1==0,flag_nan1==0])
          
        %2.SATELLITEs POSITION
        if ~isempty(SatPos) | ~isnan(SatPos)
            
           %***ROUND TO THE NEAREST THOUSAND            
           t2=roundn(SatPos,3);
                  
             if  size(SatPos,2)==1 %if the number of columns in SatPos is 1                                                      
               if  all(t2(:,1))==0                       
                 satEL_deg=SatPos;%Assigning SatPos to satEL(Satellite Elevation)
                 satEL_rad=satEL_deg.* pi / 180;%Elevation Angles converted to radian
                 
               end    %if all(t2(:,1))==0 
               
             else    
                Satpos=SatPos;%Satellite position in XYZ                   
             end  %\\if size(SatPos,2)==1 
             
        end %//if any([flag_empty2==0,flag_nan2==0])
        
        %******************************************************************
        %         EXTRACT USER INPUT RECEIVER Coordinates,.....
        %******************************************************************
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
              latD=lat;%Latitude in Decimal degrees
              lonD=lon;%Longitude in Decimal degrees 
                 
          end
          
          latR=latD.*degrad; %latitude in radian
          lonR=lonD.*degrad; %longitude in radian
          
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
           if exist('Satpos','var')                         
         
             %Call the "satAzEl.m" function 
             [satEL_deg,satEL_rad]=satAzEl([X Y Z],Satpos);%compute satellite elevation angle
             
           end
         
        end %if ~exist('satEL','var'
        
        else
            if exist('h','var')  
               latD = zeros(h); %Assign zeros of h dimension to latD
              [lonD,latR,lonR] = deal(latD); %Create copy of latD 
            else
                h = 0 ; %Assign 0 to h 
                [latD,lonD,latR,lonR] = deal(h); %Create copy of latD
            end
            %SAVE COORDs for use
            setappdata(0,'latD',latD)
            setappdata(0,'lonD',lonD)
            setappdata(0,'latR',latR)
            setappdata(0,'lonR',lonR)
            setappdata(0,'ht',h)
        
        end %if exist('Rpos','var')
                                                                       
    otherwise
             %ISSUE ERROR MESSAGE INPUT IS ONE
              beep%Give a beep sound 
              errmsg{1}=sprintf('Insuficient Data Input / Wrong Data Input format .');
              errmsg{2}='';
              errmsg{3}='Please Check file / Data format & Try Again.';
              errordlg(errmsg,'Coordinate(s) Input Error','modal')  
              
              %Return empty ([]) outputs
              STD=[]; SHD=[]; SWD=[]; ZTD=[]; ZHD=[]; ZWD=[];
              
              return                          
end %switch nargin

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%************TROPOHERIC DELAY MODELING USING SAASTAMOINEN MODEL
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%USE UNB3m MODEL TO COMPUTE ANY MISSING METEOTOLOGICAL PARAMETER(S)
%********GET METEOROLOGICAL PARAMETERS
if all([isempty(P) isempty(es) isempty(Tm) isempty(lambda)])
    
   %RETRIEVE STORED UTC TIME from "RINEX_get_epoch.m"
    UTCtime = getappdata(0,'date');  
    
   if all([exist('Rpos','var'),~isempty(UTCtime)])
       
      %Generate Meteorological parameters from UNB3m model
      [~,P,es,Tm,~,lambda] = getMETpara_UNB3m(UTCtime,latD,lonD,h);
         
       %NOTE:ALL UNITS ARE IN STANDARD FORMS.
       
   else  
        %ISSUE ERROR MESSAGE  
        beep%Give a beep sound 
        errmsg1{1}=sprintf('Insuficient Inputs for Meteorological Parameters.\n');
        errmsg1{2}=sprintf('The following parameters have not been provided :.\n');
        errmsg1{3}=sprintf('1. Surface Pressure ( P ) .\n');
        errmsg1{4}=sprintf('2. Water Vapour Partial Pressure ( e ) .\n');
        errmsg1{5}=sprintf('3. Mean Temperature of Water Vapor ( Tm ).\n');
        errmsg1{6}=sprintf('4. Water Vapour Lapse Rate ( lambda ).\n');
        errmsg1{7}=sprintf('Please provide these Meteorological parameters & Try again.');
        errordlg(errmsg1,'Meteorological parameters Input Error','modal') 
                    
        %RETURN EMPTY([]) DELAYS
        STD =[]; SHD = []; SWD = []; ZTD = []; ZHD = []; ZWD = [];
                  
        return
           
   end 
            
elseif  all([~isempty(P) isempty(es) isempty(Tm) isempty(lambda)])
        
        %RETRIEVE STORED UTC TIME from "RINEX_get_epoch.m"
        UTCtime = getappdata(0,'date');  
    
        if all([exist('Rpos','var'),~isempty(UTCtime)])
            
           %Generate Meteorological parameters from UNB3m model
           [~,~,es,Tm,~,lambda] = getMETpara_UNB3m(UTCtime,latD,lonD,h);
              
        else  
             %ISSUE ERROR MESSAGE  
             beep%Give a beep sound 
             errmsg2{1}=sprintf('Insuficient Inputs for Meteorological Parameters.\n');
             errmsg2{2}=sprintf('The following parameters have not been provided :.\n');
             errmsg2{3}=sprintf('1. Water Vapour Partial Pressure ( e ) .\n');
             errmsg2{4}=sprintf('2. Mean Temperature of Water Vapor ( Tm ).\n');
             errmsg2{5}=sprintf('3. Water Vapour Lapse Rate ( lambda ).\n');
             errmsg2{6}=sprintf('Please provide these Meteorological parameters & Try again.');
             errordlg(errmsg2,'Meteorological parameters Input Error','modal')
               
             %RETURN EMPTY([]) DELAYS
             STD =[]; SHD = []; SWD = []; ZTD = []; ZHD = []; ZWD = [];
                  
             return
           
        end  
           
elseif  all([~isempty(P) ~isempty(es) isempty(Tm) isempty(lambda)])
        
        %RETRIEVE STORED UTC TIME from "RINEX_get_epoch.m"
        UTCtime = getappdata(0,'date');  
    
        if all([exist('Rpos','var'),~isempty(UTCtime)])
           %Generate Meteorological parameters from UNB3m model
           [~,~,~,Tm,~,lambda] = getMETpara_UNB3m(UTCtime,latD,lonD,h);
              
        else  
             %ISSUE ERROR MESSAGE  
             beep%Give a beep sound 
             errmsg3{1}=sprintf('Insuficient Inputs for Meteorological Parameters.\n');
             errmsg3{2}=sprintf('The following parameters have not been provided :.\n');
             errmsg3{3}=sprintf('1. Mean Temperature of Water Vapor ( Tm ).\n');
             errmsg3{4}=sprintf('2. Water Vapour Lapse Rate ( lambda ).\n');
             errmsg3{5}=sprintf('Please provide these Meteorological parameters & Try again.');
             errordlg(errmsg3,'Meteorological parameters Input Error','modal')
               
             %RETURN EMPTY([]) DELAYS
             STD =[]; SHD = []; SWD = []; ZTD = []; ZHD = []; ZWD = [];
                  
             return
           
        end   
           
elseif  all([~isempty(P) ~isempty(es) ~isempty(Tm) isempty(lambda)])
        
        %RETRIEVE STORED UTC TIME from "RINEX_get_epoch.m"
         UTCtime = getappdata(0,'date');  
    
        if all([exist('Rpos','var'),~isempty(UTCtime)])
            
           %Generate Meteorological parameters from UNB3m model
           [~,~,~,~,~,lambda] = getMETpara_UNB3m(UTCtime,latD,lonD,h);
           
        else  
              %ISSUE ERROR MESSAGE  
               beep%Give a beep sound 
               errmsg4{1}=sprintf('Insuficient Inputs for Meteorological Parameters.\n');
               errmsg4{2}=sprintf('The following parameter has not been provided :.\n');
               errmsg4{3}=sprintf('1. Water Vapour Lapse Rate ( lambda ).\n');
               errmsg4{4}=sprintf('Please provide this Meteorological parameter & Try again.');
               errordlg(errmsg4,'Meteorological parameter Input Error','modal')
               
               %RETURN EMPTY([]) DELAYS
               STD =[]; SHD = []; SWD = []; ZTD = []; ZHD = []; ZWD = [];
                  
               return
           
        end   
              
                   
end   %if if all([isempty(T) isempty(P) isempty(es) isempty(Tm) isempty(lambda)])
       
%**********COMPUTE TOTAL TROPOSPHERIC DELAY IN THE ZENITH  DIRECTION
%Call the "ZenithTropDelay.m" function
[ZTD,ZHD,ZWD] = ZenithTropDelay(P,es,Tm,lambda);

%**********COMPUTE TOTAL TROPOSPHERIC DELAY IN THE SATELLITE DIRECTION(LOS)
 if exist('satEL_deg','var')
     
%Call the "SlantTropDelay.m" function
[STD,SHD,SWD] = SlantTropDelay(ZHD,ZWD,satEL_deg);

 else
     %Create zero matrix     
     STD = zeros(size(ZTD));
     SHD = zeros(size(ZTD));
     SWD = zeros(size(ZTD));    

 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF TropModel_Askne_Nordius.m 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%A.******SUBROUTINE TO COMPUTE ZENITH TROPOSPHERIC DELAY
function [ZTD, ZHD, ZWD] = ZenithTropDelay(P,es,Tm,lambda)
%**************************************************************************
%DESCRIPTION:
%            "ZenithTropDelay" Computes Zenith Tropospheric Delay by ...   * 
%             Askne and Nordius model.                                     *
%******SYNTAX:
%             [ZTD, ZHD, ZWD] = ZenithTropDelay(T,P,es,Tm,lambda);         *
%
%******INPUT:   
%            P = atmospheric pressure [mbar/hPa]                           *                       *
%           es = partial pressure of water vapor [mbar/hPa]                *
%           Tm = mean temperature in Kelvin                                *
%       lambda = water vapor lapse rate                                    *
%
%****OUTPUT:
%           ZHD = Zenith Hydrostatic Delay                                 *
%           ZWD = Zenith Wet Delay                                         *
%           ZTD = Zenith Total Delay                                       *

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%*****RETRIEVE STORED DATA
latD = getappdata(0,'latD'); %Get station latitude Coords [degrees]
lonD = getappdata(0,'lonD'); %Get station latitude Coords [degrees]
latR = getappdata(0,'latR'); %Get station latitude Coords [radians]
   h = getappdata(0,'ht'); %Get station height [meters]

%********COMPUTE ZENITH TROPOSPHERIC DELAYs

%*************COMPUTE ORTHOMETRIC HEIGHT(H)
%--------------------------------------------------------------------------  
%NOTE:
%     SAASTAMOINEN MODEL REQUIRE +VE ORTHOMETRIC HEIGHT

%COMPUTE CORRECTION FACTOR FOR LOCAL GRAVITY ACCELERATION 

if (~isempty(latR) && ~isempty(h))
    
 %COMPUTE UNDULATION & ORTHOMETRIC HEIGHT 
  [~,horth] = GETundulation(latD,lonD,h);  
  
   dgref = 1-0.00266.*cos(2.*latR)-0.00028e-3.*horth;

else
    dgref = 1;
    
end %//if (~isempty(latR) && ~isempty(h))

%*****DEFINE CONSTANTs
k1   = 77.604; 				   % K/hPa 
k2   = 64.79; 				   % K/hPa
Md   = 28.9644; %Molar Mass of Dry Air  in [g/mol]   
Mw   = 18.0152;%Molar Mass of Water in [g/mol]
k2p  = k2 - k1*Mw/Md;          % K/hPa  [called k2 prime]
k3   = 377600; 				   % KK/hPa
dMtr = Md*10^-3;%molar mass of dry air in kg/mol     
R    = 8.3143;%universal gas constant in J/K/mol          
Rd   = R/dMtr ; %specific gas constant for dry consituents 

%ACCELERATION GRAVITY AT THE ATMOSPHERIC COLUMN
gm = 9.784 .* dgref;%[m/s^2]

%*****COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
if size(P,2) > 1
    
   for i = 1:size(P,2)
             
       %*****COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
       ZHD(:,i) = (1.0e-6.*k1.*Rd.*P(:,i))./gm;%[m] 
   end  
         
else   
     %*****COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
     ZHD = (1.0e-6.*k1.*Rd.*P)./gm;%[m]  
          
end  %//if size(P,2) > 1

%*****COMPUTE ZENITH WET DELAY(ZWD)   
if any([size(Tm,2) > 1,size(es,2) > 1])
    
   for i = 1 : size(Tm,2)
      
       try
          ZWD(:,i) = ((1e-6.*(k2p + k3./Tm(:,i)).*Rd)./(gm.*(lambda(:,i) + 1))).*es(:,i);%[m]
          
       catch %ALTERNATIVE
           ZWD(:,i) = 1e-6.*(k2p + k3./Tm(:,i)).*Rd./(lambda(:,i) + 1)./gm.*es(:,i);%[m]
           
       end
   end  
  
else 
     %*****COMPUTE ZENITH WET DELAY(ZWD)
     try
        ZWD = ((1e-6.*(k2p + k3./Tm).*Rd)./(gm.*(lambda + 1))).*es;%[m]
        
     catch  %ALTERNATIVE
           ZWD = 1e-6.*(k2p + k3./Tm).*Rd./(lambda + 1)./gm.*es;%[m]
     end
end  

%****FINALLY, COMPUTE ZENITH TOTAL DELAY(ZTD)
if any([size(ZHD,2) > 1,size(ZWD,2) > 1]) 
    
   for k = 1: size(ZHD,2)
       
       ZTD(:,k) = ZHD(:,k) + ZWD(:,k);%#ok<*AGROW> %Zenith Total Delay [m]
   end  
   
else  
     ZTD = ZHD + ZWD;%Zenith Total Delay [m]
end  

%***********************END OF ZenithTropDelay.m ***********************    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%B.SUBROUTINE TO COMPUTE SLANT TROPOSPHERIC DELAY (SATELLITE IN VIEW)
function [STD, SHD, SWD] = SlantTropDelay(ZHD,ZWD,satELEV)

%DESCRIPTION:
%            "SlantTropDelay" Computes Slant Tropospheric Delay by ...     * 
%             Askne and Nordius model.                                     *
%******SYNTAX:                                                             *
%             [STD, SHD, SWD] = SlantTropDelay(ZHD,ZWD,satELEV);           *
%
%******INPUT:
%            ZHD = Zenith Hydrostatic Delay [m]                            *
%            ZWD = Zenith Wet Delay [m]                                    *
%        satELEV = Satellite Elevation Angle in degrees                    *

%      satELEV is in the format (nxm);
%Where:
%      n = Number of rows representing each satellite
%      m = Number of columns representing each station/receiver
%    
%****OUTPUT:
%           SHD = Slant Hydrostatic Delay                                  *
%           SWD = Slant Wet Delay                                          *
%           STD = Slant Total Delay                                        *

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%*********COMPUTE SLANT TROPOSPHERIC DELAY (SATELLITE IN VIEW)

%*****INITIALIZING OUTPUT VARIABLEs
SHD = zeros(size(satELEV,1),size(satELEV,2));%Assign zeros of nxm to SHD
[STD,SWD] = deal(SHD);%copy the contents of SHD to all the requested outputs

%COMPUTE MAPPING FUNCTION USING BLACK Using Black & Eisner(1984) MODEL
[MFh,MFw] = Black_Eisner_MF(satELEV);

%FIND NUMBER OF ROWS & COLUMNs IN satELEV
[nrow,ncol] = size(satELEV);

for i = 1:ncol %Loop over the Number of Receiver positions
    
   for j = 1:nrow %Loop over the Number of Satellite Elevations
    
       %COMPUTE Slant Hydrostatic Delay(SHD)
       SHD(j,i) = ZHD(i,1).*MFh(j,i);

       %COMPUTE Slant Wet Delay(ZWD)
       SWD(j,i) = ZWD(i,1).*MFw(j,i);

       %****COMPUTE TOTAL SLANT TROPOSPHERIC DELAY(STD)
       STD(j,i) = SHD(j,i)+ SWD(j,i);            
                                                    
   end  %for j=1:nrow
          
end     %for i=1:ncol

%***********************END OF SlantTropDelay.m ***********************    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
%%%%%%%%%%%%%%%%%%%%%%%%END OF TropModel_SAAS.m %%%%%%%%%%%%%%%%%%%%%%% 

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

%5.******SUBROUTINE TO COMPUTE METEOROLOGICAL PARAMETERS USING UNB3m MODEL
function [T, P, E,TM,beta,lambda] = getMETpara_UNB3m(UTCtime,lat,lon,h)
%**************************************************************************
%DESCRIPTION:                                                              * 
%           "getMETpara_UNB3m" computes surface Meteorological parameters  *
%           temperature,pressure,relative humidity and mean temperature... *
%           of water vapor for a given latitude, height and day of year.   *
%USAGE:                                                                    *
%      [T, P, E,TM,beta,lambda]=getMETpara_UNB3m(UTCtime,lat,lon,h)        *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%1.    UTCtime : Receiver reception time in[year,Month,Day,Hour,Min,Sec] OR*
%                                          [year,Month,Day]                *
%2.    lat     : Receiver ellipsoidal latitude in [degrees]                *
%      lon     : Receiver ellipsoidal longitude in [degrees]               *
%      h       : Receiver elellipsoidal height in [meters]                 *
%--------------------------------------------------------------------------*
%***OUTPUTs:                                                               *
%1.         T : Surface temperature in kelvin(k)                           *
%2.         P : Surface pressure in millibar(mbar)                         *
%3.         E : Surface water vapor pressure in millibar(mbar)             *
%4.        TM : Mean temperature of water vapor in kelvin(k)               *
%5.        beta : temperature `lapse' rate in kelvin per meter[(K/m)]      *
%6      lambda : water vapour `lapse rate'(dimensionless)]                 *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%GET SIZE OF USER INPUT TIME & POSITION
nrow_time=size(UTCtime,1);
nrow_pos=size(latD,1);

%1.******UTCdate
Yr  = UTCtime(:,1);%get Year
Mn  = UTCtime(:,2);%get Month
Day = UTCtime(:,3);%get Day

%*****INITIALIZE OUTPUTs 
P=zeros(nrow_pos,nrow_time);
    
[T,E,TM,beta,lambda]=deal(P);%create copy of P in T,E,TM

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
      
   for i=1:nrow_time
      
      %Compute Met Parameters 
      [T(i,1), P(i,1), E(i,1),TM(i,1),beta(i,1),lambda(i,1)] = getMETpara(UTCtime(i,:),lat(i),lon(i),h(i)); 
   end
   
  else
       for i = 1:nrow_time
        
        for j = 1:nrow_pos
            
            %Compute Met Parameters 
            [T(j,i), P(j,i), E(j,i), TM(j,i),beta(j,i),lambda(j,i)] = getMETpara(UTCtime(i,:),lat(j),lon(j),h(j));    
        end 
       end 
    
  end  %//if isequal(Identical,1) 

else
     
    for i = 1:nrow_time
        
        for j = 1:nrow_pos
            
            %Compute Met Parameters 
            [T(j,i), P(j,i), E(j,i), TM(j,i),beta(j,i),lambda(j,i)] = getMETpara(UTCtime(i,:),lat(j),lon(j),h(j));    
        end
    end
    
end %//if isequal(nrow_time,nrow_pos)
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% getMETpara_UNB3m.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%5A.SUBROUTINE TO COMPUTE METEOROLOGICAL PARAMETERs 
function [T, P, E,TM,beta,lambda] = getMETpara(UTCtime,lat,lon,hgt)
%**************************************************************************
%DESCRIPTION:                                                              * 
%           "getMETpara" computes surface Meteorological parameters        *
%           temperature,pressure,relative humidity and mean temperature... *
%           of water vapor for a given latitude, height and day of year.   *
%USAGE:                                                                    *
%     [T, P, E,TM,beta,lambda]=getMETpara(UTCtime,lat,lon,hgt)             *      
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%****INPUTs:                                                               *
%1.       UTCtime:.........UTC time in [Year,Month,Day,Hour,Minute,Seconds]*                         
%2.           lat:.........Station geodetic latitude in [degrees]          *
%3.           lon:.........station ellipsoidal longitude in [degrees]      *
%4.           hgt:.........Station ellipsoidal height in [meters]          *
%NOTE: intAVG & intAMP are in the format:
%             1       2     3       4       5                              *
%     intAVG:[P       T     RH     beta   lambda]                          *
%     intAMP:[P       T     RH     beta   lambda]                          *

%***OUTPUTs:                                                               *
%1.         T : Surface temperature in kelvin(k)                           *
%2.         P : Surface pressure in millibar(mbar)                         *
%3.         E : Surface water vapor pressure in millibar(mbar)             *
%4.        TM : Mean temperature of water vapor in kelvin(k)               *
%        beta : temperature `lapse' rate in kelvin per meter[(K/m)]        *
%      lambda : water vapour `lapse rate'(dimensionless)]                  * 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-    
%==========================================================================+
%REFERENCE:                                                                *
%1.        Leandro R.F., Santos, M.C. and Langley R.B., (2006).            *
%          UNB Neutral Atmosphere Models: Development and Performance.     * 
%          Proceedings of ION NTM 2006, pp. 564-573, Monterey, California, * 
%          January, 2006.                                                  *
%2.        Collins, P., Langley, R., & LaMance, J. (1996, June). Limiting 
%          factors in tropospheric propagation delay error modelling for 
%          GPS airborne navigation. In Proceedings of The Institute of 
%          Navigation 52nd Annual Meeting (Vol. 10). Cambridge: MA, USA.
%3.        IERS (2004). IERS Conventions (2003), eds. D.D. McCarthy and    * 
%          G. Petit, IERS Technical Note No. 32, International Earth       *
%          Rotation and Reference Systems Service, Verlag des Bundesmates  * 
%          für Kartographie und Geodasie, Frankfurt am Main.               *
% =========================================================================+
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================+
%**************************************************************************+
%*******************Initialize UNB3m LOOK-UP TABLE
%(1)AVG:AVERAGE METEOROLOGICAL PARAMETERS
AVG=[  15.0  1013.25  299.65  75.00  6.30e-3  2.77
       30.0  1017.25  294.15  80.00  6.05e-3  3.15
       45.0  1015.75  283.15  76.00  5.58e-3  2.57
       60.0  1011.75  272.15  77.50  5.39e-3  1.81
       75.0  1013.00  263.65  82.50  4.53e-3  1.55];
      %lat     Po      To      RH     beta0  lambda0
%%WHERE:
%      lat: Table range of latitude values in degrees(deg)
%       Po: average Pressure in mbar
%       To: average temperature in kelvin(k)
%       RH: average Relative Humidity in percent(%)
%    beta0: average temperature `lapse' rate in kelvin per meter[(K/m)]
%   lamda0: average water vapour `lapse rate'(dimensionless)]

AMP=[  15.0   0.00   0.00   0.00  0.00e-3  0.00
       30.0  -3.75   7.00   0.00  0.25e-3  0.33
       45.0  -2.25  11.00  -1.00  0.32e-3  0.46
       60.0  -1.75  15.00  -2.50  0.81e-3  0.74
       75.0  -0.50  14.50   2.50  0.62e-3  0.30];
     % lat     dP    dT     dRH    dbeta  dlambda
%WHERE:
%      lat: Table range of latitude values in degrees(deg)
%      dP : Seasonal variation in Pressure(mbar)
%      dT : Seasonal variation in temperature(kelvin(k))
%     dRH : Seasonal variation in Relative Humidity in percent(%)
%    dbeta: Seasonal variation in temperature `lapse' rate(K/m)]
%   dlamda: Seasonal variation in water vapour `lapse rate'(dimensionless)]     
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%GET UTC TIME COMPONENTs

if size(UTCtime,2) == 6
  %1.******UTCdate
  Yr = UTCtime(:,1);%get Hour
  Mn = UTCtime(:,2);%get Month
 Day = UTCtime(:,3);%get Day
%2.UTCtime
   H = UTCtime(:,4);%get Hour
 MIN = UTCtime(:,5);%get Minute
SECs = UTCtime(:,6);%get Seconds

elseif size(UTCtime,2) == 3
       %1.******UTCdate
       Yr = UTCtime(:,1);%get Hour
       Mn = UTCtime(:,2);%get Month
      Day = UTCtime(:,3);%get Day      
     
      %2.UTCtime
      H = zeros(size(Yr,1),1);%Assign zeros to Hour
    MIN = zeros(size(Yr,1),1);%Assign zeros to Minute
   SECs = zeros(size(Yr,1),1);%Assign zeros to Seconds
   
end

%***UNB3m REQUIRES +VE ORTHOMETRIC HEIGHT
[~,horth] = GETundulation(lat,lon,hgt);

%**************DEFINE CONSTANTs  
g   = 9.80665;%Surface acceleration of gravity in m/s^2
R   = 287.054;% The gas constant for dry air in (J/kg/K); 

%****DEFINE pi
pi = 3.14159265358979356d0;

%Coefficient for changing DOY to radian 
doy2rad = 2*pi/365.25d0; 
latRAD=lat.*pi/180;%Latitude in deg converted to radians
        
%***COMPUTE DAY OF YEAR(DOY)                  
%Call the "utc2JulianDay_DoY.m" function
[~, ~,DOY]=utc2JulianDay_DoY(Yr,Mn,Day,H,MIN,SECs);           
   
%************DEAL WITH SOUTHERN HEMISPHERE & YEARLY VARIATION
%NOTE:  
      %Since no data from the southern hemisphere were used in developing
      %these mapping functions, the inversion of the seasons has been
      %accounted for simply by adding half a year to the phase for
      %southern latitudes. 
if (lat < 0) %if Southern Hemisphere
   DOY = DOY + 365.25/2; %Add half a year to the phase
end  
%***COMPUTE COS OF PHASE ANGLE
COSphase = cos((DOY - 28) * doy2rad );
 
%*********INTERPOLATE METEOROLOGICAL PARAMETERs
%1.Compute Average & seasonal Values
[intAVG,intAMP] =interpolate_UNB3m(lat,AVG,AMP);

%COMPUTE SURFACE TROPOSPHERIC VALUEs
   P0 = intAVG(1,1) - intAMP(1,1) * COSphase;
    T0 = intAVG(1,2) - intAMP(1,2) * COSphase;
   RH0 = intAVG(1,3) - intAMP(1,3) * COSphase;
  BETA = intAVG(1,4) - intAMP(1,4) * COSphase;
LAMBDA = intAVG(1,5) - intAMP(1,5) * COSphase;
    
%*****TRANSFORM FROM RELATIVE HUMIDITY(RH) TO WVP (IERS Conventions 2003) 
ES = 0.01 * exp(1.2378847e-5 * (T0 ^ 2) - 1.9121316e-2 * ...
    T0 + 33.93711047 - 6.3431645e3 * (T0 ^ -1));%Saturation Vapour Pressure
FW = 1.00062 + 3.14e-6 * P0 + 5.6e-7 * ((T0 - 273.15) ^ 2);%Enhancement factor
E0 = (RH0 / 100) * ES * FW;%Water Vapor Pressure in [mbar or hpa]

%***COMPUTE POWER VALUE FOR PRESSURE & WATER VAPOUR
EP = g / (R * BETA);

%*****SCALE SURFACE VALUEs TO USER'S HEIGHT
T = T0 - BETA * horth; %Surface Temperature at user Location
P = P0 * ( T / T0 ) ^ EP; %Surface Pressure at user Location
E = E0 * ( T / T0 ) ^ ( EP * (LAMBDA+1) );%water vapor pressure at user Location
beta=BETA;
lambda=LAMBDA;
%*****COMPUTE ACCELERATION GRAVITY AT THE ATMOSPHERIC COLUMN
dgref = 1 - 2.66e-3*cos(2*latRAD) - 2.8e-07 * horth;
gm   = 9.784 * dgref;

%*********COMPUTE MEAN TEMPERATURE OF WATER VAPOR(TM)
TM  = T * (1 - (BETA * R / (gm*(LAMBDA+1))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF getMETpara.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%5B.SUBROUTINE TO INTERPOLATE UNB3m AVEARAGE & AMPLITUDE VALUEs
function [AVG,AMP] =interpolate_UNB3m(UserLat,avgValues,ampValues)
                                                                 
%***************************************************************************
%DESCRIPTION:                                                              * 
%            The function  "interpolate_UNB3m" Interpolates and computes   * 
%            averaged and seasonal variation Meteorological parameters...  *   
%            used for tropospheric delay prediction, given receiver...     * 
%            latitude based on the UNB3m Model.Parameters above |lat|<=15° *   
%            and |lat|>=75° are extracted directly  while Parameters for...* 
%            latitudes 15°<|lat|<75° are linearly interpolated between...  *
%            values of two closest latitudes.                              *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%1.       UserLat: Station/User geodetic latitude(s) in degrees(deg)(n x 1)*
%2.     avgValues: Average Meteorological Values                           *
%3.     ampValues: Amplitude/Seasonal Meteorological Values                *

%%*********************LOOK-UP TABLE  of UNB3m model                       *
%1.*****AVERAGE METEOROLOGICAL PARAMETERS                                  *                                              
%avg=[  15.0  1013.25  299.65  75.00  6.30e-3  2.77                        *
%       30.0  1017.25  294.15  80.00  6.05e-3  3.15                        *
%       45.0  1015.75  283.15  76.00  5.58e-3  2.57                        *
%       60.0  1011.75  272.15  77.50  5.39e-3  1.81                        *
%       75.0  1013.00  263.65  82.50  4.53e-3  1.55  ];                    *
%       lat     Po       To     RH     beta0   lambda0                     *

%2.*****AMPLITUDE/SEASONAL METEOROLOGICAL PARAMETERS                       *
%amp=[ 15.0   0.00   0.00   0.00  0.00e-3  0.00                            *
%      30.0  -3.75   7.00   0.00  0.25e-3  0.33                            *
%      45.0  -2.25  11.00  -1.00  0.32e-3  0.46                            *
%      60.0  -1.75  15.00  -2.50  0.81e-3  0.74                            *
%      75.0  -0.50  14.50   2.50  0.62e-3  0.30 ];                         *
%      lat     dP     dT     dRH   dbeta  dlambda                          *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%OUTPUT:                                                                   *
%       AVG - Averaged meteorological parameters (nx5)                     *
%       AMP - Seasonal  variation/Amplitude meteorological parameters (nx5)*
%WHERE:                                                                    *
%      nx5-> Means n Rows,five(5) Columns (e.g. 3x5)                       *
%      n->Indicates the Number of User/Receiver Input Coordinates(Latitude)*
%          That is, each row in (AVG/AMP) represents USER coordinate       *
%          The Columns arrangements are as follows:                        *
%          1       2     3       4       5                                 *
%     AVG:[P       T     RH     beta   lambda]                             *
%     AMP:[P       T     RH     beta   lambda]                             *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%==========================================================================+
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************
%***************************************************************************
%*********ASSIGNMENT
avg=avgValues;
amp=ampValues;

UserLat=abs(UserLat);%User Latitude Coordinates

%******INITIALIZE OUTPUT
AVG=zeros(length(UserLat),length(avg)-1);%Assigning zeros to AVG
AMP=deal(AVG);%copy the contents AVG to AMP

for i=1:length(UserLat) %Loop Over User given Latitude coordinates
    
    for j=1:size(avg,1) %Loop Over Look-up Table
        
       if (j==1 & UserLat(i)<= avg(j,1))
           
         AVG(i,:) = avg(j,2:end);
         AMP(i,:) =amp(j,2:end);
           
       elseif  (j==size(avg,1) & UserLat(i)>= avg(j,1))
           
             %*****AVERAGE SURFACE TROPO VALUEs
              AVG(i,:) = avg(j,2:end);
              AMP(i,:) =amp(j,2:end);
              
       else       
            if (UserLat(i)== avg(j,1))
                
              AVG(i,:) = avg(j,2:end);
              AMP(i,:) =amp(j,2:end);
                               
            elseif  (UserLat(i)>avg(j,1) && avg(j,1)<avg(j+1,1))
                  
                  %********COMPUTE SURFACE TROPO VALUEs BY INTERPOLATION
                  %***COMPUTE DIFFERENCEs
                  Davg= avg(j+1,2:end)-avg(j,2:end);%Difference in average Met parameters... 
                                                    %between two closest latitudes .
                  Damp= amp(j+1,2:end)-amp(j,2:end);%Difference in seasonal variation of Met ... 
                                                    %parameters between two closest latitudes.
                  DlatUT=UserLat(i)-avg(j,1);%Difference in latitude values between user ...
                                             %Lat and closest latitude of Met parameters.
                  DlatT=avg(j+1,1)-avg(j,1);%Difference in latitude values in the look-up Table                       
                
                  %****PERFORM INTERPOLATION
                  AVG(i,:) = ((DlatUT/DlatT).*Davg)+avg(j,2:end); %Average Meteorological parameters
                  AMP(i,:) = ((DlatUT/DlatT).*Damp)+ amp(j,2:end);%Seasonal variation Meteorological parameters
                                            
            end %if (UserLat(i)== avg(j,1))
            
       end %if (j==1 & latR(i)<= avg(j,1))
       
       if (j+1)>size(avg,1)
             
          break
       end  %if j+1>length(TableLat)
       
    end %for j=1:size(avg,1)
    
end %for i=1:length(UserLat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF interpolate_UNB3m.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%6.SUBROUTINE TO COMPUTE JULIAN,MODIFIED DAY(JD,MJD) &  DAY OF YEAR(DoY)
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
 
%6.1.*********SUBROUTINE TO COMPUTE DAY OF YEAR(DoY)
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
     jDay = JulianDay(Yr,Mn,D,0,0,0);%Midnight Morning of this Day
     jDayJan1   = JulianDay(Yr,1,1,0,0,0);%Midnight Morning of January 1st
     dayNumber = jDay - jDayJan1 + 1;
     DoY=dayNumber;%Day of the Year     
end %try
%******************************************END OF DayofYear.m 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%6.2*********SUBROUTINE TO CONVERT TWO DIGITs TO YEAR 4 DIGITS
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

%7.*******SUBROUTINE TO COMPUTE GEOID UNDULATION & ORTHOMETRIC HEIGHT
function [undu,horth] = GETundulation(lat,lon,h)

%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            This subroutine determines Geoid Undulation and computes      * 
%            Orthometric Height using either goGPS Geoid model(EGM2008) OR *
%            MATLAB EGM96 Geopotential Model                               *                          
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%      lat   = Latitude of point[degrees]                                  *
%      lon   = Longitude of point[degrees]                                 *
%      h     = Ellipsoidal height of point (m)                             *

%OUTPUT:                                                                   *
%       undu  = Geoid undulation                                           *
%       horth = Orthometric Height                                         *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%--------------------------------------------------------------------------
%USE goGPS GEOID MODEL(EGM2008 Geopotential Model) AS DEFAULT  GEOID MODEL
%--------------------------------------------------------------------------  
%RETRIEVE STORED GEOID MODEL
 geoid = getappdata(0,'geoid');
 
 if isempty(geoid)%if geoid model was not stored
   
    %LOAD IT FROM goGPS geoid folder
     
   try %goGPS v0.5.2 beta1
      gs = Go_State.getInstance;
      geoid = gs.getRefGeoid(); 
      
   catch %goGPS v0.4.3 
         try
              try
                 %LOAD GEOID MODEL
                 load ../data/geoid/geoid_EGM2008_05.mat  %goGPS v0.4.3
              catch
                   load ../data/reference/geoid/geoid_EGM2008_05.mat %goGPS v0.5.2 beta1
              end
    
             %geoid grid and parameters
             geoid.grid = N_05x05;
             geoid.cellsize = 0.5;
             geoid.Xll = -180;
             geoid.Yll = -90;
             geoid.ncols = 720;
             geoid.nrows = 360;

             clear N_05x05
         catch 
              %geoid unavailable
              geoid.grid = 0;
              geoid.cellsize = 0;
              geoid.Xll = 0;
              geoid.Yll = 0;
              geoid.ncols = 0;
              geoid.nrows = 0;
         end 
         
   end %try
    
 end %//if isempty(geoid)
     
%**************COMPUTE GEOID UNDULATION AND ORTHOMETRIC HEIGHT
try    
   if (exist('geoid','var') && isfield(geoid,'ncols') && geoid.ncols ~= 0)
         
      %GEOID UNDULATION INTERPOLATION
      undu = grid_bilin_interp(lon, lat, geoid.grid, geoid.ncols, geoid.nrows, geoid.cellsize, geoid.Xll, geoid.Yll, -9999);
      horth   = h - undu;%Orthometric height
        
   else %IF goGPS GEOID MODEL FAILS,USE MATLAB EGM96 Geopotential Model 
        undu = (geoidheight(lat,lon))';%Undulation
        horth   = h - undu; %Orthometric height 
   end 

catch %IF ANYTHING GOES WRONG,USE MATLAB EGM96 Geopotential Model 
      undu = (geoidheight(lat,lon))';%Undulation
      horth   = h - undu; %Orthometric height 
           
end
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GETundulation.m 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  

%7.1.SUB-ROUTINE TO INTERPOLATE(BILINEAR) GEOID GRIDS
function [N] = grid_bilin_interp(lon, lat, grid, ncols,...
                                            nrows, cellsize, Xll, Yll, nodata)

%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            The Function applies a bilinear interpolation of the four     * 
%            nearest nodes of a georeferenced grid in correspondence of a  * 
%            point of given coordinates.                                   *
%--------------------------------------------------------------------------*
% SYNTAX:                                                                  *
%      [N] = grid_bilin_interp(lon, lat, grid, ncols,nrows,cellsize,Xll,...*
%                                                            Yll, nodata); *
%--------------------------------------------------------------------------*
% INPUT:
%       lon = X(longitude) coordinate of the interpolation point
%       lat = Y(latitude) coordinate of the interpolation point
%       grid = matrix containing the grid
%       ncols = number of columns of the grid
%       nrows = number of rows of the grid
%       cellsize = ground size of a cell
%       Xll = X coordinate of the center of the lower left cell
%       Yll = Y coordinate of the center of the lower left cell
%       nodata = value used for cells not containing data
%
% OUTPUT:
%        N = interpolated value [which is the geoid undulation]            *
%
%--------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%--------------------------------------------------------------------------
%Initialize/Preallocate output(N-undulation) with zeros
n = size(lat,1);
N = zeros(n,1);

%Loop through all of the lon/lat input pairs
for i = 1 : size(lat,1)
    
     Lon = lon(i);
     Lat = lat(i);
     
     %preparation of the grid axes
     X = (Xll : cellsize : Xll + (ncols - 1) * cellsize)';
     Y = (Yll : cellsize : Yll + (nrows - 1) * cellsize)';

     if (Lon <= X(1) | Lon >= X(end) | Lat <= Y(1) | Lat >= Y(end))
        interp_value = nodata; %#ok<*NASGU>
        return
     end 
    
     %detection of the grid node nearest to the interpolation point
     [mX, posX] = min(abs(X - Lon));
     [mY, posY] = min(abs(Y - Lat));

     %definition of the four grid nodes that sorround the interpolation point
     % (i,j) image coordinates (upper-left origin)
     % (X,Y) ground coordinates (bottom-left origin)
     if (X(posX) > Lon) | (mX ==0)
        j_left = posX - 1;
        j_right = posX;
        X_left = X(posX - 1);
     else 
         j_left = posX;
         j_right = posX + 1;
         X_left = X(posX);
     end

     if (Y(posY) > Lat) | (mY ==0)
        i_up = nrows + 1 - posY;
        i_down = i_up + 1;
        Y_down = Y(posY - 1);
     else 
         i_down = nrows + 1 - posY;
         i_up = i_down - 1;
         Y_down = Y(posY);
     end 

     %if one of the interp_value values of the four sorrounding points is a nodata value, do not interpolate and return nodata value
     if (grid(i_up,j_left) == nodata | grid(i_up,j_right) == nodata | grid(i_down,j_left) == nodata | grid(i_down,j_right) == nodata)
        interp_value = nodata;
        return
     end 

     %computation of the parameters of the bilinear function
     %f(X, Y) = a*X*Y + b*X + c*Y + d

     A = [0 0 cellsize 1; ...
         cellsize^2 cellsize cellsize 1; ...
         0 0 0 1; ...
         0 cellsize 0 1];

     B = [grid(i_up,j_left);...
          grid(i_up,j_right);...
          grid(i_down,j_left);...
          grid(i_down,j_right)];

     bilin_param = A\B;

     i_approx = Lat - Y_down;
     j_approx = Lon - X_left;

     %computation of the interpolated value
     interp_value = bilin_param(1) * j_approx * i_approx + bilin_param(2) * j_approx + bilin_param(3) * i_approx + bilin_param(4);
     interp_value = double(interp_value);

     N(i,1) = interp_value;%Geoid Undulation

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF grid_bilin_interp.m 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-   


%8.SUBROUTINE TO COMPUTE BLACK & EISNER MAPPING FUNCTION
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [MFh,MFw]=Black_Eisner_MF(satELEV)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%           "MF_Black_Eisner" is a subroutine that Computes the wet and dry*
%            Tropospheric mapping function Using Black & Eisner(1984)      *            
%            mapping function                                              *

%USAGE:                                                                    *                                                           *
%      [MFh,MFw]=Black_Eisner_MF(satELEV)                                  *
%INPUT:                                                                    *
%      satELEV : Satellite Elevation in [degrees]                          *

%OUTPUT:                                                                   *
%1.     MFh ==> Hydrostatic/dry mapping function                           *
%2.     MFw ==> Non Hydrostatic/wet mapping function                       *
%==========================================================================
%REFERENCE:                                                                * 
%1.        Black, H. D., & Eisner, A. (1984). Correcting satellite Doppler *
%                                          data for tropospheric effects   *
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************

%1.*****INITIALIZING OUTPUT VARIABLEs
satEL=zeros(size(satELEV,1),size(satELEV,2));%Assign zeros of nxm to SHD

[MFh,MFw]=deal(satEL);%copy the contents of satEl to all the requested outputs

%2.**************COMPUTE MAPPING FUNCTIONs(MFh & MFw)

%2.1.FIND NUMBER OF ROWS & COLUMNs IN satELEV
[nrow,ncol]=size(satELEV);

%2.2.LOOP OVER USER POSITIONs & SATELLITE ELEVATIONs
for i=1:ncol %Loop over the Number of Receiver positions
    
   for j=1:nrow %Loop over the Number of Satellite Elevations; 
                
       %2.2.1.CONVERT ELEVATION ANGLEs TO RADIAN  
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       %First, Check if the elevation is less than 0, set it to .1 deg
       %REMARKs:
       %       Valid elevations are between -pi/2 and pi/2. Elevations below 
       %       0.1 will have the a delay mapped to 0.1 degree.No mapping  
       %       below zero degree elevation is modeled.
       
       EL_zero = find(satELEV(j,i) < .0017);
       if ~isempty(EL_zero)
          satELEV(EL_zero) = ones(size(satELEV(EL_zero)))*.0017;
       end   % if ~isempty(EL_zero)
       
       %CONVERT ELEVATION ANGLE IN DEGREES TO RADIANS
       satEL(j,i) = abs(satELEV(j,i)) * pi / 180;  
      
        %NOW THE MAPPING FUNCTIONs(MF)
      try
         MFh(j,i)=1.001./(sqrt(0.002001+sin(satEL(j,i)).^2));%Hydrostatic Mapping Function
         MFw(j,i)=1.001./(sqrt(0.002001+sin(satEL(j,i)).^2));%Wet Mapping Function         
      catch 
           MFh(j,i)=1./(sqrt(1-((cos(satEL(j,i))./1.001).^2)));%Hydrostatic Mapping Function
           MFw(j,i)=1./(sqrt(1-((cos(satEL(j,i))./1.001).^2)));%Wet Mapping Function
      end %\\try  

   end   %for j=1:nrow

end  %for i=1:ncol

%%%%%%%%%%%%%%%%%%%%%%%%END OF MF_Black_Eisner.m %%%%%%%%%%%%%%%%%%%%%%% 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-