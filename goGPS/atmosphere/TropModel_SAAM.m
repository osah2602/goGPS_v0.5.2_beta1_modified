%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            "TropModel_SAAM" is a subroutine that Computes the wet,dry    *
%            and Total Tropospheric Delays Using Saastamoinen Refined Model* 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%USAGE:                                                                    *
%      General:                                                            *
%             [STD,SHD,SWD,ZTD,ZHD,ZWD]=TropModel_SAAM(ReceiverPos,SatPos,.*
%                                      Temperature,Pressure,WaterVaporPres)*
%      Others:                                                             *
%I.         [STD,SHD,SWD,ZTD,ZHD,ZWD]=TropModel_SAAM(ReceiverPos,SatPos,...*
%                                                                [T  P  e])*
%II.        [STD,SHD,SWD,ZTD,ZHD,ZWD]=TropModel_SAAM(ReceiverPos,SatPos,...*
%                                                  Temperature,Pressure,[])*
%III.       [STD,SHD,SWD,ZTD,ZHD,ZWD]=TropModel_SAAM(ReceiverPos,SatPos,...*
%                                                   Temperature,[],[])     *
%IV.        [STD,SHD,SWD,ZTD,ZHD,ZWD]=TropModel_SAAM(ReceiverPos,SatPos,...
%                                                                [],[],[]) *   
%V.         [STD,SHD,SWD,ZTD,ZHD,ZWD]=TropModel_SAAM(ReceiverPos,SatPos)   *                     
%VI.        [STD,SHD,SWD,ZTD,ZHD,ZWD]=TropModel_SAAM(height,satElevation,..*
%                                      Temperature,Pressure,WaterVaporPres)*
%VII.       [STD,SHD,SWD,ZTD,ZHD,ZWD]=TropModel_SAAM(height,satElevation,..*
%                                                                 [],[],[])*
%VIII.      [STD,SHD,SWD,ZTD,ZHD,ZWD]=TropModel_SAAM(height,satElevation)  *
%IX.        [STD,SHD,SWD,ZTD,ZHD,ZWD]=TropModel_SAAM([],satElevation)      *
%X.         [STD,SHD,SWD,ZTD,ZHD,ZWD]=TropModel_SAAM([],satElevation,...   *
%                                      Temperature,Pressure,WaterVaporPres)*
%XI.        [STD,SHD,SWD,ZTD,ZHD,ZWD]=TropModel_SAAM([],[],Temperature...  *
%                                                ,Pressure,WaterVaporPres) *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%**THE FUNCTION CALLS THE FOLLOWING SUBROUTINE(s):                         *
%1.[X,Y,Z] = geo2xyz(latitude,longitude,height,RefEllipsoid);              *
%2.[latRAD,longRAD,h,latDEC,longDEC] = xyz2LLH(X,Y,Z,RefEllipsoid);        *
%3.[a,finv] = Elipsoidpara(RefEllipsoid);                                  *
%4.[Az_rad,El_rad,Az_deg,El_deg,D]=satAzimuthElevation(UserXYZ,SatXYZ,...  *
%                                                            RefEllipsoid);*
%5.[ZTD, ZHD, ZWD] =ZenithTropDelay(T,P, es);                              *
%6.[STD, SHD, SWD] =SlantTropDelay(T,P,es,satELEV);                        *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%      Generally, the function accepts five(5) sets of inputs:             *
%A.    POSITIONs:                                                          *
%      =-=-=-=-=                                                           *
%1.    ReceiverPos--> Receiver position in either Latitude,Longitude &     *
%                     height or ECEF(XYZ) Coordinates                      *
%2.    SatPos  -----> Satellite position(s) in ECEF(XYZ) Coordinates       *

%B.    METEOROLOGICAL PARAMETER:                                           *
%      =-=-=-=-=-=-=-=-=-=-=-=-=                                           *
%3.    Temperature(T)--> Atmospheric Temperature in Degree Celsius(C)      *
%4.    Pressure(P)  ---> Atmospheric Pressure in millibars(mbar /hPa)      *
%5.    WaterVaporPres(e)--> Water Vapor Partial Pressure in [mbar]         *

%OTHER CONSIDERATIONs:                                                     *
%*******Any missing Meteorological parameter, standard Meteorological value* 
%       is/are assigned as default value(s).i.e: ....                      *
%       [P: 1013.25(mbar), T: 291.16(K), RH:50(%)]                         *

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
% =========================================================================
%REFERENCE:                                                                *
%1.        Bosy J.,Borkowski J.:Troposphere Modeling in Local GPS Network  * 
%2.        J. Sanz Subirana,et al., GNSS Data Processing, Vol. I:          *
%          Fundamentals and Algorithms(ESA TM-23/1, May 2013)              *
%3.        "GPS Theory and application",edited by B.Parkinson,J.Spilker,   *
%                                                         P.Enge, AIAA,1996*
%4.        GPS Theory and Practice,Hofmann-Wellenhof, Lichtenegger, and... *
%                            Collins,Fifth, revised edition  pages 106-115.* 
%5         "Global Positioning System, Mishra & Enge", pg 172              *
%6.        Modeling of Tropospheric Delays Using ANFIS,Wayan Suparta;...   *
%                                                Kemal Maulana Alhasa(2016)*
%7.        GPS Theory,Algorithms and Applications 2ed,Guochang Xu,June 2007*
%                                                 pg 55-61 eqns(5.95-5.108)*
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************
%***************************************************************************

function [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_SAAM(ReceiverPos,SatPos,Temperature,...
                                                   Pressure,WaterVaporPres)
%**********CHECK INPUTs & REFORMAT INPUT DATA
switch nargin
    
    case {5,4,3,2,1} %When all inputs are provided        
        
        if (any(nargin==[5,4,3,2,1]))
            
          if nargin ==5 %for five(5) or all inputs...              
            %******Assignments   
            T=Temperature;%Atmospheric Temperature in Degree Celcius
            P=Pressure;%Atmospheric Pressure in millibars(mbar)/Hecto pascal(hpa)
            es=WaterVaporPres;%partial water vapor pressure(e) 
            
          elseif  nargin ==4 %For four(4) inputs...                
                %******Assignments               
                T=Temperature;%Atmospheric Temperature in Degree Celcius
                P=Pressure;%Atmospheric Pressure in millibars(mbar)/Hecto pascal(hpa)
                es=[];%Assign empty matrix([]) to partial water vapor pressure(e)  
                
          elseif nargin==3 %For three(3) inputs...                              
                %CHECK THE NUMBER OF COLUMNs IN Temperature
                Tcol=size(Temperature,2); 
                
                switch Tcol 
                    
                    case 1 %if num of col is 1
                        T=Temperature;%Atmospheric Temperature in Degree Celcius
                        P=[];%Assign empty matrix([]) to Atmospheric Pressure 
                        es=[];%Assign empty matrix([]) to Relative Humidity 
                        
                    case 2 %if num of col is 2
                        T=Temperature(:,1);%Assign 1st Column of Temperature
                                            %Temperature in Degree Celcius
                        P=Temperature(:,2);%Assign 2nd Column of Temperature
                                           %        to Atmospheric Pressure 
                        es=[];%Assign empty matrix([]) to partial water vapor pressure(e)  
                        
                    case 3 %if num of col is 3
                        T=Temperature(:,1);%Assign 1st Column of Temperature
                                            %Temperature in Degree Celcius
                        P=Temperature(:,2);%Assign 2nd Column of Temperature
                                           %        to Atmospheric Pressure 
                        es=Temperature(:,3);%Assign 3rd Column of Temperature
                                            %to partial water vapor pressure(e)                                            
                    case 0
                         T=[];%Assign empty matrix([]) to Temperature in Degree Celcius
                         P=[];%Assign empty matrix([]) to Atmospheric Pressure 
                        es=[];%Assign empty matrix([]) to partial water vapor pressure(e)
                          
                    otherwise
                             T=Temperature(:,1);%Assign 1st Column of Temperature
                                            %Temperature in Degree Celcius
                             P=Temperature(:,2);%Assign 2nd Column of Temperature
                                           %        to Atmospheric Pressure 
                             es=Temperature(:,3);%Assign 3rd Column of Temperature
                                            %to partial water vapor pressure(e)                    
                end %switch Tcol  
                
          elseif  nargin==2 %For two(2) inputs...              
                %*****Assignments
                T=[];%Assign empty matrix([]) to Temperature in Degree Celcius
                P=[];%Assign empty matrix([]) to Atmospheric Pressure 
                es=[];%Assign empty matrix([]) to partial water vapor pressure(e) 
                
          elseif  nargin==1 %for single input            
                %*****Assignments                
                T=[];%Assign empty matrix([]) to Temperature in Degree Celcius
                P=[];%Assign empty matrix([]) to Atmospheric Pressure 
                es=[];%Assign empty matrix([]) to partial water vapor pressure(e)      
            SatPos=[];%Assign empty matrix([]) to SatPos
                                               
          end  %if nargin ==5
          
        end  %if (any(nargin==[5,4,3,2]))              
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
%************TROPOHERIC DELAY MODELING USING SAASTAMOINEN MODEL(modified)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%**********CHECK FOR MIN & MAX HEIGHTS(0 <= h <= 5km 
%Saastamoinen model requires positive ellipsoidal height

if  exist('h','var') %Check if user provided heights
    
  if all(h(:,1)~=0 ) %if heights are not zeros
   
    h(h < 0) = 0;
    h(h > 5000) = 5000;%Ellipsoidal height should not be > 5000m(5km)
    
  end
  
else
     h = 0;%assign zero(0) to heights if not provided      
end
  
%***************COMPUTE ANY MISSING METEOTOLOGICAL PARAMETER(S)  
%********GET METEOROLOGICAL PARAMETERS
if all([isempty(T) isempty(P) isempty(es)])
   %Use Standard atmosphere @ MSL(h=0)- Berg, 1948 (Bernese)
   Ps  = 1013.25;%pressure [mbar] @ SEALEVEL
   Ts  = 291.15;%temperature [K] @ SEALEVEL
   RHs = 50.0;%Relative Humidity(%)@ SEALEVEL
    
   %***COMPUTE ATMOSPHERIC PARAMETERS FROM STANDARD PARAMETERS
   P  = Ps .* (1-0.0000226.*h).^5.225;%pressure [mbar] @ altitude(h)
   T  = Ts - 0.0065.*h; %temperature [K] @ altitude(h)
   RH = RHs * exp(-0.0006396.*h);%Relative Humidity(%)@ altitude(h)
    
elseif  (all([~isempty(T) isempty(P) isempty(es)]))
         T   = T + 273.15; %Convert temperture in Celcius(�C) to kelvin(K)
         Ps  = 1013.25;%pressure [mbar] @ SEALEVEL
	     RHs = 50.0;%Relative Humidity(RH)[%] @ SEALEVEL

         %***COMPUTE PRESSURE & RH 
         P  = Ps .* (1-0.0000226.*h).^5.225; %pressure [mbar] @ altitude(h)
	     RH = RHs * exp(-0.0006396.*h); %Relative Humidity(%)@ altitude(h)
     
elseif  (all([~isempty(T) ~isempty(P) isempty(es)]))
        T   = T + 273.15; %Convert temperture in Celcius(�C) to kelvin(K)
        RHs = 50.0;%Relative Humidity(RH)[%] @ SEALEVEL
      
        %***COMPUTE PRESSURE & RH 
	    RH = RHs * exp(-0.0006396.*h); %Relative Humidity(%)@ altitude(h)          
else   
      T = T + 273.15; %Convert temperture in Celcius(�C) to kelvin(K)  
    
end    %if (all([~isempty(T) isempty(P) isempty(es)]))
       
%COMPUTE PARTIAL PRESSURE OF WATER VAPOR(e)
if exist('RH','var')
try    
    es = 0.01 .* RH .* exp(-37.2465 + 0.213166.*T - 0.000256908.*T.^2);
catch
    try
       es = 6.11.*(RH/100).*10.^(7.5.*(T-273.15)./(237.3+(T-273.15)));
    catch
        try
           es = 6.108.*(RH./100).*exp((17.15.*T-4684)./(T-38.45));
        catch
            try
                es = 6.107.*(RH./100).*exp((19.82.*T-5414.03)./T);
            catch
                 %Clausius-Clapeyron eqn
                 es = (RH/100).*6.11.*exp((2.5e6/461).*((1/273.15)-(1./T)));
            end
        end
    end
end

end %if exist('RH','var')

%**********COMPUTE TOTAL TROPOSPHERIC DELAY IN THE ZENITH  DIRECTION
%Call the "ZenithTropDelay.m" function
[ZTD, ZHD, ZWD] = ZenithTropDelay(T,P,es);

%**********COMPUTE TOTAL TROPOSPHERIC DELAY IN THE SATELLITE DIRECTION(LOS)
 if exist('satEL_deg','var')
     
%Call the "SlantTropDelay.m" function
[STD, SHD, SWD] = SlantTropDelay(T,P,es,satEL_deg);

 else
     %Create zero matrix 
     STD = zeros(size(ZTD));
     SHD = zeros(size(ZTD));
     SWD = zeros(size(ZTD));  
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF TropModel_SAAM.m 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%A.******SUBROUTINE TO COMPUTE ZENITH TROPOSPHERIC DELAY
function [ZTD, ZHD, ZWD] = ZenithTropDelay(T,P,es)

%DESCRIPTION:
%            "ZenithTropDelay" Computes Zenith Tropospheric Delay by ...   * 
%             Saastamoinen model.                                          *
%******SYNTAX:
%             [ZTD, ZHD, ZWD] = ZenithTropDelay(T,P, es);                  *
%
%******INPUT:
%            P = atmospheric pressure [mbar/hPa]                           *
%            T = atmospheric temperature in degrees Kelvin [k]                               *
%            es = partial pressure of water vapor [mbar/hPa]               *
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

%*****COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
ZHD = (0.002277 .* P)./dgref;  %Zenith Hydrostatic or dry delay (m)

%*****COMPUTE ZENITH WET DELAY(ZWD)
ZWD = (0.002277.*((1255 ./ T) + 0.05).* es)./dgref;% Zenith wet delay

%****FINALLY, COMPUTE ZENITH TOTAL DELAY(ZTD)
try
   ZTD = ZHD + ZWD;%Zenith Total Delay (m)
catch    
       ZTD = (0.002277 .*(P + ((1255 ./ T) + 0.05).* es))./dgref; %Zenith  Total Delay
end 

%***********************END OF ZenithTropDelay.m ***********************    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%B.SUBROUTINE TO COMPUTE SLANT TROPOSPHERIC DELAY (SATELLITE IN VIEW)
function [STD, SHD, SWD] = SlantTropDelay(T,P,es,satELEV)

%DESCRIPTION:
%            "SlantTropDelay" Computes Slant Tropospheric Delay by ...     * 
%             Saastamoinen model.                                          *
%******SYNTAX:
%             [STD, SHD, SWD] = SlantTropDelay(T,P,e,satELEV);             *
%
%******INPUT:
%            P = atmospheric pressure [mbar/hPa]                           *
%            T = atmospheric temperature [k]                               *
%            es = Water Vapor Partial Pressure [mbar/hPa]                  *
%      satELEV = Satellite Elevation Angle in degrees

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
%*****RETRIEVE STORED DATA
latD = getappdata(0,'latD'); %Get station latitude Coords [degrees]
lonD = getappdata(0,'lonD'); %Get station latitude Coords [degrees]
latR = getappdata(0,'latR'); %Get station latitude Coords [radians]
   h = getappdata(0,'ht'); %Get station height [meters]
   
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

%**********COMPUTE ZENITH DELAYs
%Call the "ZenithTropDelay.m" function
[~, ZHD, ZWD] = ZenithTropDelay(T,P,es);

%FIND NUMBER OF ROWS & COLUMNs IN satELEV
[nrow,ncol] = size(satELEV);

for i = 1:ncol %Loop over the Number of Receiver positions
    
   for j = 1:nrow %Loop over the Number of Satellite Elevations
    
      %1.COMPUTE ZENITH ANGLEs(Z) FROM ELEVATION ANGLEs & CONVERT TO RADIAN
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      %First, Check if the elevation is less than 0, set it to .1 deg
      %REMARKs:
      %       Valid elevations are between -pi/2 and pi/2. Elevations below 
      %       0.1 will have the a delay mapped to 0.1 degree.No mapping  
      %       below zero degree elevation is modeled.
						
      EL_zero = find(satELEV(j,i) < .0017);
      if ~isempty(EL_zero)
        satELEV(EL_zero) = ones(size(satELEV(EL_zero)))*.0017;
      end  % if ~isempty(EL_zero)
      
      Za(j,i) = (90-satELEV(j,i));%Zenith Angle(s)/distance in degrees(deg)
      Zd(j,i) = (90-satELEV(j,i)).* pi / 180; %Zenith Angle(s)/distance converted to radian

      %*****INTERPOLATE CORRECTION FACTORS(B,deltaR)
      %Call the "interpolate_Cfactors.m" function 
      [B, deltaR] = interpolate_Cfactors(horth(i,1),Za(j,i));
      
      %*****COMPUTE MAPPING FUNCTIONs
      MF1(j,i) = 1./cos(Zd(j,i));%Cosecant MF
      MF2(j,i) = -0.002277.*B.*(tan(Zd(j,i)).^2./cos(Zd(j,i)));
      MF3(j,i) = deltaR;
      
      %*****COMPUTE SLANT HYDROSTATIC DELAY(shd)      
      SHD(j,i) = ZHD(i).*MF1(j,i); %(m)
            
      %****COMPUTE SLANT TOTAL DELAY(STD)
      try
         STD(j,i) = MF1(j,i).*(ZHD(i)+ZWD(i))+ MF2(j,i)+MF3(j,i);      
      catch
           STD(j,i) = (((0.002277./cos(Zd(j,i))) .* (P(i,1) + ((1255 ./ T(i,1)) + 0.05).* es(i,1)- B.*(tan(Zd(j,i))).^2 ))./dgref(i,1))+deltaR; %(m)
      end
      
      %***********COMPUTE SLANT WET DELAY(swd)
      try
         SWD(j,i) = STD(j,i)- SHD(j,i);
      catch
           SWD(j,i) = (((0.002277./cos(Zd(j,i))) .* (((1255 ./ T(i,1)) + 0.05).* es(i,1)- B.*(tan(Zd(j,i))).^2 ))./dgref(i,1))+deltaR;%(m)
      end
                                   
   end %for j=1:nrow
          
end %for i=1:ncol
    
%***********************END OF SlantTropDelay.m ***********************    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
%%%%%%%%%%%%%%%%%%%%%%%%END OF TropModel_SAASM.m %%%%%%%%%%%%%%%%%%%%%%% 

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

%1.******SUBROUTINE TO INTERPOLATE SAAS CORRECTION FACTORS
function [B, deltaR] = interpolate_Cfactors(Height,ZenithAngle)
% *************************************************************************
% *************************************************************************
%DESCRIPTION:                                                              *
%            "interpolate_Cfactors" is a subroutine that  Linearly....
%             interpolates Saastamoinen Refined Model Correction factors...
%             of B and deltaR given user Height(s) and Zenith Angle(s)
%             i.e.B is a function of height [B(h)] and deltaR is a 
%             function of height and Zenith Angle[deltaR(Z,h)]
%NOTE:
%     Values of B and deltaR  Correction terms for Height Zenith Angle 
%     respectively are given in a table form by Saastamoinen

%USAGE:
%      [B, deltaR] = interpolate_Cfactors(Height,ZenithAngle)

%INPUT:
%1.    Height--> Station ellipsoidal height in meters(m)
%2.    ZenithAngle--> Zenith Angle of the Satellite in degrees(deg)

%OUTPUT:
%1.     B --> Correction term for Height in millibar(mb)
%2.     deltaR --> Correction term for Zenith Angle

% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%**************************************************************************
%**************************************************************************

%*********ASSIGNMENTs
uZA=ZenithAngle;
h=Height;

h(h < 0) = 0;%Ellipsoidal height should not be < 0
h(h>5000)=5000;%Ellipsoidal height should not be > 5000m(5km)

%*******TABLE VALUEs FOR correction terms B & deltaR
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%Ref:
%1.  GPS Theory and Practice,Hofmann-Wellenhof, Lichtenegger, and... 
%                            Collins,Fifth, revised edition  pages 106-115.
%2.  GPS Theory, Algorithms and Applications,Guochang Xu 2Ed ,pages 56-57.
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

 %*****HEIGHT CORRECTION VALUEs
Bc = [1.156;1.079;1.006;0.938;0.874;0.813;0.757;0.654;0.563];%Corrections Terms for height(mb)
hc = [0;500;1000;1500;2000;2500;3000;4000;5000];%Heights(m)  nx1 matrix

%****ZENITH ANGLE CORRECTION VALUEs 
tZA = [60;66;70;73;75;76;77;78;78.5;79;79.5;79.75;80];%Table Zenith Angle(deg) 
hcc = [0 500 1000 1500 2000 3000 4000 5000];%1xn matrix

Rc =[0.003,0.003,0.002,0.002,0.002,0.002,0.001,0.001;...
	0.006,0.006,0.005,0.005,0.004,0.003,0.003,0.002;...
	0.012,0.011,0.010,0.009,0.008,0.006,0.005,0.004;...
	0.020,0.018,0.017,0.015,0.013,0.011,0.009,0.007;...
	0.031,0.028,0.025,0.023,0.021,0.017,0.014,0.011;...
	0.039,0.035,0.032,0.029,0.026,0.021,0.017,0.014;...
	0.050,0.045,0.041,0.037,0.033,0.027,0.022,0.018;...
	0.065,0.059,0.054,0.049,0.044,0.036,0.030,0.024;...
	0.075,0.068,0.062,0.056,0.051,0.042,0.034,0.028;...
	0.087,0.079,0.072,0.065,0.059,0.049,0.040,0.033;...
	0.102,0.093,0.085,0.077,0.070,0.058,0.047,0.039;...
	0.111,0.101,0.092,0.083,0.076,0.063,0.052,0.043;...
	0.121,0.110,0.100,0.091,0.083,0.068,0.056,0.047];%Corrections Terms for ZA
   
%***********LINEAR INTERPOLATION FOR SAAS CORRECTION TERMS(B & deltaR)

%1.******INTERPOLATING FOR B

try
    %*********INITIALIZING OUTPUT VARIABLEs  
 	B = zeros(length(h),1);
    
     for i=1:length(h)
         
         for j=1:length(hc) 
        
            if h(i)==hc(j)
              B(i,1)=Bc(j);
            else 
                if (h(i)>hc(j) && h(i)<hc(j+1)) 
                  B(i,1) = Bc(j) + (((h(i) - hc(j))/(hc(j+1)-hc(j)))*(Bc(j+1)-Bc(j)));
                                                            
                end %if (h(ii)>hc(iii) && h(ii)<hc(iii+1))                                         
            end %h(ii)==hc(iii)
         
         if (j+1)>length(hc)
             
             break
         end %if j+1>length(hc)
         
         end  %for j=1:length(hc) 
     end %for i=1:length(h)
     
catch 
      %USING MATLAB IN-BUILT INTERPOLATION ALGORITHM (interp1)
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
      %USAGE:vq = interp1(x,v,xq) returns interpolated values of a 1-D     +
      %     function at specific query points using linear interpolation.  +
      %     Vector x contains the sample points, and v contains the corresponding 
      %     values, v(x). Vector xq contains the coordinates of the query points
      %     In this case x=hc, v=Bc, xq=h, interpolation method=linear     +
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
      B=interp1(hc,Bc,h,'linear');%Linear interpolation
      
end %try
     
%2.******INTERPOLATING FOR deltaR
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%NOTE:                                                                     +
%Cases considered in this interpolation                                    +
%Case1:When user given height equals(=)table height & User given Zenith Angle
%      equals(=) table Zenith Angle(ZA)
%Case2:When user given height equals(=)table height & User given Zenith Angle
%      is greater/less(<>) than table Zenith Angle                         +
%Case3:When user given height is greater/less(<>) than table height & ...  +
%      User given Zenith Angle(Z) equals(=) table Zenith Angle(ZA)         +
%Case4:When user given height is greater/less(<>) than table height & ...  +
%      User given Zenith Angle is greater/less(<>) than table Zenith Angle +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%*********INITIALIZING OUTPUT VARIABLEs  
deltaR=zeros(size(uZA,1),1);
Rc1=zeros(size(tZA,1),1);
Rc2=zeros(size(tZA,1),1);

for ii=1:length(uZA) %Loop through User given Zenith Angles
    
   for k=1:length(tZA)-1%Loop through table Zenith Angles
       
      for l=1:length(h)%Loop through user given heights
          
         for m=1:length(hcc)-1%Loop through table heights
             
             if uZA(ii)<60 || uZA(ii)>80
                 
                 deltaR(ii,1)= 0;
             end 
             
                if (h(l)==hcc(m) && (uZA(ii)==tZA(k)))%Case1
                  
                    deltaR(ii,1)= Rc(k,m);
                  
                elseif ( h(l)==hcc(m) && (uZA(ii)>tZA(k) && tZA(k)<tZA(k+1)) )%Case2
                      
                    deltaR(ii,1)=Rc(k,m)+(((uZA(ii) - tZA(k))/(tZA(k+1)-tZA(k)))*(Rc(k+1,m)-Rc(k,m)));
                      
                else   
                     if (h(l)>hcc(m) && (hcc(m)<hcc(m+1) && uZA(ii)==tZA(k)) )%Case3
                      
                         deltaR(ii,1) = Rc(k,m) + (((h(l) - hcc(m))/(hcc(m+1)-hcc(m)))*(Rc(k,m+1)-Rc(k,m)));
                       
                     elseif ((h(l)>hcc(m) && hcc(m)<hcc(m+1)) && (uZA(ii)>tZA(k)&& tZA(k)< tZA(k+1)))%Case4
                                                    
                           Rc1(k,1)=Rc(k,m)+ (((uZA(ii) - tZA(k))/(tZA(k+1)-tZA(k)))*(Rc(k+1,m)-Rc(k,m)));
                           Rc2(k,1)=Rc(k,m+1)+ (((uZA(ii) - tZA(k))/(tZA(k+1)-tZA(k)))*(Rc(k+1,m+1)-Rc(k,m+1)));
                          deltaR(ii,1)=Rc1(k,1) + (((h(ii) - hcc(m))/(hcc(m+1)-hcc(m)))*(Rc2(k,1)-Rc1(k,1)));
                          
                     end  %if (h(l)>hcc(m) && (hcc(m)<hcc(m+1) && uZA(ii)==tZA(k)) )
             
                end %\\f (h(l)==hcc(m) && uZA(ii)==tZA(k))
                                             
         end %\\for m=1:length(hcc)
         
      end %\\ for l=1:length(h)
      
   end %\\for k=1:length(tZA)
   
end %\\for ii=1:length(uZA)

%%%%%%%%%%%%%%%%%%%%%%END OF interpolate_Cfactors.m %%%%%%%%%%%%%%%%%%%%%%%

%2.******SUBROUTINE TO CONVERT XYZ COORDs TO GEOGRAPHIC COORDs
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
      
%3.******SUBROUTINE TO CONVERT XYZ COORDs TO GEOGRAPHIC COORDs
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
 
%4.******SUBROUTINE TO COMPUTE ELLIPSOIDAL PARAMETERs
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

%5.******SUBROUTINE TO COMPUTE SATELLITE AZIMUTH & ELEVATION ANGLE

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

%*******SUBROUTINE TO COMPUTE GEOID UNDULATION & ORTHOMETRIC HEIGHT
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
    %LOAD GEOID MODEL[%USING goGPS DEFAULT GEOID MODEL (EGM2008)]  
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
        
   else  %IF goGPS GEOID MODEL FAILS,USE MATLAB EGM96 Geopotential Model 
        undu = (geoidheight(lat,lon))';%Undulation
        horth   = h - undu; %Orthometric height 
   end 

catch %IF ANYTHING GOES WRONG,USE MATLAB EGM96 Geopotential Model 
      undu = (geoidheight(lat,lon))';%Undulation
      horth   = h - undu; %Orthometric height 
           
end
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GETundulatio.m 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
 
%SUB-ROUTINE TO INTERPOLATE(BILINEAR) GEOID GRIDS
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