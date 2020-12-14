function [pwv,Tm,PI] = compute_PWV(UTCtime,ReceiverPos,ZTD,ZHD)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            This function computes/retrieve GNSS Precipitable/Integrated  *  
%            Water Vapour (PWV/IWV) for input Zenith Total Delay(ZTD),...  *
%            Zenith Hydrostatic Delay(ZHD) & Optionally Zenith Wet Delay...*
%            (ZWD)and and station information(UTCtime,ReceiverPos)         *
% SYNTAX:                                                                  *
%        [pwv] = compute_PWV(UTCtime,ReceiverPos,ZTD,ZHD,ZWD)              *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
% INPUT:                                                                   *
%       Generally, the function accepts four(4) sets of inputs:            *
%1.    UTCtime : Receiver reception time in[year,Month,Day,Hour,Min,Sec]   *                                                     
%2.ReceiverPos : Receiver position of the form Latitude,Longitude & Height *
%                [lat lon h] where latitude & longitude coordinates can    *
%                either be in DMS or Decimal Degress(Deg) 7 Height in [m]. *
%3.        ZTD : Zenith Total Delay(ZTD) in [m]                            *
%4.        ZHD : Zenith Hydrostatic Delay(ZHD) in [m]                      *
%%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
%OUTPUT:                                                                   *
%       pwv = Precipitable or Integrated Water Vapour (PWV/IWV) in [ mm ]  *
%%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
%--- * --. --- --. .--. ... * ---------------------------------------------+
%               ___ ___ ___                                                +
%     __ _ ___ / __| _ | __                                                +
%    / _` / _ \ (_ |  _|__ \                                               +
%    \__, \___/\___|_| |___/                                               +
%    |___/                    v 0.4.3                                      +
%--------------------------------------------------------------------------+
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini                 *
% Portions of code contributed by Laboratorio di Geomatica, Polo Regionale * 
% di Como,Politecnico di Milano, Italy.                                    *
%--------------------------------------------------------------------------+
%Codes Written BY:
%      OSAH SAMUEL, MSC GEOMATIC ENGINEERING ==============================+
%      Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com
%      Phone:+233(0)246137410 / +233(0)509438484
%**************************************************************************
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%******PRELIMINARY

%(A)*************RETRIEVE ALL STORED DATA FROM gui_GPS.fig
%                -----------------------------------------
%GET SELECTED MEAN TEMPERATURE(Tm) MODEL
Tm_model = getappdata(0,'Tm_model');%selected Tm model as string
TmVAL    = getappdata(0,'TmVAL');%selected Tm model as value/numeric

%CHECK IF SELECTED Tm MODEL IS ANY OF THE GPT MODELS & GET GRID VALUES & RESOLUTION
if any([strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35),...
        any(TmVAL==[11,12])])
    
   %GET GRID VALUES & GRID RESOLUTION
   grid_GPT     = getappdata(0,'grid_TM') ;
   grid_res_GPT = getappdata(0,'grid_res_TM');

   %GET TIME VARIATION STATUS
   Timevar_GPT  = getappdata(0,'Timevar_Tm');
   
else %IF NOT ANY OF THE GPT MODELs
     %*****ASSIGN EMPTY([]) MATRIX
        grid_GPT       = [];%VMF GRID FILE
        grid_res_GPT   = []; %VMF GRID RESOLUTION
        Timevar_GPT    = 1;%set time variation to static
    
end

%CHECK IF SELECTED Tm MODEL IS GTVR-Tm model & GET GRID VALUES 
if any([strncmpi(Tm_model,'GTVR-Tm model [Yao et al 2018]',30),TmVAL == 9]) 
    
   %GET SAVED GTVR-Tm GRID VALUES
   gridV_GTVR = getappdata(0,'gridV_GTVR');

else
    gridV_GTVR = [];
    
end

%CHECK IF SELECTED Tm MODEL IS TVGG-Tm model & GET GRID VALUES 
if  any([strncmpi(Tm_model,'TVGG-Tm model [Jiang et al 2018]',32),TmVAL == 10])
         
    %GET SAVED TVGG-Tm model GRID VALUES
    gridV_TVGG = getappdata(0,'gridV_TVGG');
    
else
    gridV_TVGG = [];
    
end

%CHECK IF SELECTED Tm MODEL IS GTrop-Tm model & GET GRID VALUES / COEFF
if  any([strncmpi(Tm_model,'GTrop model [Sun et al 2019]',28),TmVAL == 15])
         
    %GET SAVED GTrop model GRID VALUES
    gridV_GTrop = getappdata(0,'gridV_GTrop');
    
else
    gridV_GTrop = [];
    
end


%******************************GET TEMPERATURE(T) 
%NOTE:
%     SOME OF THE MEAN TEMPERATURE (Tm) MODELS REQUIRE TEMPERATURE,T TO
%     PERFORM.FOR THOSE MODELS WE NEED TO GET TEMPERATURE AVLUES FROM THE SOURCE OF
%     METEOROLOGICAL PARAMETERS SELECTED TO COMPUTE THE TROPOSPHERIC DELAY
%     CORRECTIONS OR GET THEM FROM SELF GENERATING MET PARAMETER TROPO MODELS
%     SUCH AS THE GPT MODELS, UNB3m,EGNOS & MOPS MODELS
%--------------------------------------------------------------------------

%*********CHECK TEMPERATURE(TEMP) DEPENDENT MODELS
if any([strncmpi(Tm_model,'Bevis et al (1992) model',24),strncmpi(Tm_model,'Bevis et al (1995) model',24),strncmpi(Tm_model,'Mendes et al (2000) model',25),...
        strncmpi(Tm_model,'Schueler et al (2001) model',27),strncmpi(Tm_model,'Yao et al (2014) model',22),strncmpi(Tm_model,'TVGG-Tm model [Jiang et al 2018]',32)])


   %***********GET STORED STATION MET PARAMETER
   %IF USER OBSERVED/MEASURED TEMPERATURE AT GNSS SITE WITH MET DEVICE
   Ts = getappdata(0,'Ts');%get STATION TEMPERATURE METEOROLOGICAL SENSOR

   if any([~isempty(Ts),all(Ts(:) ~=0),any(find(Ts(:) ~=0)),any(Ts(:) ~=0)])
    
      T = Ts; %ASSIGNMENT
   
   else 

        %(2.0)GET OPTION FOR TROPO MODELLING
        option_Cmodels = getappdata(0,'option_Cmodels');%option button for combine tropo models
        option_Smodels = getappdata(0,'option_Smodels');%option button for separate tropo models

        %2.*******************GET TROPOSPHERIC DELAY CORRECTION MODELs
        %                     ========================================
        %(1.0)**********COMBINE TROPO MODELS
        if option_Cmodels == 1 %if Combine model option is selected [Dry + Wet Delay Models]
   
          %GET SELECTED MODEL 
          Tropo_Cmodel   = getappdata(0,'Tropo_Cmodel');%get selected item (string)from pop_lTropo_combine
          Tropo_Cmodel_v = getappdata(0,'Tropo_Cmodel_v');%get selected item (Numeric)from pop_lTropo_combine
   
          %CHECK IF SELECTED MODEL IS ONE OF THE GLOBAL PRESSURE & TEMPERATURE(GPT) MODELS
          if any([any([strncmpi(Tropo_Cmodel,'GPT2 (5° x 5°)',14),strncmpi(Tropo_Cmodel,'GPT2w (1° x 1°)',14),strncmpi(Tropo_Cmodel,'GPT2w (5° x 5°)',14),...
                       strncmpi(Tropo_Cmodel,'GPT3 (1° x 1°)',14),strncmpi(Tropo_Cmodel,'GPT3 (5° x 5°)',14)]),any(Tropo_Cmodel_v == [12,13,14,15,16])]) 
      
             T = getappdata(0,'T_GPT');%Retrieve stored T values
    
              %IF SELECTED MODEL UNB3m MODEL  
          elseif strncmpi(Tropo_Cmodel,'UNB3m',5)%if UNB3m
          
                 T = getappdata(0,'T_UNB');%Retrieve stored T values
          
                 %CONVERT T in Kelvin TO degree Celsius[°C]
                 %NOTE:UNB3m Model T OUTPUT IS IN degree Kelvin
                 T = T - 273.15;
     
               %IF SELECTED MODEL EGNOS MODEL       
          elseif strncmpi(Tropo_Cmodel,'EGNOS',5) %if EGNOS
           
                 T = getappdata(0,'T_EGNOS');%Retrieve stored T values
   
                 %CONVERT T in Kelvin TO degree Celsius[°C]
                 %NOTE:EGNOS Model T OUTPUT IS IN degree Kelvin
                 T = T - 273.15;
              
               %IF SELECTED MODEL MOPS MODEL 
          elseif strncmpi(Tropo_Cmodel,'MOPS',4) %if MOPS
          
                 T = getappdata(0,'T_MOPS');%Retrieve stored T values
          
                 %CONVERT T in Kelvin TO degree Celsius[°C]
                 %NOTE:MOPS Model T OUTPUT IS IN degree Kelvin
                 T = T - 273.15;
          
               %IF SELECTED MODEL IS THE VMF ZENITH TROPOSPHERIC DELAY(ZPD)       
          elseif any([strncmpi(Tropo_Cmodel,'VMF gridded ZTD',15),any([strncmpi(Tropo_Cmodel,'GTrop [Sun et al 2019]',22),...
                                                                strfind(Tropo_Cmodel,'GTrop')])])%if VMF ZTD OR GTrop
          
                 %COMPUTE SURFACE TEMPERATURE USING GPT MODEL IF TIME & STATION 
                 %POSITION INPUTS ARE NOT EMPTY=============================
                 if all([~isempty(UTCtime),~isempty(ReceiverPos)])
      
                    if exist('getMETpara_GPT.m','file')
          
                       try
                          %Call the "getMETpara_GPT.m" external fxn 
                          [~,T] = getMETpara_GPT(UTCtime,ReceiverPos);
                       catch    
                            %Call the "getTEMP_GPT.m" internal fxn 
                            T = getTEMP_GPT(UTCtime,ReceiverPos);
                       end   
         
                    else    
                        %Call the "getTEMP_GPT.m" internal fxn 
                        T = getTEMP_GPT(UTCtime,ReceiverPos);
          
                    end  %if exist('getMETpara_GPT.m','file') 
             
                 end   
         
          else %IF NONE OF GPT,UNB3m,EGNOS & MOPS MODELS
        
               %GET STORED MET PARAMETERS
               metPARA = getappdata(0,'metPARA');
       
               if ~isempty(metPARA) %IF MET SOURCE IS NOT EMPTY
           
                  T = metPARA(:,1); %TEMPERATURE IS AT THE 1ST COLUMN
               end  
       
          end %//if any([any([strncmpi(Tropo_Cmodel,'GPT2 (5° x 5°)',14),strncmpi(Tropo_Cmodel,'GPT2w (1° x 1°)',14),strncmpi(Tropo_Cmodel,'GPT2w (5° x 5°)',14),...
              %               strncmpi(Tropo_Cmodel,'GPT3 (1° x 1°)',14),strncmpi(Tropo_Cmodel,'GPT3 (5° x 5°)',14)]),any(Tropo_Cmodel_v == [12,13,14,15,16])]) 

        end  %//if option_Cmodels == 1

        %(2.0)**************SEPARATE TROPO MODELs

        if option_Smodels == 1 %if separate model option is selected [Separate  Models]
    
           %GET DRY & WET MODELs
           dryModel = getappdata(0,'dryModel');%get selected Hydrostatic model(string)from pop_ITropo_hydrostatic
           wetModel = getappdata(0,'wetModel');%get selected Wet model(string)from pop_ITropo_wet        
    
           dryModel_v = getappdata(0,'dryModel_v');%get selected Hydrostatic model(Numeric)from pop_ITropo_hydrostatic
           wetModel_v = getappdata(0,'wetModel_v');%get selected Wet model(Numeric)from pop_ITropo_wet
   
           %CHECK IF SELECTED MODEL IS ONE OF THE GLOBAL PRESSURE & TEMPERATURE(GPT) MODELS
           if any([any([any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                             strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),any(dryModel_v==[10,11,12,13,14])]),...        
                   any([any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                             strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)]),any(wetModel_v==[15,16,17,18,19])])])  
            
              T = getappdata(0,'T_GPT');%Retrieve stored T values
       
                %IF SELECTED MODEL IS EITHER UNB3m,EGNOS OR MOPS MODEL     
           elseif any([any([strncmpi(dryModel,'UNB3m',5),strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4)]),...  
                       any([strncmpi(wetModel,'UNB3m',5),strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4)])])
           
                  T = getappdata(0,'T_uem');%Retrieve stored T values
       
                 %IF SELECTED MODEL IS THE VMF ZENITH TROPOSPHERIC DELAY(ZPD) OR GTrop      
           elseif  any([strncmpi(dryModel,'VMF gridded ZHD',15),strncmpi(wetModel,'VMF gridded ZWD',15),...
                        any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),strfind(dryModel,'GTrop')]),... 
                        any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),strfind(wetModel,'GTrop')])])    %if VMF ZHD,VMF ZWD,GTrop
           
                   %COMPUTE SURFACE TEMPERATURE USING GPT MODEL IF TIME & STATION 
                   %POSITION INPUTS ARE NOT EMPTY=============================
                   if all([~isempty(UTCtime),~isempty(ReceiverPos)])
      
                      if exist('getMETpara_GPT.m','file')
          
                         try
                            %Call the "getMETpara_GPT.m" external fxn 
                            [~,T] = getMETpara_GPT(UTCtime,ReceiverPos);
                         catch    
                              %Call the "getTEMP_GPT.m" internal fxn 
                              T = getTEMP_GPT(UTCtime,ReceiverPos);
                         end   
         
                      else   
                          %Call the "getTEMP_GPT.m" internal fxn 
                          T = getTEMP_GPT(UTCtime,ReceiverPos);
          
                      end  %if exist('getMETpara_GPT.m','file') 
             
                   end 
           
           else %IF NONE OF GPT,UNB3m,EGNOS & MOPS MODELS
        
                %GET STORED MET PARAMETERS
                metPARA = getappdata(0,'metPARA');
       
                if ~isempty(metPARA) %IF MET SOURCE IS NOT EMPTY
           
                   T = metPARA(:,1); %TEMPERATURE IS AT THE 1ST COLUMN
                end   
       
           end  
   
        end %//if option_Smodels == 1
    
   end %//if any([~isempty(Ts),all(Ts(:) ~=0),any(find(Ts(:) ~=0)),any(Ts(:) ~=0)])

   %************FINAL CHECK ON TEMPERATURE(T)
   if any([~exist('T','var'),isempty(T)])
    
      %COMPUTE SURFACE TEMPERATURE USING GPT MODEL IF TIME & STATION 
      %POSITION INPUTS ARE NOT EMPTY=============================
      if all([~isempty(UTCtime),~isempty(ReceiverPos)])
      
         if exist('getMETpara_GPT.m','file')
          
            try
               %Call the "getMETpara_GPT.m" external fxn 
               [~,T] = getMETpara_GPT(UTCtime,ReceiverPos);
            catch  
                 %Call the "getTEMP_GPT.m" internal fxn 
                 T = getTEMP_GPT(UTCtime,ReceiverPos);
            end 
         
         else 
             %Call the "getTEMP_GPT.m" internal fxn 
             T = getTEMP_GPT(UTCtime,ReceiverPos);
          
         end %if exist('getMETpara_GPT.m','file')
      
      else 
          if ~exist('T','var')
           
            T = [];
          end 
         
      end %//if all([~isempty(UTCtime),~isempty(ReceiverPos)])
   
   end %//if any([~exist('T','var'),isempty(T)])
   
else %IF NOT ANY OF THE T DEPENDENT MODELS
    
    T = [];%ASSIGN EMPTY([]) MATRIX TO T
    
end %//any([strncmpi(Tm_model,'Bevis et al (1992) model',24),strncmpi(Tm_model,'Bevis et al (1995) model',24),strncmpi(Tm_model,'Mendes et al (2000) model',25),...
    %       strncmpi(Tm_model,'Schueler et al (2001) model',27),strncmpi(Tm_model,'Yao et al (2014) model',22),strncmpi(Tm_model,'TVGG-Tm model [Jiang et al 2018]',32)])
        
 
%B.********COMPUTE PI, THE WATER VAPOUR CONVERSION FACTOR
%Call the "computePI.m" Function
%[PI,Tm] = computePI(Tm_model,Time,ReceiverPos,Temp,GPT_grid,gridRES,Timevar);
[PI,Tm] = computePI(Tm_model,UTCtime,ReceiverPos,T,gridV_TVGG,gridV_GTVR,gridV_GTrop,grid_GPT,grid_res_GPT,Timevar_GPT);
          
%C.********COMPUTE PRECIPITABLE OR INTEGRATED WATER VAPOUR(PWV/IWV)

%(1)Compute Zenith Wet Delay
ZWD = (ZTD - ZHD).*1000; %ZWD IN [ mm ];

%(2)Compute (PWV/IWV)
if any([size(PI,2) > 1,size(ZWD,2) > 1])
    
   pwv = zeros(size(ZWD,2));%Initialize output
   
   for i = 1 : size(PI,2) 
            
       pwv(:,i) = (PI(:,i).*ZWD(:,i)); %[ mm ]
   end
   
else
    
    pwv = PI.*ZWD; %[ mm ]
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF compute_PWV.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%             -------------------------------------------------------
%A.***********SUBROUTINE TO COMPUTE SURFACE TEMPERATURE(T)BY GPT MODEL 
%             --------------------------------------------------------
function [T,P,undu] = getTEMP_GPT(Time,Lat,Lon,h)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%           "getTEMP_GPT" computes surface Meteorological parameters       *
%           such as pressure,temperature,and geoid undulation for a given  *
%           latitude,height and day of year(doy)                           *
%           It is based on on Spherical Harmonics up to degree and order 9 *

%USAGE:                                                                    *
%      [T,P,undu]=getTEMP_GPT(Time,ReceiverPos)                            *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%INPUT:                                                                    *
%1.     Time : Time of observation                                         *
% Time format:[year month day hour minute seconds] OR [year month day]     *                  
%2.      lat : Latitude of Station in degrees OR [D M S]                   *
%3.      lon : Longitude of Station in degrees OR [D M S]                  *
%4.      h   : Station Height(Ellipsoidal) in meters                       *

%OUTPUT:                                                                   *
%1.      T : Surface temperature in degree Celcius(°C)                     *
%2.      P : Surface pressure in millibar(mbar)                            *
%3.    undu: geoid undulation in meters                                    *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================
%REFERENCE: 
% J. Böhm, R. Heinkelmann, H. Schuh, Short Note: A Global Model of Pressure 
% and Temperature for Geodetic Applications, Journal of Geodesy, 
% doi:10.1007/s00190-007-0135-3, 2007.
%==========================================================================

%********CHECK INPUT
if nargin == 4 %IF ALL 4 INPUTS ARE PROVIDED

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


   %****CHECK LONGITUDE INPUT 
   ncol_lon= size(Lon,2); %Get number of longitude entry  
   
   switch ncol_lon
                                
       case  3     %if input is in DMS,Convert to degrees
           lon = dms2degrees(Lon);%Lon in degrees      
       case  2   %if input is in DM,Convert to degrees    
           lon = dm2degrees(Lon);%Lon in degrees
        
       otherwise     
                lon = Lon;
   end 
   
   
elseif nargin == 2 %IF ONLY 2 INPUTS ARE PROVIDED ASSUMING POSITION IS IN [LAT LON h]
       
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
       
end

%****CONVERT LATITUDE & LONGITUDE COORDS IN DEGREES TO RADIAN
lat=deg2rad(lat); %Latitude coord in radian
lon=deg2rad(lon); %Longitude coord in radian


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

%******************************************END OF getMETpara_GPT.m 
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





