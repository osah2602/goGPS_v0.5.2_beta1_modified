function [STD,SHD,SWD,ZTD,ZHD,ZWD] = tropo_error_correction(UTCtime,ReceiverPos,SatPos)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            This function Computes the pseudo-range correction due to     *  
%            tropospheric refraction caused by the Earth's neutral         *
%            atmosphere.It uses the various Zenith Tropopheric Delay and   *
%            Mapping Function Models.
% SYNTAX:                                                                  *
%        [STD,SHD,SWD,ZTD,ZHD,ZWD] = tropo_error_correction(UTCtime,...    *
%                                                      ReceiverPos,SatPos )*
%     OR:                                                                  *
%        [STD,SHD,SWD,ZTD,ZHD,ZWD] = tropo_error_correction(UTCtime,...    *
%                                                 ReceiverPos,satElevation)*
%%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
% INPUT:                                                                   *
%       Generally, the function accepts five(4) sets of inputs:            *
%1.    UTCtime : Receiver reception time in[year,Month,Day,Hour,Min,Sec]   *                                                     
%2.ReceiverPos : Receiver position in either Latitude,Longitude &          *
%                height or ECEF([X Y Z) Coordinates                        *
%3.     SatPos : Satellite position(s) in ECEF([X Y Z]) Coordinates        *
%   OR :                                                                   *
%3.     SatPos : Satellite Elevation in decimal degrees                    *
%4.    METpara : Meteorological parameters in the form [P T e Tm lambda] or*
%                [P T e] or just indicate in text format source of MET     * 
%                PARAMETERs.EXAMPLE OF MET PARAMETER SOURCE IN THIS ROUTINE
%                ARE:'standard','UNB3m','EGNOS','GPT','GPT2','GPT2w','GPT3'*
%                SO in a nutshell, you provide either [P T e Tm lambda] or *
%                [P T e] or any of these:                                  *
%                'standard','UNB3m','EGNOS','GPT','GPT2','GPT2w','GPT3'    *
%EX: tropo_error_correction(UTCtime,ReceiverPos,SatPos,[P T e Tm lambda])  *
%    tropo_error_correction(UTCtime,ReceiverPos,SatPos,'UNB3m')            *
%    tropo_error_correction(UTCtime,ReceiverPos,SatPos,'GPT3'),etc.        *

%Other inputs such as METEOROLOGICAL PARAMETERs, "troposheric and mapping  *
%function model type" are called from the gui_goGPS main figure.           *
%==>{METEOROLOGICAL PARAMETERs will include}:                              *                                       
%MAIN:
%4.    Temperature(T): Atmospheric Temperature in Degree Celsius(C)        *
%5.       Pressure(P): Atmospheric Pressure in millibars(mbar /hPa)        *
%6.   RelHumidity(RH): Relative Humidity in(%) eg: 50                      *
%Others include the ff depending on the model type:                        *
%%7.              Tm : Mean temperature of water vapor in kelvin(k)        *
%8.           lambda : water vapour `lapse rate'(dimensionless)]           * 
%%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
%OUTPUT:                                                                   *
%1.     STD => Slant Total Tropospheric Delay in meters                    *
%2.     SHD => Slant Hydrostaic Tropospheric Delay in meters               *
%3.     SWD => slant Wet Tropospheric Delay  in meters                     *
%4.     ZTD => Zenith Total Tropospheric Delay in meters                   *
%5.     ZHD => Zenith Hydrostaic Tropospheric Delay in meters              *
%6.     ZWD => Zenith Wet Tropospheric Delay  in meters                    *
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
%Edited & Modified BY:
%                    OSAH SAMUEL, MSC GEOMATIC ENGINEERING(PhD STUDENT) ===+
%                    Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com
%                    Phone:+233(0)246137410 / +233(0)509438484
%**************************************************************************
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%******PRELIMINARY
%GET GEODETIC COORDs                             
lat  = ReceiverPos(:,1);%LATITUDE
lon  = ReceiverPos(:,2);%LONGITUDE
h    = ReceiverPos(:,3);%ELLIPSOIDAL HEIGHT 

%(A)*************RETRIEVE ALL STORED DATA FROM gui_GPS.fig
%                -----------------------------------------
%(1.0)GET OPTION FOR TROPO MODELLING
option_Cmodels = getappdata(0,'option_Cmodels');%option button for combine tropo models
option_Smodels = getappdata(0,'option_Smodels');%option button for separate tropo models

%                     ========================================
%2.*******************GET TROPOSPHERIC DELAY CORRECTION MODELs
%                     ========================================

%                         --------------------
%(1.0)********************COMBINE TROPO MODELS
%                         --------------------
if option_Cmodels == 1 %if Combine model option is selected [Dry + Wet Delay Models]
   
   %GET SELECTED MODEL 
   Tropo_Cmodel   = getappdata(0,'Tropo_Cmodel');%get selected item (string)from pop_lTropo_combine
   TROPO_Cmodel   = getappdata(0,'TROPO_Cmodel');%Tropo models for goGPS output file 
   Tropo_Cmodel_v = getappdata(0,'Tropo_Cmodel_v');%get selected item (Numeric)from pop_lTropo_combine
   
   %CHECK IF SELECTED MODEL IS ONE OF THE GPT MODELS
   if any([any([strncmpi(Tropo_Cmodel,'GPT2 (5° x 5°)',14),strncmpi(Tropo_Cmodel,'GPT2w (1° x 1°)',14),strncmpi(Tropo_Cmodel,'GPT2w (5° x 5°)',14),...
                strncmpi(Tropo_Cmodel,'GPT3 (1° x 1°)',14),strncmpi(Tropo_Cmodel,'GPT3 (5° x 5°)',14)]),any(Tropo_Cmodel_v == [12,13,14,15,16])]) 
      
      %GET STORED GPT GRID VALUES & RESOLUTION
      grid_c     = getappdata(0,'grid_c');%VMF GRID FILE
      grid_res_c = getappdata(0,'grid_res_c');%VMF GRID RESOLUTION 
      Timevar_c  = getappdata(0,'Timevar_c');%time variation
      
   %IF SELECTED MODEL IS GTrop MODEL
   elseif any([strncmpi(Tropo_Cmodel,'GTrop [Sun et al 2019]',22),Tropo_Cmodel_v == 18])
          
          %GET STORED GTrop GRID VALUES/COEFFICINETS
          GTropCoeff = getappdata(0,'gridV_GTrop');
          
   elseif any([strncmpi(Tropo_Cmodel,'VMF ZTD',7),Tropo_Cmodel_v == 17])
          
          %*********GET STORED VMF GRID FILES [STRUCT GRID FILES]
          VMFgrids = getappdata(0,'VMFgrids');
          
          %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
          setappdata(0,'VMF_grid_found',getappdata(0,'VMF_grid_found'))    
    
   end
    
 %GET REFERENCE GEOID MODEL(EGM2008 Geopotential Model)
 %-------------------------------------------------------------------------  
%NOTE:
%SAASTAMOINEN,MARINI,Askne & Nordius,UNB3m,EGNOS & MOPS MODELS REQUIRE +VE
%ORTHOMETRIC HEIGHT--------------------------------------------------------
if any([strncmpi(Tropo_Cmodel,'Saastamoinen',12),strncmpi(Tropo_Cmodel,'Saastamoinen(Refined)',21),...
        strncmpi(Tropo_Cmodel,'Askne & Nordius',15),strncmpi(Tropo_Cmodel,'Marini',6),...
        strncmpi(Tropo_Cmodel,'UNB3m',5),strncmpi(Tropo_Cmodel,'EGNOS',5),strncmpi(Tropo_Cmodel,'MOPS',4)])   
    
   geoid = getappdata(0,'geoid');
end
     
end %//if option_Cmodels == 1 

%                          ---------------------                   
%(2.0)*********************SEPARATE TROPO MODELs
%                          --------------------- 

if option_Smodels == 1 %if separate model option is selected [Separate  Models]
    
   %GET DRY & WET MODELs
   dryModel = getappdata(0,'dryModel');%get selected Hydrostatic model(string)from pop_ITropo_hydrostatic
   wetModel = getappdata(0,'wetModel');%get selected Wet model(string)from pop_ITropo_wet
   
   DryModel = getappdata(0,'DryModel');%Hydrostatic model for goGPS output file 
   WetModel = getappdata(0,'WetModel');%Wet model for goGPS output file 
     
   dryModel_v = getappdata(0,'dryModel_v');%get selected Hydrostatic model(Numeric)from pop_ITropo_hydrostatic
   wetModel_v = getappdata(0,'wetModel_v');%get selected Wet model(Numeric)from pop_ITropo_wet
   
   %GET GRIDFILE IF USER SELECTS ANY OF THE HYDROSTATIC & WET GPT MODELs 
   %*****HYDROSTATIC
   if any([any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),any(dryModel_v==[10,11,12,13,14])])
           
      %RETRIEVE SAVED HYDROSTATIC GPT GRID FILE & RESOLUTION 
      grid_h       = getappdata(0,'grid_h');%VMF GRID FILE
      grid_res_h   = getappdata(0,'grid_res_h'); %VMF GRID RESOLUTION
      Timevar_dry  = getappdata(0,'Timevar_dry');%time variation
      
   else %IF NOT ANY OF THE GPT MODELs
        %*****ASSIGN EMPTY([]) MATRIX
        grid_h       = [];%VMF GRID FILE
        grid_res_h   = []; %VMF GRID RESOLUTION
        Timevar_dry  = 1;%set time variation to static    
   end
   
   
   %*********WET   
   if any([any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
           strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)]),any(wetModel_v==[15,16,17,18,19])])  
 
      %RETRIEVE SAVED WET GPT GRID FILE & RESOLUTION 
      grid_w      = getappdata(0,'grid_w');%VMF GRID FILE
      grid_res_w  = getappdata(0,'grid_res_w');%VMF GRID RESOLUTION
      Timevar_wet = getappdata(0,'Timevar_wet');%time variation
      
   else %IF NOT ANY OF THE GPT MODELs
        %*****ASSIGN EMPTY([]) MATRIX
        grid_w       = [];%VMF GRID FILE
        grid_res_w   = []; %VMF GRID RESOLUTION
        Timevar_wet  = 1;%set time variation to static  
      
   end
   
   %IF SELECT MODEL IS GTrop
   if any([any([any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),dryModel_v == 16]),...
                any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),wetModel_v == 21])]),... 
           all([any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),dryModel_v == 16]),...
                any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),wetModel_v == 21])])])
  
      %GET STORED GTrop GRID VALUES/COEFFICINETS
      GTropCoeff = getappdata(0,'gridV_GTrop'); 
      
   else %IF NOT GTrop MODEL
       
        %*****ASSIGN EMPTY([]) MATRIX
        GTropCoeff = [];
   end
   
   %IF SELECT MODEL IS VMF
   if any([any([any([strncmpi(dryModel,'VMF ZHD',7),dryModel_v == 15]),...
                any([strncmpi(wetModel,'VMF ZWD',7),wetModel_v == 20])]),... 
           all([any([strncmpi(dryModel,'VMF ZHD',7),dryModel_v == 15]),...
                any([strncmpi(wetModel,'VMF ZWD',7),wetModel_v == 20])])])
  
      %*********GET STORED VMF GRID FILES [STRUCT GRID FILES]
      VMFgrids = getappdata(0,'VMFgrids');
      
      %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
      setappdata(0,'VMF_grid_found',getappdata(0,'VMF_grid_found')) 
      
   else %IF NOT VMF MODEL
       
        %*****ASSIGN EMPTY([]) MATRIX
        VMFgrids = [];
   end
   
   
   %GET REFERENCE GEOID MODEL(EGM2008 Geopotential Model)
   %-----------------------------------------------------------------------  
   %NOTE"
   %SAASTAMOINEN,DAVIS ET AL,Askne & Nordius,UNB3m,EGNOS & MOPS MODELS 
   %REQUIRE +VE ORTHOMETRIC HEIGHT-----------------------------------------
   if any([strncmpi(dryModel,'Saastamoinen',12),strncmpi(dryModel,'Davis et al)',12),strncmpi(dryModel,'Askne & Nordius',15),...
           strncmpi(dryModel,'UNB3m',5),strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),strncmpi(wetModel,'Saastamoinen',12),...
           strncmpi(wetModel,'Askne & Nordius',15),strncmpi(wetModel,'UNB3m',5),strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4)])  
    
      geoid = getappdata(0,'geoid');
   end 
   
   %                       -----------------
   %(3.0)******************MAPPING FUNCTIONS  
   %                       -----------------
   
   %(3.1)GET OPTIONs FOR USE OF MODELs / OTHER MAPPING FUNCTION
   option_model_MF = getappdata(0,'option_model_MF');%option button for models MF
   option_different_MF = getappdata(0,'option_different_MF');%option button for different MF models

%(3.2)GET SELECTED MAPPING FUNCTION(MF) MODELs
if option_different_MF == 1 %IF User Selects Option for Use of Different MF
   
   %GET SELECTED MF MODELS
   MFh_model = getappdata(0,'MFh_model');%Hydrostatic mapping function Model
   MFw_model = getappdata(0,'MFw_model');%Wet mapping function Model
   
   MFh_Model = getappdata(0,'MFh_Model');%Hydrostatic mapping function Model for goGPS Output file
   MFw_Model = getappdata(0,'MFw_Model');%Wet mapping function Model for goGPS Output file
   
   %*********WHEM MF MODEL IS ANY OF THE VMF MODELs
   %HYDROSTATIC
   if any([strncmpi(MFh_model,'VMF1(1° x 1°)',12) strncmpi(MFh_model,'VMF1(5° x 5°)',12) ...
           strncmpi(MFh_model,'VMF3(1° x 1°)',12) strncmpi(MFh_model,'VMF3(5° x 5°)',12)])
       
      %RETRIEVE SAVED HYDROSTATIC VMF GRID FILES & GRID RESOLUTIONS                                                             
      grid_MFh    = getappdata(0,'grid_MFh');%VMF GRID FILE
      gridRES_MFh = getappdata(0,'gridRES_MFh');%VMF GRID RESOLUTION
      
       %*********GET STORED VMF GRID FILES [STRUCT GRID FILES]
       VMFgrids = getappdata(0,'VMFgrids');
      
   else %IF MF MODEL IS NEITHER OF THE VMF MODELs
        %ASSIGN EMPTY([]) MATRIX
        grid_MFh    = [];%VMF GRID FILE
        gridRES_MFh = [];%VMF GRID RESOLUTION
             
   end
   
   %WET
   if any([strncmpi(MFw_model,'VMF1(1° x 1°)',12) strncmpi(MFw_model,'VMF1(5° x 5°)',12) ...
           strncmpi(MFw_model,'VMF3(1° x 1°)',12) strncmpi(MFw_model,'VMF3(5° x 5°)',12)])
    
      %RETRIEVE SAVED WET VMF GRID FILES & GRID RESOLUTIONS                                                              
      grid_MFw    = getappdata(0,'grid_MFw');%VMF GRID FILE
      gridRES_MFw = getappdata(0,'gridRES_MFw');%VMF GRID RESOLUTION 
      
      %*********GET STORED VMF GRID FILES [STRUCT GRID FILES]
      VMFgrids = getappdata(0,'VMFgrids');
      
   else %IF MF MODEL IS NEITHER OF THE VMF MODELs
        %ASSIGN EMPTY([]) MATRIX
        grid_MFw    = [];%VMF GRID FILE
        gridRES_MFw = [];%VMF GRID RESOLUTION
        
   end
   
   if ~exist('VMFgrids','var')
       
       VMFgrids = [];%ASSIGN EMPTY([]) MATRIX
   end
   
end %//if option_different_MF == 1

end %//if option_Smodels == 1


%                ------------------------------------------------
%(4.0)***********GET SELECTED SOURCE OF METEOROLOGICAL PARAMETERs
%                ------------------------------------------------

metENABLE=getappdata(0,'metENABLE');%GET ENABLE STATE(on or off)

%CHECK IF MET PARAMETERS popup menu IS ENABLE 'on'
if strcmpi(metENABLE,'on')
    
   METpara=getappdata(0,'METpara');
   metVAL=getappdata(0,'metVAL');%SELECTED MET PARAMETER OPTION(Numeric) 

%CHECK IF USER SELECTS ANY OF THE GPT MODELs AS SOURCE
if any([any(strcmpi(METpara,{'GPT2','GPT2w','GPT3'})),any(metVAL==[3,4,5])])
   grid_Met     = getappdata(0,'grid_Met');%GET SAVE GRIDFILE FOR MET PARAMETER
   grid_res_Met = getappdata(0,'grid_res_Met');%GET SAVE GRID RESOLUTION FOR MET PARAMETER
   Timevar_Met  = getappdata(0,'Timevar_met');%time variation
   
   switch Timevar_Met
       case 0
              Timevar_met = 'time variation (annual and semiannual terms)';
                         
       case 1
             Timevar_met = 'no time variation but static quantities';
   end 
      
end

%CHECK IF USER SELECTS TO IMPORT / ENTER MET PARAMETERS
if any(metVAL==[7,8])
   %GET STORED MET DATA FORMAT FROM [importMETPARA_csv  metpara_table sub-routines. REFER TO goGUIcallback_ATMOSmodelling]
   METformat = getappdata(0,'metPARAformat');
   
   %GET STORED MET SOURCE
   METsource = getappdata(0,'METsource'); 
end
%*****GET / COMPUTE METEOROLOGICAL PARAMETER AS SPECIFY BY USER
if all([all([~isempty(METpara),~ischar(METpara)]),any(metVAL==[7,8])]) 
   P0 = mean(METpara(:,1));%Pressure; mean in case # row > 1
   T0 = mean(METpara(:,2));%Temperature; mean in case  # row > 1
   
   %CHECK IF THE LAST COLUMN IN MET DATA IS RH AND COMPUTE e 
   if any([strcmpi(METformat,'P,T,RH') strcmpi(METformat,'T,P,RH')] ) 
   
      RH0 = mean(METpara(:,3));%Relative humidity; mean in case  # row > 1
      
     %CHECK IF USER ENTERED METEOROLOGICAL PARAMETERS IS FROM A WEATHER
     %STATION AND INTERPOLATE TO GNSS SITE
     if ~isempty(METsource) & strcmpi(METsource,'WEATHER STATION')
      
        %***COMPUTE ATMOSPHERIC PARAMETERS FROM SYNOPTIC STATION - Berg, 1948
        
        %SAASTAMOINEN,DAVIS ET AL, MARINI & Askne & Nordius MODELS REQUIRE 
        %+VE ORTHOMETRIC HEIGHT--------------------------------------------
        if option_Cmodels == 1 %if Combine model option is selected
            
           if any([strncmpi(Tropo_Cmodel,'Saastamoinen',12),strncmpi(Tropo_Cmodel,'Saastamoinen(Refined)',21),...
                   strncmpi(Tropo_Cmodel,'Askne & Nordius',15),strncmpi(Tropo_Cmodel,'Marini',6)])   
               
              %COMPUTE UNDULATION & ORTHOMETRIC HEIGHT 
              [undu,horth] = GETundulation(geoid,lat,lon,h);  
          
              %***COMPUTE ATMOSPHERIC PARAMETERS FROM SYNOPTIC STATION
              P  = P0 .* (1-0.0000226.*horth).^5.225;%pressure [mbar] @ GNSS altitude(h)
              T  = (T0 + 273.15) - 0.0065.*horth ;%temperature [K] @ GNSS altitude(h)
              RH = RH0 .* exp(-0.0006396.*horth);%Relative Humidity(%)@ GNSS altitude(h) 
      
              %CONVERT TEMPERATE BACK TO °C 
              T = T - 273.15; %Temperture in kelvin(K) converted to °C  
              
           else
               %***COMPUTE ATMOSPHERIC PARAMETERS FROM SYNOPTIC STATION
               P  = P0 .* (1-0.0000226.*h).^5.225;%pressure [mbar] @ GNSS altitude(h)
               T  = (T0 + 273.15) - 0.0065.*h ;%temperature [K] @ GNSS altitude(h)
               RH = RH0 .* exp(-0.0006396.*h);%Relative Humidity(%)@ GNSS altitude(h) 
               
               %CONVERT TEMPERATE BACK TO °C 
               T = T - 273.15; %Temperture in kelvin(K) converted to °C 
           end
           
        elseif option_Smodels == 1 %if separate model option is selected
                
               if any([strncmpi(dryModel,'Saastamoinen',12),strncmpi(dryModel,'Davis et al)',12),strncmpi(dryModel,'Askne & Nordius',15),...
                       strncmpi(wetModel,'Saastamoinen',12),strncmpi(wetModel,'Askne & Nordius',15)])
                   
                  %COMPUTE UNDULATION & ORTHOMETRIC HEIGHT           
                  [undu,horth] = GETundulation(geoid,lat,lon,h);  
          
                  %***COMPUTE ATMOSPHERIC PARAMETERS FROM SYNOPTIC STATION
                  P  = P0 .* (1-0.0000226.*horth).^5.225;%pressure [mbar] @ GNSS altitude(h)
                  T  = (T0 + 273.15) - 0.0065.*horth ;%temperature [K] @ GNSS altitude(h)
                  RH = RH0 .* exp(-0.0006396.*horth);%Relative Humidity(%)@ GNSS altitude(h) 
      
                  %CONVERT TEMPERATE BACK TO °C 
                  T = T - 273.15; %Temperture in kelvin(K) converted to °C 
                
               else
                    %***COMPUTE ATMOSPHERIC PARAMETERS FROM SYNOPTIC STATION
                    P  = P0 .* (1-0.0000226.*h).^5.225;%pressure [mbar] @ GNSS altitude(h)
                    T  = (T0 + 273.15) - 0.0065.*h ;%temperature [K] @ GNSS altitude(h)
                    RH = RH0 .* exp(-0.0006396.*h);%Relative Humidity(%)@ GNSS altitude(h)  
                  
                    %CONVERT TEMPERATE BACK TO °C 
                    T = T - 273.15; %Temperture in kelvin(K) converted to °C 
               end     
        end
     else
         %ASSIGNMENT
         T  = T0; %Temperture in [°C ]
         P  = P0; %pressure in [mbar]
         RH = RH0; %Relative Humidity in [%]
     end
     
      %COMPUTE PARTIAL PRESSURE OF WATER VAPOR(e)
      e = 0.01 .* RH .* exp(-37.2465 + 0.213166.*(T + 273.15) - 0.000256908.*(T + 273.15).^2);
     
   else
       %CHECK IF USER ENTERED METEOROLOGICAL PARAMETERS IS FROM A WEATHER
       %STATION AND INTERPOLATE TO GNSS SITE
       if ~isempty(METsource) & strcmpi(METsource,'WEATHER STATION')
           
          %***COMPUTE ATMOSPHERIC PARAMETERS FROM SYNOPTIC STATION - Berg, 1948
          
          %SAASTAMOINEN,DAVIS ET AL, MARINI & Askne & Nordius MODELS REQUIRE 
          %+VE ORTHOMETRIC HEIGHT--------------------------------------------
          if option_Cmodels == 1 %if Combine model option is selected
            
             if any([strncmpi(Tropo_Cmodel,'Saastamoinen',12),strncmpi(Tropo_Cmodel,'Saastamoinen(Refined)',21),...
                     strncmpi(Tropo_Cmodel,'Askne & Nordius',15),strncmpi(Tropo_Cmodel,'Marini',6)])   
              
                %COMPUTE UNDULATION & ORTHOMETRIC HEIGHT 
                [undu,horth] = GETundulation(geoid,lat,lon,h);  
          
                 %***COMPUTE ATMOSPHERIC PARAMETERS FROM SYNOPTIC STATION
                 P  = P0 .* (1-0.0000226.*horth).^5.225;%pressure [mbar] @ GNSS altitude(h)
                 T  = (T0 + 273.15) - 0.0065.*horth ;%temperature [K] @ GNSS altitude(h)
              
                 %CONVERT TEMPERATE BACK TO °C 
                 T = T - 273.15; %Temperture in kelvin(K) converted to °C  
             else 
                 %***COMPUTE ATMOSPHERIC PARAMETERS FROM SYNOPTIC STATION
                 P  = P0 .* (1-0.0000226.*h).^5.225;%pressure [mbar] @ GNSS altitude(h)
                 T  = (T0 + 273.15) - 0.0065.*h ;%temperature [K] @ GNSS altitude(h)
                 
                 %CONVERT TEMPERATE BACK TO °C 
                 T = T - 273.15; %Temperture in kelvin(K) converted to °C 
             end 
           
          elseif option_Smodels == 1 %if separate model option is selected
                
                 if any([strncmpi(dryModel,'Saastamoinen',12),strncmpi(dryModel,'Davis et al)',12),strncmpi(dryModel,'Askne & Nordius',15),...
                         strncmpi(wetModel,'Saastamoinen',12),strncmpi(wetModel,'Askne & Nordius',15)])
                   
                     %COMPUTE UNDULATION & ORTHOMETRIC HEIGHT           
                     [undu,horth] = GETundulation(geoid,lat,lon,h);  
          
                     %***COMPUTE ATMOSPHERIC PARAMETERS FROM SYNOPTIC STATION
                     P  = P0 .* (1-0.0000226.*horth).^5.225;%pressure [mbar] @ GNSS altitude(h)
                     T  = (T0 + 273.15) - 0.0065.*horth ;%temperature [K] @ GNSS altitude(h)
                  
                     %CONVERT TEMPERATE BACK TO °C 
                     T = T - 273.15; %Temperture in kelvin(K) converted to °C 
                
                 else 
                      %***COMPUTE ATMOSPHERIC PARAMETERS FROM SYNOPTIC STATION
                      P  = P0 .* (1-0.0000226.*h).^5.225;%pressure [mbar] @ GNSS altitude(h)
                      T  = (T0 + 273.15) - 0.0065.*h ;%temperature [K] @ GNSS altitude(h)
                  
                      %CONVERT TEMPERATE BACK TO °C 
                      T = T - 273.15; %Temperture in kelvin(K) converted to °C 
                 end    
          end 
          
       else   
            %ASSIGNMENT
            T  = T0; %Temperture in [°C ]
            P  = P0; %pressure in [mbar]
         
       end  
       
       e = mean(METpara(:,3));%Water Vapour Pressure; mean in case # row > 1
       
   end  
       
   metPARA = [T P e];

   metPARA_type='User Input Meteorological Parameters';

else  
     if strcmpi(METpara,'Standard')%if user opt for use of standard met para,
                                    %assignment empty matrix([]).The various
                                    %subroutines will compute the Met parameters 
          
        %Use Standard atmosphere @ MSL(h=0)- Berg, 1948 (Bernese)
        try
            %pressure [mbar]
            Ps = goGNSS.STD_PRES;
            %temperature [K]
            Ts = goGNSS.STD_TEMP;
            %humidity [%]
            RHs = goGNSS.STD_HUMI;
        catch
             Ps = 1013.25;%pressure [mbar] @ SEALEVEL
             Ts = 291.15;%temperature [K] @ SEALEVEL
            RHs = 50.0;%Relative Humidity(%)@ SEALEVEL
        end
    
       %SAASTAMOINEN,DAVIS ET AL, MARINI & Askne & Nordius MODELS REQUIRE 
       %+VE ORTHOMETRIC HEIGHT--------------------------------------------
        if option_Cmodels == 1 %if Combine model option is selected
             
           if any([strncmpi(Tropo_Cmodel,'Saastamoinen',12),strncmpi(Tropo_Cmodel,'Saastamoinen(Refined)',21),...
                   strncmpi(Tropo_Cmodel,'Askne & Nordius',15),strncmpi(Tropo_Cmodel,'Marini',6)])   
              
               %COMPUTE UNDULATION & ORTHOMETRIC HEIGHT 
               [undu,horth] = GETundulation(geoid,lat,lon,h);  
          
              %***COMPUTE ATMOSPHERIC PARAMETERS FROM STANDARD PARAMETERS
              P  = Ps .* (1-0.0000226.*horth).^5.225;%pressure [mbar] @ GNSS altitude(h)
              T  = Ts - 0.0065.*horth ;%temperature [K] @ GNSS altitude(h)
              RH = RHs .* exp(-0.0006396.*horth);%Relative Humidity(%)@ GNSS altitude(h) 
      
              %CONVERT TEMPERATE BACK TO °C 
              T = T - 273.15; %Temperture in kelvin(K) converted to °C  
              
           else
               %***COMPUTE ATMOSPHERIC PARAMETERS FROM STANDARD PARAMETERS
               P  = Ps .* (1-0.0000226.*h).^5.225;%pressure [mbar] @ GNSS altitude(h)
               T  = Ts - 0.0065.*h ;%temperature [K] @ GNSS altitude(h)
               RH = RHs .* exp(-0.0006396.*h);%Relative Humidity(%)@ GNSS altitude(h) 
               
               %CONVERT TEMPERATE BACK TO °C 
               T = T - 273.15; %Temperture in kelvin(K) converted to °C 
           end 
           
        elseif option_Smodels == 1 %if separate model option is selected
                
               if any([strncmpi(dryModel,'Saastamoinen',12),strncmpi(dryModel,'Davis et al)',12),strncmpi(dryModel,'Askne & Nordius',15),...
                       strncmpi(wetModel,'Saastamoinen',12),strncmpi(wetModel,'Askne & Nordius',15)])
                   
                  %COMPUTE UNDULATION & ORTHOMETRIC HEIGHT           
                  [undu,horth] = GETundulation(geoid,lat,lon,h);  
          
                  %***COMPUTE ATMOSPHERIC PARAMETERS FROM STANDARD PARAMETERS
                  P  = Ps .* (1-0.0000226.*horth).^5.225;%pressure [mbar] @ GNSS altitude(h)
                  T  = Ts - 0.0065.*horth ;%temperature [K] @ GNSS altitude(h)
                  RH = RHs .* exp(-0.0006396.*horth);%Relative Humidity(%)@ GNSS altitude(h) 
      
                  %CONVERT TEMPERATE BACK TO °C 
                  T = T - 273.15; %Temperture in kelvin(K) converted to °C 
                
               else
                    %***COMPUTE ATMOSPHERIC PARAMETERS FROM STANDARD PARAMETERS
                    P  = Ps .* (1-0.0000226.*h).^5.225;%pressure [mbar] @ GNSS altitude(h)
                    T  = Ts - 0.0065.*h ;%temperature [K] @ GNSS altitude(h)
                    RH = RHs .* exp(-0.0006396.*h);%Relative Humidity(%)@ GNSS altitude(h)  
                  
                    %CONVERT TEMPERATE BACK TO °C 
                    T = T - 273.15; %Temperture in kelvin(K) converted to °C 
               end 
               
        end
  
       %COMPUTE PARTIAL PRESSURE OF WATER VAPOR(e)
       e = 0.01 .* RH .* exp(-37.2465 + 0.213166.*(T + 273.15) - 0.000256908.*(T + 273.15).^2);
                                        
       metPARA = [T P e];
  
       metPARA_type='Standard Meteorological Parameters';
     
     elseif strcmpi(METpara,'GPT')%if user opt for use of GPT met para,
         
            %RETRIEVE SAVED SUPPLEMENTARY GRID FILE & RESOLUTION
            grid_met_ss = getappdata(0,'grid_met_s'); %GRIDFILE 
            grid_res_met_ss = getappdata(0,'grid_res_met_s');% GRID RESOLUTION 
            grid_MET_model = getappdata(0,'grid_METmodel');%GRID MODEL
         
            [P,T]=getMETpara_GPT(UTCtime,ReceiverPos);
            
            %GET SUPPLEMENTARY MET PARAMETERS
            if strcmpi(grid_MET_model,'GPT2w') %if model is 'GPT2w'
               [~,~,e,Tm,lambda]=getMETpara_GPT2w(UTCtime,ReceiverPos,grid_met_ss,grid_res_met_ss); 
        
            elseif strcmpi(grid_MET_model,'GPT3') %if model is 'GPT3'
                   [~,~,e,Tm,lambda]=getMETpara_GPT3(UTCtime,ReceiverPos,grid_met_ss,grid_res_met_ss); 
            end 
     
            %NOTE:GPT T OUTPUT IS IN degree Celsius[°C]
            metPARA=[T P e Tm lambda];
            metPARA_type='GPT Model Generated Meteorological Parameters based on Spherical Harmonics';
            
     elseif  strcmpi(METpara,'GPT2')%if user opt for use of GPT2 met para,
              [P,T,e]=getMETpara_GPT2(UTCtime,ReceiverPos,grid_Met,Timevar_Met);
              
              %NOTE:GPT2 T OUTPUT IS IN degree Celsius[°C]
              metPARA=[T P e];
                 
              metPARA_type=strjoin({'GPT2 Model Generated Meteorological Parameters','based on a',...
                                    strcat(num2str(5),'°'),'x',strcat(num2str(5),'°'),'grid','with',Timevar_met});
                                
     elseif  strcmpi(METpara,'GPT2w')%if user opt for use of GPT2w met para,
             [P,T,e,Tm,lambda]=getMETpara_GPT2w(UTCtime,ReceiverPos,grid_Met,grid_res_Met,Timevar_Met); %#ok<*ASGLU>
              
              %NOTE:GPT2w T OUTPUT IS IN degree Celsius[°C]
              metPARA=[T P e Tm lambda];

              metPARA_type = strjoin({'GPT2w Model Generated Meteorological Parameters','based on a',...
                                      strcat(num2str(grid_res_Met),'°'),'x',strcat(num2str(grid_res_Met),'°'),'grid','with',Timevar_met});                 
     elseif  strcmpi(METpara,'GPT3')%if user opt for use of GPT2w met para,
             [P,T,e,Tm,lambda]=getMETpara_GPT3(UTCtime,ReceiverPos,grid_Met,grid_res_Met,Timevar_Met); %#ok<*ASGLU>
              
              %NOTE:GPT3 T OUTPUT IS IN degree Celsius[°C]
              metPARA=[T P e Tm lambda];

              metPARA_type=strjoin({'GPT3 Model Generated Meteorological Parameters','based on a',...
                                    strcat(num2str(grid_res_Met),'°'),'x',strcat(num2str(grid_res_Met),'°'),'grid','with',Timevar_met});  
                                
       elseif strcmpi(METpara,'UNB3m')%if user opt for use of UNB3m met para,
              [T, P, e,Tm,~,lambda]=getMETpara_UNB3m(UTCtime,ReceiverPos);
              
              %CONVERT T in Kelvin TO degree Celsius[°C]
              %NOTE:UNB3m Model T OUTPUT IS IN degree Kelvin
              T = T - 273.15;
              metPARA=[T P e Tm lambda];
              metPARA_type='UNB3m Model Generated Meteorological Parameters';
              
              
     elseif strcmpi(METpara,'MET file')%if user opt for use of MET file
                                       %MET file contains Met parameters
                                       %observed/measured at GNSS site
                                       
            %**************EXTACT STATION MET PARAMETERS
            %IF USER OBSERVED/MEASURED TEMPERATURE AT GNSS SITE WITH MET DEVICE  
            gs = Go_State.getInstance();
            state = gs.getCurrentSettings();

            if exist('ini_settings_file', 'var')
               state.importIniFile(ini_settings_file);
            end 
                
             md = Meteo_Data(state.getMetFile());%METEO DATA
             
           try  
              if (md.isValid())%IF FILE IS VALID,EXTRACT METEO PARAMETERS
                
                P  = md.getPressure(GPS_Time(datenum(UTCtime(:,:))));%get PRESSURE
                T  = md.getTemperature(GPS_Time(datenum(UTCtime(:,:))));%get TEMPERATURE
                RH = md.getHumidity(GPS_Time(datenum(UTCtime(:,:))));%get Relative humidity 
               
                %COMPUTE PARTIAL PRESSURE OF WATER VAPOUR(e)
                e  = 0.01 .* RH .* exp(-37.2465 + 0.213166.*(T + 273.15) - 0.000256908.*(T + 273.15).^2);
                
                metPARA=[T P e ];
                metPARA_type='Meteorological Parameters observed at GNSS site';
               
              else %IF FILE IS INVALID, %USE GPT2W 1° x 1° MET DATA AS ALTERNATIVE  
                
                  %RETRIEVE SAVED GRID VALUEs & RESOLUTION 
                  grid_METfile     = getappdata(0,'grid_METfile');
                  grid_res_METfile = getappdata(0,'grid_res_METfile');
                
                  Timevar_METfile = 0;%Use default time variation of zero (0)
                
                  %Call the "getMETpara_GPT2w.m" fxn
                  [P,T,e,Tm,lambda]=getMETpara_GPT2w(UTCtime,ReceiverPos,grid_METfile,grid_res_METfile,Timevar_METfile); %#ok<*ASGLU>
              
                  %NOTE:GPT2w T OUTPUT IS IN degree Celsius[°C]
                  metPARA=[T P e Tm lambda];

                  metPARA_type = strjoin({'GPT2w Model Generated Meteorological Parameters','based on a',...
                    
                  strcat(num2str(grid_res_METfile),'°'),'x',strcat(num2str(grid_res_METfile),'°'),'grid','with',Timevar_METfile});                 
                
              end 
            
           catch %SHOULD THERE BE ANY ERROR 
               
                 %RETRIEVE SAVED GRID VALUEs & RESOLUTION 
                 grid_METfile     = getappdata(0,'grid_METfile');
                 grid_res_METfile = getappdata(0,'grid_res_METfile');
                
                 Timevar_METfile = 0;%Use default time variation of zero (0)
                
                 %Call the "getMETpara_GPT2w.m" fxn
                 [P,T,e,Tm,lambda]=getMETpara_GPT2w(UTCtime,ReceiverPos,grid_METfile,grid_res_METfile,Timevar_METfile); %#ok<*ASGLU>
              
                 %NOTE:GPT2w T OUTPUT IS IN degree Celsius[°C]
                 metPARA=[T P e Tm lambda];

                 metPARA_type = strjoin({'GPT2w Model Generated Meteorological Parameters','based on a',...
                    
                 strcat(num2str(grid_res_METfile),'°'),'x',strcat(num2str(grid_res_METfile),'°'),'grid','with',Timevar_METfile});                 
                
           end %//try  
                                       

     end %if strcmpi(METpara,'Standard')
        
end %//if all([all([~isempty(METpara),~ischar(METpara)]),any(metVAL==[7,8])])

%(4.1)****************COMPUTE SUPPLEMENTARY MET PARAMETERS
%NOTE:---------------------------------------------------------------------
%GPT,GPT2 MODELS & USER INPUT MET PARAMETERS DOESN'T PROVIDE ALL THE NEEDED
%METEOROLOGICAL PARAMETERS(e.g.:Tm,lambda).SOME OF THE TROPO
%MODELS(e.g.:'Askne & Nordius' WET MODEL) REQUIRE Mean temperature of water 
%vapor(Tm) & water vapour `lapse rate' or decrease factor(lambda). WE NEED
%TO SUPPLEMENT THE MET PARAMETERS ANYTIME ANY OF THESE MODELS IS SELECTED
%--------------------------------------------------------------------------
if size(metPARA,2) == 3 %if # of Columns in metPARA is 3
    
   %CHECK IF Askne & Nordius TROPO MODEL IS SELECTED
   if option_Cmodels == 1 %if Combine model option is selected
       
      if strncmpi(Tropo_Cmodel,'Askne & Nordius',15)
         
        %RETRIEVE SAVED SUPPLEMENTARY GRID FILE & RESOLUTION
        grid_met_s = getappdata(0,'grid_met_s'); %GRIDFILE 
        grid_res_met_s = getappdata(0,'grid_res_met_s');% GRID RESOLUTION 
        grid_METmodel = getappdata(0,'grid_METmodel');%GRID MODEL
        
        if strcmpi(grid_METmodel,'GPT2w') %if model is 'GPT2w'
           [~,~,~,Tm,lambda]=getMETpara_GPT2w(UTCtime,ReceiverPos,grid_met_s,grid_res_met_s); 
        
        elseif strcmpi(grid_METmodel,'GPT3') %if model is 'GPT3'
               [~,~,~,Tm,lambda]=getMETpara_GPT3(UTCtime,ReceiverPos,grid_met_s,grid_res_met_s); 
        end
     
        metPARA = [metPARA Tm lambda];
        
      end   
                   
   elseif  option_Smodels == 1 %if separate model option is selected
             
           if strncmpi(wetModel,'Askne & Nordius',15)
         
                %RETRIEVE SAVED SUPPLEMENTARY GRID FILE & RESOLUTION
                grid_met_s = getappdata(0,'grid_met_s'); %GRIDFILE 
                grid_res_met_s = getappdata(0,'grid_res_met_s');% GRID RESOLUTION 
                grid_METmodel = getappdata(0,'grid_METmodel');%GRID MODEL
        
               if strcmpi(grid_METmodel,'GPT2w') %if model is 'GPT2w'
                  [~,~,~,Tm,lambda]=getMETpara_GPT2w(UTCtime,ReceiverPos,grid_met_s,grid_res_met_s); 
        
               elseif strcmpi(grid_METmodel,'GPT3') %if model is 'GPT3'
                      [~,~,~,Tm,lambda]=getMETpara_GPT3(UTCtime,ReceiverPos,grid_met_s,grid_res_met_s); 
               end 
     
               metPARA = [metPARA Tm lambda];
               
           end         
   end 
   
end %//if size(metPARA,2) == 3

else %IF MET PARAMETERS popup menu IS ENABLE 'off' 
     %ASSIGN EMPTY MATRIX([])
     metPARA = [];
     metPARA_type = [];
     
end %//if strcmpi(metENABLE,'on')
  
%                   ------------------------------------
%(B).***************COMPUTE TROPOSHERIC ERROR CORRECTIONS
%                   ------------------------------------

%              ----------------------------------------------
%(B.1)*********WHEN USER OPTs FOR COMBINE TROPOSPHERIC MODELS
%              ----------------------------------------------
if option_Cmodels==1 & option_Smodels==0 %if user selects option for combine tropo model
   
   metPARA_typeC = metPARA_type;%MET STRING FOR COMBINE TROPO MODELS
   %--------------------------------------------------------------------------------
   %NOTE: 
   %     "metPARA_type" IS A STRING INDICATING THE SOURCE OF METEOROLOGICAL PARAMETER
   %      GENERATED FROM THE SOURCE OF METEOROLOGICAL PARAMETERS POP UP MENU IN goGPS
   %      HOWEVER, SOME TROPO MODELS GENERATE THEIR OWN MET PARAMETERS EXAMPLES INCLUDE:
   %      UNB3m,ENGNOS,MOPS,AND THE GPT2,GPT2w & GPT3 MODELS.SHOULD ANY OF
   %      THESE MODELS BE SELECTED, THE "metPARA_type" STRING IS BEING
   %      MODIFIED.
   %----------------------------------------------------------------------------------
   if any([strncmpi(Tropo_Cmodel,'GPT2 (5° x 5°)',14),strncmpi(Tropo_Cmodel,'GPT2w (1° x 1°)',14),...
           strncmpi(Tropo_Cmodel,'GPT2w (5° x 5°)',14),strncmpi(Tropo_Cmodel,'GPT3 (1° x 1°)',14),...
           strncmpi(Tropo_Cmodel,'GPT3 (5° x 5°)',14)])
       
      switch Timevar_c
          case 0
                Timevar_C = 'time variation (annual and semiannual terms)';
                         
          case 1
                Timevar_C = 'no time variation but static quantities';
      end 
      
   end 
   
   %****************NOW,WHICH MODEL IS BEEN SELECTED ?? 
   %                -----------------------------------
   if all([strncmpi(Tropo_Cmodel,'Saastamoinen',12),Tropo_Cmodel_v == 2]) %if Saastamoinen
      [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_SAAS(ReceiverPos,SatPos,metPARA);
       
      MF_Model = 'Saastamoinen Mapping Function'; %Mapping Function model for goGPS Output
      
   elseif all([strncmpi(Tropo_Cmodel,'Saastamoinen (Refined)',22),Tropo_Cmodel_v == 3])%if Saastamoinen (Refined)
      
          [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_SAAM(ReceiverPos,SatPos,metPARA);
          
          MF_Model = 'Saastamoinen Mapping Function'; %Mapping Function model for goGPS Output
          
   elseif strncmpi(Tropo_Cmodel,'Hopfield',8)%if Hopfield
       
          [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_Hopfield(ReceiverPos,SatPos,metPARA);
            
          MF_Model = 'Hopfield Mapping Function'; %Mapping Function model for goGPS Output
          
   elseif strncmpi(Tropo_Cmodel,'Modified Hopfield(Goads & Goodman)',34) %if Modified Hopfield(Goads & Goodman)
          [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_GG(ReceiverPos,SatPos,metPARA);
           
          MF_Model = 'Goads & Goodman Mapping Function'; %Mapping Function model for goGPS Output
          
   elseif strncmpi(Tropo_Cmodel,'Marini',6) %Marini 
          [STD,ZTD,ZHD,ZWD] = TropModel_Marini(ReceiverPos,SatPos,metPARA); 
          
          SHD = zeros(size(STD));%Create zero matrix
          SWD = deal(SHD); %Create copy of SHD
          
          MF_Model = 'Marini Mapping Function'; %Mapping Function model for goGPS Output
          
   elseif strncmpi(Tropo_Cmodel,'Askne & Nordius',15)%if Askne & Nordius
          [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_Askne_Nordius(ReceiverPos,SatPos,metPARA); 
          
           MF_Model = 'Black & Eisner Mapping Function'; %Mapping Function model for goGPS Output 
          
   elseif strncmpi(Tropo_Cmodel,'Black',5) %if Black
          [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_Black(ReceiverPos,SatPos,metPARA);
          
          MF_Model = 'Black Mapping Function'; %Mapping Function model for goGPS Output
          
   elseif strncmpi(Tropo_Cmodel,'UNB3m',5)%if UNB3m
          [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_UNB3m(UTCtime,ReceiverPos,SatPos);
          
          %FOR goGPS OUTPUT
          MF_Model     = 'Neill Mapping Function(NMF)';
          metPARA_typeC = 'UNB3m Model Generated Meteorological Parameters';
          
   elseif  strncmpi(Tropo_Cmodel,'EGNOS',5) %if EGNOS
           [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_EGNOS(UTCtime,ReceiverPos,SatPos);
           
           %FOR goGPS OUTPUT
           MF_Model     = 'Black & Eisner Mapping Function';
           metPARA_typeC = 'EGNOS Model Generated Meteorological Parameters';
           
   elseif  strncmpi(Tropo_Cmodel,'MOPS',4) %if MOPS
           [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_MOPS(UTCtime,ReceiverPos,SatPos);
           
           %FOR goGPS OUTPUT
           MF_Model     = 'Black & Eisner Mapping Function';
           metPARA_typeC = 'MOPS Model Generated Meteorological Parameters';
                              
   elseif  strncmpi(Tropo_Cmodel,'GPT2 (5° x 5°)',14) % if GPT2 5° x 5°
           [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_GPT2(UTCtime,ReceiverPos,SatPos,grid_c,Timevar_c);         
           
           %FOR goGPS OUTPUT
           MF_Model     = 'Vienna Mapping Function 1(VMF1)';
           TROPO_Cmodel = 'GPT2 based on a 5° x 5° external grid file';
           metPARA_typeC = strjoin({'GPT2 Model Generated Meteorological Parameters','based on a',...
                                        strcat(num2str(5),'°'),'x',strcat(num2str(5),'°'),'grid','with',Timevar_C});
                                                                                                
   elseif  any([strncmpi(Tropo_Cmodel,'GPT2w (1° x 1°)',14),strncmpi(Tropo_Cmodel,'GPT2w (5° x 5°)',14)]) % if GPT2w 1° or 5°
           [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_GPT2w(UTCtime,ReceiverPos,SatPos,grid_c,grid_res_c,Timevar_c);
           
           %FOR goGPS OUTPUT
           MF_Model     = 'Vienna Mapping Function 1(VMF1)';
           TROPO_Cmodel = strjoin({'GPT2w based on a',strcat(num2str(grid_res_c),'°'),'x',strcat(num2str(grid_res_c),'°'),'external grid file,'});
           metPARA_typeC = strjoin({'GPT2w Model Generated Meteorological Parameters','based on a',...
                                   strcat(num2str(grid_res_c),'°'),'x',strcat(num2str(grid_res_c),'°'),'grid','with',Timevar_C});
                                            
   elseif  any([strncmpi(Tropo_Cmodel,'GPT3 (1° x 1°)',14),strncmpi(Tropo_Cmodel,'GPT3 (5° x 5°)',14)]) % if GPT3 1° or 5°
           [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_GPT3(UTCtime,ReceiverPos,SatPos,grid_c,grid_res_c,Timevar_c);
   
           %FOR goGPS OUTPUT
           MF_Model     = 'Vienna Mapping Function 3(VMF3)';
           TROPO_Cmodel = strjoin({'GPT3 based on a',strcat(num2str(grid_res_c),'°'),'x',strcat(num2str(grid_res_c),'°'),'external grid file,'});
           metPARA_typeC = strjoin({'GPT3 Model Generated Meteorological Parameters','based on a',...
                                   strcat(num2str(grid_res_c),'°'),'x',strcat(num2str(grid_res_c),'°'),'grid','with',Timevar_C});
                                    
   elseif strncmpi(Tropo_Cmodel,'VMF gridded ZTD',15) %if Vinna Mapping Function Zenith tropospheric delays
          [STD,SHD,SWD,ZTD,ZHD,ZWD,~,~,VMF_model] = TropModel_VMF(UTCtime,lat,lon,h,VMFgrids,SatPos);   
          
          %FOR goGPS OUTPUT
          MF_Model     = 'Vienna Mapping Functions(VMF)';
          
          if any([strcmpi(VMF_model,'VMF1'),strfind(VMF_model,'VMF1')])
              
             TROPO_Cmodel = 'Vienna Mapping Function 1(VMF1 ) Gridded Zenith Total Delays';
          elseif any([strcmpi(VMF_model,'VMF3'),strfind(VMF_model,'VMF3')])
              
                  TROPO_Cmodel = 'Vienna Mapping Function 3(VMF3 ) Gridded Zenith Delays';
                 
          else   
              TROPO_Cmodel     = 'Vienna Mapping Function (VMF) Gridded Zenith Total Delays';
          end   
          
          metPARA_typeC = 'Numerical Weather Models (NWMs) Generated Meteorological Parameters';
          
          
   elseif any([strncmpi(Tropo_Cmodel,'GTrop [Sun et al 2019]',22),strfind(Tropo_Cmodel,'GTrop'),Tropo_Cmodel_v == 18])
          
          [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModel_GTrop(UTCtime,lat,lon,h,GTropCoeff,SatPos);
          
          %FOR goGPS OUTPUT
          MF_Model     = 'Neill Mapping Function(NMF)';
          TROPO_Cmodel = ' Global Tropospheric Model(GTrop) based on a 1° x 1° spatial resolution from ECMWF ERA-Interim reanalysis data';
          metPARA_typeC = 'GTrop Model Generated Meteorological Parameters';
          
   elseif  strncmpi(Tropo_Cmodel,'no model',8)%if no model
       
           if ~isempty(SatPos)
              STD = zeros(size(SatPos,1),size(ReceiverPos,1));
           else
               STD = zeros(size(ReceiverPos));
           end
           
           [SHD,SWD,ZTD,ZHD,ZWD] = deal(STD);
           
           %FOR goGPS OUTPUT
           MF_Model     = 'No Mapping Functions';
           metPARA_typeC = 'No Meteorological Parameters Generated';
           
   end %//if strncmpi(Tropo_Cmodel,'Saastamoinen',12)
   
%**************END OF MODELLING TROPO DELAYS USING COMBINE(DRY+WET) MODELS 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
   
%              ----------------------------------------------
%(B.2)*********WHEN USER OPTs FOR SEPARATE TROPOSPHERIC MODELS
%              ----------------------------------------------
   
elseif option_Smodels==1 & option_Cmodels==0 %if user selects option for separate tropo model
    
    %ASSIGNMENT
    metPARA_typeD = metPARA_type;%MET STRING FOR DRY MODELS
    metPARA_typeW = metPARA_type;%MET STRING FOR WET MODELS
   %-----------------------------------------------------------------------
   %NOTE: 
   %     "metPARA_type" IS A STRING INDICATING THE SOURCE OF METEOROLOGICAL PARAMETER
   %      GENERATED FROM THE SOURCE OF METEOROLOGICAL PARAMETERS POP UP MENU IN goGPS
   %      HOWEVER, SOME TROPO MODELS GENERATE THEIR OWN MET PARAMETERS EXAMPLES INCLUDE:
   %      UNB3m,ENGNOS,MOPS,AND THE GPT2,GPT2w & GPT3 MODELS.SHOULD ANY OF
   %      THESE MODELS BE SELECTED, THE "metPARA_type" STRING IS BEING
   %      MODIFIED.
   %-----------------------------------------------------------------------
    if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
            strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)])
              
       switch Timevar_dry
           case 0
                 Timevar_Dry = 'time variation (annual and semiannual terms)';
                         
           case 1
                 Timevar_Dry = 'no time variation but static quantities';
       end               
                    
    end 
    
    if any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
            strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])
              
       switch Timevar_wet
           case 0
                 Timevar_Wet = 'time variation (annual and semiannual terms)';
                         
           case 1
                 Timevar_Wet = 'no time variation but static quantities';
       end                
                    
    end          
                    
      %1.*****HYDROSTATIC/DRY MODELS
      %       ----------------------
       if strncmpi(dryModel,'UNB3m',5)%if UNB3m
          
         %FOR goGPS OUTPUT  
          metPARA_typeD = 'UNB3m Model Generated Meteorological Parameters';
          
       elseif  strncmpi(dryModel,'EGNOS',5) %if EGNOS
               
               %FOR goGPS OUTPUT
               metPARA_typeD = 'EGNOS Model Generated Meteorological Parameters';
           
       elseif  strncmpi(dryModel,'MOPS',4) %if MOPS
               
               %FOR goGPS OUTPUT
               metPARA_typeD = 'MOPS Model Generated Meteorological Parameters';
                              
       elseif  strncmpi(dryModel,'GPT2 (5° x 5°)',14) % if GPT2 5° x 5°
               
               %FOR goGPS OUTPUT
               DryModel     = 'GPT2 based on a 5° x 5° external grid file';
               metPARA_typeD = strjoin({'GPT2 Model Generated Meteorological Parameters','based on a',...
                                       strcat(num2str(5),'°'),'x',strcat(num2str(5),'°'),'grid','with',Timevar_Dry});
                               
       elseif  any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14)]) % if GPT2w 1° or 5°
            
               %FOR goGPS OUTPUT
               DryModel     = strjoin({'GPT2w based on a',strcat(num2str(grid_res_h),'°'),'x',strcat(num2str(grid_res_h),'°'),'external grid file'});
               metPARA_typeD = strjoin({'GPT2w Model Generated Meteorological Parameters','based on a',...
                                       strcat(num2str(grid_res_h),'°'),'x',strcat(num2str(grid_res_h),'°'),'grid',Timevar_Dry});
                                   
       elseif any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]) % if GPT3 1° or 5°
              
              %FOR goGPS OUTPUT
              DryModel     = strjoin({'GPT3 based on a',strcat(num2str(grid_res_h),'°'),'x',strcat(num2str(grid_res_h),'°'),'external grid file'});
              metPARA_typeD = strjoin({'GPT3 Model Generated Meteorological Parameters','based on a',...
                                      strcat(num2str(grid_res_h),'°'),'x',strcat(num2str(grid_res_h),'°'),'grid','with',Timevar_Dry});
                                                                          
       elseif strncmpi(dryModel,'VMF gridded ZHD',15) %if Vinna Mapping Function Zenith tropospheric delays
           
              %FOR goGPS OUTPUT
              DryModel      = 'Vienna Mapping Function (VMF ) Gridded Zenith Hydrostatic Delays';  
              metPARA_typeD = 'Numerical Weather Models (NWMs) Generated Meteorological Parameters';
          
       elseif any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),dryModel_v == 16])
              
              %FOR goGPS OUTPUT
              DryModel     = 'Global Tropospheric Model(GTrop) based on a 1° x 1° spatial resolution from ECMWF ERA-Interim reanalysis data'; 
              metPARA_typeD = 'GTrop Model Generated Meteorological Parameters from ECMWF ERA-Interim reanalysis data';
                    
       end
       
       
      %2.****************WET MODELS
      %                  ----------
       if strncmpi(wetModel,'UNB3m',5)%if UNB3m
          
          %FOR goGPS OUTPUT
          metPARA_typeW = 'UNB3m Model Generated Meteorological Parameters';
          
       elseif  strncmpi(wetModel,'EGNOS',5) %if EGNOS
               
               %FOR goGPS OUTPUT
               metPARA_typeW = 'EGNOS Model Generated Meteorological Parameters';
           
       elseif  strncmpi(wetModel,'MOPS',4) %if MOPS
           
               %FOR goGPS OUTPUT
               metPARA_typeW = 'MOPS Model Generated Meteorological Parameters';
                              
       elseif  strncmpi(wetModel,'GPT2 (5° x 5°)',14) % if GPT2 5° x 5°
            
               %FOR goGPS OUTPUT
               WetModel     = 'GPT2 based on a 5° x 5° external grid file';
               metPARA_typeW = strjoin({'GPT2 Model Generated Meteorological Parameters','based on a',...
                                        strcat(num2str(5),'°'),'x',strcat(num2str(5),'°'),'grid','with',Timevar_Wet});
                               
       elseif  any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14)]) % if GPT2w 1° or 5°
            
               %FOR goGPS OUTPUT
               WetModel     = strjoin({'GPT2w based on a',strcat(num2str(grid_res_w),'°'),'x',strcat(num2str(grid_res_w),'°'),'external grid file'});
               metPARA_typeW = strjoin({'GPT2w Model Generated Meteorological Parameters','based on a',...
                                       strcat(num2str(grid_res_w),'°'),'x',strcat(num2str(grid_res_w),'°'),'grid','with',Timevar_Wet});
                                   
       elseif any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)]) % if GPT3 1° or 5°
           
              %FOR goGPS OUTPUT
              WetModel     = strjoin({'GPT3 based on a',strcat(num2str(grid_res_w),'°'),'x',strcat(num2str(grid_res_w),'°'),'external grid file'});
              metPARA_typeW = strjoin({'GPT3 Model Generated Meteorological Parameters','based on a',...
                                      strcat(num2str(grid_res_w),'°'),'x',strcat(num2str(grid_res_w),'°'),'grid','with',Timevar_Wet});
                                                                          
       elseif strncmpi(wetModel,'VMF gridded ZWD',15) %if Vinna Mapping Function Zenith tropospheric delays
           
              %FOR goGPS OUTPUT
              WetModel      = 'Vienna Mapping Function (VMF ) Gridded Zenith Wet Delays'; 
              metPARA_typeW = 'Numerical Weather Models (NWMs) Generated Meteorological Parameters';
              
       elseif any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),wetModel_v == 21])
              
              %FOR goGPS OUTPUT
              WetModel     = 'Global Tropospheric Model(GTrop) based on a 1° x 1° spatial resolution from ECMWF ERA-Interim reanalysis data'; 
              metPARA_typeW = 'GTrop Model Generated Meteorological Parameters from ECMWF ERA-Interim reanalysis data';       
          
       end 
       
       %IF BOTH DRY & WET TROPO MODELS ARE NEITHER OF THESE,THEN SOURCE OF
       %MET PARAMETERS ARE SAME. [THIS NEEDED FOR goGPS OUTPUT FILE]
       if all([~any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                     strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),...
                     strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),strncmpi(dryModel,'VMF ZHD',7),strncmpi(dryModel,'GTrop [Sun et al 2019]',22)]),...
               ~any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                     strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),...
                     strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4),strncmpi(wetModel,'VMF ZWD',7),strncmpi(wetModel,'GTrop [Sun et al 2019]',22)])])
                  
           metPARA_typeW_D = 1; %FLAG TO INDICATE MET SOURCE ARE THE SAME
           
       else
           metPARA_typeW_D = 0; %FLAG TO INDICATE MET SOURCE ARE not THE SAME
           
       end
                                    
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
       %                   END OF "metPARA_type" MODIFICATION
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       
       %                    ---------------------------
       %********************GET MAPPING FUNCTION OPTION
       %                    ---------------------------
       
       %(A)WHEN USER SELECTS OPTION FOR USE OF DIFFERENT MF MODEL
       if option_different_MF==1 & option_model_MF==0 %if user selects option for use of different MF model
          
          %                 FOR goGPS OUTPUT FILE
          %================================================================
          %***********HYDROSTATIC MF FOR goGPS OUTPUT
          if strncmpi(MFh_model,'VMF1(1° x 1°)',12)
             MFh_Model='Vienna Mapping Function 1(VMF1) based on a 1° x 1° grid';
        
          elseif strncmpi(MFh_model,'VMF1(5° x 5°)',12)
                 MFh_Model='Vienna Mapping Function 1(VMF1) based on a 5° x 5° grid';
            
          elseif strncmpi(MFh_model,'VMF3(1° x 1°)',12)
                 MFh_Model='Vienna Mapping Function 3(VMF3) based on a 1° x 1° grid';
            
          elseif strncmpi(MFh_model,'VMF3(5° x 5°)',12)
                 MFh_Model='Vienna Mapping Function 3(VMF3) based on a 5° x 5° grid';
            
          elseif strncmpi(MFh_model,'NMF',3)
                 MFh_Model='Neill Mapping Function(NMF)';
            
          elseif strncmpi(MFh_model,'GMF',3)
                 MFh_Model='Global Mapping Function(GMF)';
          end  
           
          %*************WET MF FOR goGPS OUTPUT
          if strncmpi(MFw_model,'VMF1(1° x 1°)',12)
             MFw_Model='Vienna Mapping Function 1(VMF1) based on a 1° x 1° grid';
        
          elseif strncmpi(MFw_model,'VMF1(5° x 5°)',12)
                 MFw_Model='Vienna Mapping Function 1(VMF1) based on a 5° x 5° grid';
            
          elseif strncmpi(MFw_model,'VMF3(1° x 1°)',12)
                 MFw_Model='Vienna Mapping Function 3(VMF3) based on a 1° x 1° grid';
            
          elseif strncmpi(MFw_model,'VMF3(5° x 5°)',12)
                 MFw_Model='Vienna Mapping Function 3(VMF3) based on a 5° x 5° grid';
            
          elseif strncmpi(MFw_model,'NMF',3)
                 MFw_Model='Neill Mapping Function(NMF)';
            
          elseif strncmpi(MFw_model,'GMF',3)
                 MFw_Model='Global Mapping Function(GMF)';
          end  
          %================================================================
          
           %**********COMPUTATION OF TROPOSPHERIC CORRECTIONS
           %NOTE:
           %     FOR THE COMPUTATION OF ONLY ZENITH DELAYS, SATELLITE
           %     POSITION(SatPos)IS SET TO EMPTY([]) AS WELL AS MF MODELS
           if isempty(SatPos)
               MFh_model = [];
               MFw_model = [];
           end
          
           [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropModels(UTCtime,ReceiverPos,SatPos,dryModel,wetModel,...
                                                   MFh_model,MFw_model,metPARA,grid_h,grid_res_h,...
                                                  grid_w,grid_res_w,grid_MFh,gridRES_MFh,...
                                                  grid_MFw,gridRES_MFw,Timevar_dry,Timevar_wet,GTropCoeff,VMFgrids);
        
         
       %(B)WHEN USER SELECTS OPTION FOR USE OF MODEL'S  MF                                     
       elseif option_model_MF==1 & option_different_MF==0 %if user selects option for use of tropo models  MF                                    
              
              %                 FOR goGPS OUTPUT FILE
              %============================================================
              %***********HYDROSTATIC DELAYS MF FOR goGPS OUTPUT
              if strncmpi(dryModel,'Saastamoinen',12)
                 MFh_Model = 'Saastamoinen Mapping Function';
                 
              elseif strncmpi(dryModel,'Hopfield',8)  
                     MFh_Model = 'Hopfield Mapping Function'; 
                
              elseif strncmpi(dryModel,'Modified Hopfield(Goads & Goodman)',34)  
                     MFh_Model = 'Goads & Goodman Mapping Function'; 
                     
              elseif any([strncmpi(dryModel,'Askne & Nordius',15),strncmpi(dryModel,'EGNOS',5),...
                          strncmpi(dryModel,'MOPS',4)]) 
                 
                     MFh_Model = 'Black & Eisner Mapping Function'; 
                     
              elseif  strncmpi(dryModel,'Black)',5)  
                      MFh_Model = 'Black Mapping Function';  
                      
              elseif strncmpi(dryModel,'Davis et al)',12)
                     MFh_Model = 'Davis Mapping Function';
                     
              elseif any([strncmpi(dryModel,'UNB3m',5),any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),dryModel_v == 16])])  
                     MFh_Model = 'Neill Mapping Function(NMF)'; 
                         
              elseif any([strncmpi(dryModel,'GPT2 (5° x 5°)',15),strncmpi(dryModel,'GPT2w (1° x 1°)',15),...
                          strncmpi(dryModel,'GPT2w (5° x 5°)',15)])
                   
                     MFh_Model = 'Vienna Mapping Function 1(VMF1)';
               
              elseif any([strncmpi(dryModel,'GPT3 (1° x 1°)',15),strncmpi(dryModel,'GPT3 (5° x 5°)',15)])
            
                     MFh_Model = 'Vienna Mapping Function 3(VMF3)';
         
              elseif strncmpi(dryModel,'VMF gridded ZHD',15)
                     MFh_Model = 'Vienna Mapping Functions(VMF)'; 
                     
                
              end  
                      
              %*************WET DELAYS MF FOR goGPS OUTPUT 
              if strncmpi(wetModel,'Saastamoinen',12)
                 MFw_Model = 'Saastamoinen Mapping Function'; 
                 
              elseif strncmpi(wetModel,'Hopfield',8)  
                     MFw_Model = 'Hopfield Mapping Function'; 
                
              elseif strncmpi(wetModel,'Modified Hopfield(Goads & Goodman)',34)  
                     MFw_Model = 'Goads & Goodman Mapping Function'; 
                     
              elseif any([strncmpi(wetModel,'Askne & Nordius',15),strncmpi(wetModel,'EGNOS',5),...
                          strncmpi(wetModel,'MOPS',4)]) 
                 
                     MFw_Model = 'Black & Eisner Mapping Function'; 
                     
              elseif  strncmpi(wetModel,'Ifadis',6)  
                      MFw_Model = 'Ifadis Mapping Function';  
                      
              elseif strncmpi(wetModel,'Chao)',4)
                     MFw_Model = 'Chao Mapping Function';
                     
             elseif any([strncmpi(wetModel,'Berman 70',9),strncmpi(wetModel,'Berman 74',9),strncmpi(wetModel,'Berman TMOD',11),...
                         strncmpi(wetModel,'Mendes',6),strncmpi(wetModel,'Callahan',8)])
                     
                     MFw_Model = 'Black & Eisner Mapping Function';
                     
                     
              elseif any([strncmpi(wetModel,'UNB3m',5),any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),wetModel_v == 21])]) 
                     MFw_Model = 'Neill Mapping Function(NMF)'; 
                         
              elseif any([strncmpi(wetModel,'GPT2 (5° x 5°)',15),strncmpi(wetModel,'GPT2w (1° x 1°)',15),...
                          strncmpi(wetModel,'GPT2w (5° x 5°)',15)])
                   
                     MFw_Model = 'Vienna Mapping Function 1(VMF1)';
               
              elseif any([strncmpi(wetModel,'GPT3 (1° x 1°)',15),strncmpi(wetModel,'GPT3 (5° x 5°)',15)])
            
                     MFw_Model = 'Vienna Mapping Function 3(VMF3)';
         
              elseif strncmpi(wetModel,'VMF gridded ZWD',15)
                     MFw_Model = 'Vienna Mapping Functions(VMF)';  
                
              end  
              %============================================================

              %**********COMPUTATION OF TROPOSPHERIC CORRECTIONS
              [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropoModels(UTCtime,ReceiverPos,SatPos,dryModel,wetModel,...
                                                       metPARA,grid_h,grid_res_h,grid_w,grid_res_w,...
                                                       Timevar_dry,Timevar_wet,GTropCoeff,VMFgrids);
                                                     
       end      
       
end %//if option_Cmodels==1 & option_Smodels==0

%********SAVE VARIABLES 
%1.combine TROPO MODELS
if exist('TROPO_Cmodel','var')
   setappdata(0,'TROPO_Cmodel',TROPO_Cmodel)
   
   %MAPPING FUNCTIONS
   setappdata(0,'MF_Model',MF_Model)
   
   %STORE SOURCE METEOROLOGICAL PARAMETERS
   setappdata(0,'metPARA_typeC', metPARA_typeC)
end

%2.Separate TROPO MODELS
if any([exist('DryModel','var'),exist('WetModel','var')])
   setappdata(0,'DryModel',DryModel)
   setappdata(0,'WetModel',WetModel)

   %MAPPING FUNCTIONS
   setappdata(0,'MFh_Model',MFh_Model)
   setappdata(0,'MFw_Model',MFw_Model)
   
   %STORE SOURCE METEOROLOGICAL PARAMETERS
   setappdata(0,'metPARA_typeD', metPARA_typeD)
   setappdata(0,'metPARA_typeW', metPARA_typeW)
   setappdata(0,'metPARA_typeW_D', metPARA_typeW_D)
end

%SAVE MET PARAMETERS
setappdata(0,'metPARA', metPARA)

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%            END OF TROPOSPHERIC CORRECTION COMPUTATIONS
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 


%*******SUBROUTINE TO COMPUTE GEOID UNDULATION & ORTHOMETRIC HEIGHT
function [undu,horth] = GETundulation(geoid,lat,lon,h)

%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            This subroutine determines Geoid Undulation and computes      * 
%            Orthometric Height using either goGPS Geoid model(EGM2008) OR *
%            MATLAB EGM96 Geopotential Model                               *                          
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%      geoid = goGPS Geoid model(EGM2008)
%      lat   = Latitude of point[degrees]
%      lon   = Longitude of point[degrees]
%      h     = Ellipsoidal height of point (m)

%OUTPUT: 
%       undu  = Geoid undulation
%       horth = Orthometric Height
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%--------------------------------------------------------------------------

if isempty(geoid)
    
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

end %if isempty(geoid)

%USE goGPS DEFAULT GEOID MODEL
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
%   [N] = grid_bilin_interp(X_approx, Y_approx, grid, ncols,nrows,...      *
%                                              cellsize,Xll, Yll, nodata); *
%--------------------------------------------------------------------------*
% INPUT:
%       X_approx = X coordinate of the interpolation point
%       Y_approx = Y coordinate of the interpolation point
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
        interp_value = nodata;
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
        interp_value = nodata; %#ok<*NASGU>
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

     N(i,1) = interp_value;

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF grid_bilin_interp.m 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF tropo_error_correction.m 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-         	