%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            "TropoModels" is a subroutine that Computes the wet,dry, and   * 
%            Total Tropospheric Delays Using  the various Tropopheric      *  
%            Delay Models                                                  *
%USAGE:                                                                    *
%[STD,SHD,SWD,ZTD,ZHD,ZWD] = TropoModels(Time,ReceiverPos,ElevationAngle,dryModel,wetModel,...
%                                        metPARA,grid_dryModel,gridRES_dryModel,...
%                                        grid_wetModel,gridRES_wetModel,GTropCoeff,VMFgrids) 
%    OR

%[STD,SHD,SWD,ZTD,ZHD,ZWD] = TropoModels(Time,ReceiverPos,ElevationAngle,dryModel,wetModel,...
%                                        metPARA,grid_dryModel,gridRES_dryModel,...
%                                        grid_wetModel,gridRES_wetModel,GTropCoeff,VMFgrids) 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%*******THE FUNCTION CALLS THE FOLLOWING SUBROUTINE(s):                    *
%*****INTERNAL                                                             *
%1.[X,Y,Z] = geo2xyz(latitude,longitude,height,RefEllipsoid);              *
%2.[latRAD,longRAD,h,latDEC,longDEC] = xyz2LLH(X,Y,Z,RefEllipsoid);        *
%3.[a,finv] = Elipsoidpara(RefEllipsoid);                                  *
%4.[Az_rad,El_rad,Az_deg,El_deg,D]=satAzimuthElevation(UserXYZ,SatXYZ,...  *
%                                                            RefEllipsoid);*
%5.[AVG,AMP] =interpolate_table(UserLat,avgValues,ampValues)               *
%6.[T,P,E,TM,lambda,ZHD,ZWD,ZTD]=getMETparaZTD_table(UTCtime,lat,hgt,model)*
%7.[p,T,e,Tm,la,dT,undu,ZHD,ZWD,ZTD]=getMETparaZTD_GPT2w(UTCtime,lat,...   *
%                                                        lon,hell,varargin)*
%8.[JD, MJD,DoY]=utc2JulianDay_DoY(Year,Month,Day,Hour,Minute,Seconds)     *  
%9.[FileName, FilePath] = searchGRID(varargin)                             *
%10.[ZHD,ZWD,ZTD] =ZenithTropDelay(T,P,es,Tm,lambda,dryModel,wetModel )    *
%11.[STD,SHD,SWD] =SlantTropDelay(T,P,es,Tm,lambda,satELEV,dryModel,...    *
%                                               wetModel,MFh_type,MFw_type)*
%*****EXTERNAL                                                             *
%1.[MFh, MFw]=MF(UTCtime,lat,lon,hgt,satEL,MFh_type,MFw_type,T,P,e)        *                                                                  
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%      Generally, the function accepts twelve(16) sets of inputs:          *
%1.        Time : Receiver reception time in[year,Month,Day,Hour,Min,Sec]  *
%--------------------------------------------------------------------------*
%*.    POSITIONs:                                                          *
%      =-=-=-=-=                                                           *
%2. ReceiverPos : Receiver position in either Latitude,Longitude &         *
%                     height or ECEF(XYZ) Coordinates                      *
%3.      SatPos : Satellite position(s) in ECEF(XYZ) Coordinates           *
%--------------------------------------------------------------------------
%*.       MODELS:                                                          *
%         =-=-=-=                                                          *
%4.     dryModel: Dry Tropospheric model in string(e.g.'Hopfield')         *
%5.     wetModel: Wet Tropospheric model in string(e.g.'Ifadis')           *
%--------------------------------------------------------------------------*
%*.    METEOROLOGICAL PARAMETER:                                           *
%      =-=-=-=-=-=-=-=-=-=-=-=-=                                           *
%6     metPARA contains the following as in the order:                     *
%6.1  Temperature(T): Atmospheric Temperature in Degree Celsius(C)         *
%6.2    Pressure(P) : Atmospheric Pressure in millibars(mbar /hPa)         *
%6.3              e : Water vapor pressure in [mbar / hpa]                 *
%OTHERS (Optional, more specifically for 'Askne & Nordius' wet model) may  *                                     *
%       also include:                                                      *
%6.4.            Tm : Mean temperature of water vapor in kelvin(k)         *
%6.5.        lambda : water vapour `lapse rate'(dimensionless)]            * 

%Meteorological parameters SHOULD BE ENTERED IN A MATRIX FORM AS:          *
%[T,P,RH,Tm,lambda]                                                        *

%NOTE1:1//. MET parameters should be entered in manner as indicated above  *
%           i.e. metPARA = [T,P,e,Tm,lambda]                               *
%NOTE2:                                                                    *
%      IF BOTH HYDROSTATATIC AND WET DELAY MODELs ARE GPT MODELs...        *
%      (i.e. GPT2,GPT2w,GPT3), OR the UNB3m,EGNOS & MOPS MODELS,then       *  
%      metPARA can be assign with empty([]) matrix: i.e. metPARA = []      *
%      OTHERWISE, provide the MET parameters                               *
%--------------------------------------------------------------------------*
%* GRID VALUES/FILES FOR GRID TROPO MODELS(GPT2,GPT2w,GPT3)*               *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-                *
%7. grid_dryModel : grid values in cells extracted from the grid file for  *
%                   HYDROSTATATIC DELAY MODEL                              *
%8.gridRES_dryModel: grid resolution (°) (possible: 1 or 5) for DRY model. *
%                     indicate 1 for 1°x1° external grid file and 5        *
%                     for 5°x5° external grid                              *
%9. grid_wetModel : grid values in cells extracted from the grid file for  *
%                    WET DELAY MODEL                                       *
%10.gridRES_wetModel : grid resolution (°) (possible: 1 or 5) for the WET  * 
%                      model.indicate 1 for 1°x1° external grid file and 5 *
%                      for 5°x5° external grid                             *
%11.Timevar_dry:0 i.e. 0 means GPT model with time variation, 1 means static    
%12.Timevar_wet:0 i.e. 0 means GPT model with time variation, 1 means static 
%13.GTrop_coeff : GTrop model coefficients / grid values extracted from ...*  
%                GTropCoefficient.mat 1°x 1° external grid file            * 
%NOTE1:                                                                    *          
%      The GTrop model coefficients / grid values are extracted            *
%      using the "SearchGTropgrid.m" sub-routine                           *
%14.VMFgrids     : Structure array of VMF grid files which can be :        *
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
%--------------------------------------------------------------------------*
%NOTE3:                                                                    *
%      For GPT MODELs(i.e. GPT2,GPT2w,GPT3):                               *
%1.****IF BOTH HYDROSTATATIC AND WET DELAY MODELs ARE GPT MODELs,provide   * 
%      inputs for the ff:                                                  *
%      grid_dryModel,gridRES_dryModel                                      *
%      grid_wetModel;gridRES_wetModel                                      *
%2.*** IF HYDROSTATATIC MODEL IS  GPT MODEL AND WET DELAY MODEL IS NOT     *
%      provide the ff inputs                                               *
%      grid_dryModel,gridRES_dryModel AND.....                             *
%      grid_wetModel = [];gridRES_wetModel = []                            *
%3.*** IF WET MODEL IS  GPT MODEL AND HYDROSTATATIC DELAY MODEL IS NOT     *
%      provide the ff inputs:                                              *
%      grid_wetModel,gridRES_wetModel     AND ....                         *
%      grid_dryModel = []; gridRES_dryModel = []                           *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
%NOTE:  In This case this sub-routine must go together with "readGPTgrid.m"*
%       which reads the grid and saves it in cell structures.              *
%       +The "readGPTgrid.m" should be ran first to provide the grid input *
%       file (grid)                                                        *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-**
%4.***ANY MODEL OTHER THE GPT MODELS,the corresponding grid model and grid *
%     resolution should be empty i.e. either or all of these should be     *
%     %empty([]):                                                          *
%     grid_wetModel;gridRES_wetModel                                       *
%     grid_dryModel; gridRES_dryModel   
%--------------------------------------------------------------------------*
%OTHER CONSIDERATIONs:                                                     *
%--------------------------------------------------------------------------+
%1.*******For only zenith tropospheric delays, set SatPos to empty([])     *

%2*******Any input parameter that is not relevant in the computation of    *
%        Tropospheric Delays based on the model type should be entered as  *
%        empty([]) matrix.For example is time is not needed based on the   *
%        model type(tropo model),then it should be entered as [].          *
%--------------------------------------------------------------------------+
%3.*******The function also accept Elevation angle as input for Satellite  *
%       position(s). i.e: SatPos : elevation angle in decimal degrees      * 
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
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%**************************************************************************+
%**************************************************************************+
function [STD,SHD,SWD,ZTD,ZHD,ZWD] = TropoModels(Time,ReceiverPos,SatPos,dryModel,wetModel,...
                                                 metPARA,grid_dryModel,gridRES_dryModel,...
                                                 grid_wetModel,gridRES_wetModel,Timevar_dry,Timevar_wet,GTropCoeff,VMFgrids)
                                                             
                                                                                                                                              
%**********CHECK INPUTs & REFORMAT INPUT DATA
switch nargin
   
    case {14,13,12,6} %Various inputs format
         
        %(1)*****************CHECK TIME INPUT
        %------------------------------------------------------------------
        %SOME MODELS REQUIRE TIME INPUT TO PERFORM / FUNCTION.IF TIME INPUT
        %IS EMPTY([])& ANY OF THESE MODELS IS PROVIDED,TERMINATE PROCESS
        %------------------------------------------------------------------
        if isempty(Time)
           %(1.1)*********CHECK DRY & WET TROPOSPHERIC MODELS
           %**********HYDROSTATIC
           if ~isempty(dryModel)
               
              if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                      strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),strncmpi(dryModel,'EGNOS',5),...
                      strncmpi(dryModel,'MOPS',4),strncmpi(dryModel,'VMF gridded ZHD',15),any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),strfind(dryModel,'GTrop')])])
                            
                 %ISSUE ERROR MESSAGE  
                 beep%Give a beep sound 
                 errmsg0{1}=sprintf('Observation / Reception Time is not provided i.e it''s Empty ([ ]).\n');
                 errmsg0{2}='Please provide Observation / reception Time & Try again.';
                 errordlg(errmsg0,'Time Input Error','modal')  
             
                 %RETURN EMPTY([]) DELAYS
                 STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
                 return
              end
              
           end %//if ~isempty(dryModel)
           
           %************WET
           if ~isempty(wetModel)
               
              if any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                      strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),strncmpi(wetModel,'EGNOS',5),...
                      strncmpi(wetModel,'MOPS',4),strncmpi(wetModel,'VMF gridded ZWD',15),any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),strfind(wetModel,'GTrop')])])
              
                 %ISSUE ERROR MESSAGE  
                 beep%Give a beep sound 
                 errmsg00{1}=sprintf('Observation / Reception Time is not provided i.e it''s Empty ([ ]).\n');
                 errmsg00{2}='Please provide Observation / reception Time & Try again.';
                 errordlg(errmsg00,'Time Input Error','modal')  
             
                 %RETURN EMPTY([]) DELAYS
                 STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
                 return
             
              end
              
           end %//if ~isempty(wetModel)
 
        end   %//if isempty(Time)
                   
      %(2)*****************CHECK RECEIVER/SITE POSITION/COORDINATE  INPUT
      %-------------------------------------------------------------------
      %SOME MODELS REQUIRE RECEIVER/SITE COORDs TO PERFORM / FUNCTION.IF 
      %RECEIVER POSITION INPUT IS EMPTY([])& ANY OF THESE MODELS IS 
      %PROVIDED,TERMINATE PROCESS
      %--------------------------------------------------------------------
      if isempty(ReceiverPos)
         %(2.1)*********CHECK DRY & WET TROPOSPHERIC MODELS
         %******HYDROSTATIC
         if ~isempty(dryModel)
            if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                    strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),strncmpi(dryModel,'Askne & Nordius',15),...
                    strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),strncmpi(dryModel,'Saastamoinen',12),strncmpi(dryModel,'Davis et al)',12),...
                    strncmpi(dryModel,'VMF gridded ZHD',15),any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),strfind(dryModel,'GTrop')])])
                 
               %ISSUE ERROR MESSAGE 
               beep%Give a beep sound 
               errmsg10{1}=sprintf('Reciever / Site Position is not provided i.e it''s Empty ([ ]).\n');
               errmsg10{2}='Please provide Receiver / Site position & Try again.';
               errordlg(errmsg10,'Reciever Position(s) Input Error','modal')  
           
               %RETURN EMPTY([]) DELAYS
               STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
               return
          
            end 
            
         end %//if ~isempty(dryModel)
         
         %***********WET
         if ~isempty(wetModel)
             
            if any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                    strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),...
                    strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4),strncmpi(wetModel,'Saastamoinen',12),strncmpi(wetModel,'Askne & Nordius',15),...
                    strncmpi(wetModel,'VMF gridded ZWD',15),any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),strfind(wetModel,'GTrop')])])
             
               %ISSUE ERROR MESSAGE 
               beep%Give a beep sound 
               errmsg11{1}=sprintf('Reciever / Site Position is not provided i.e it''s Empty ([ ]).\n');
               errmsg11{2}='Please provide Receiver / Site position & Try again.';
               errordlg(errmsg11,'Reciever Position(s) Input Error','modal')  
           
               %RETURN EMPTY([]) DELAYS
               STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
               return
          
            end
            
         end %//if ~isempty(wetModel)
              
         else 
            %IF RECEIVER POSITION IS NOT EMPTY,BUT ONLY HEIGHTS ARE PROVIDED

             if size(ReceiverPos,2)==1 
               
                %CHECK MODELS THAT REQUIRE ONLY STATION HEIGHT TO PERFORM
                %******HYDROSTATIC
                
                if ~isempty(dryModel)
                    
                   if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                           strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),strncmpi(dryModel,'EGNOS',5),...
                           strncmpi(dryModel,'MOPS',4),strncmpi(dryModel,'Saastamoinen',12),strncmpi(dryModel,'Davis et al)',12),strncmpi(dryModel,'Askne & Nordius',15),...
                           strncmpi(dryModel,'VMF gridded ZHD',15),any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),strfind(dryModel,'GTrop')])]) 
                      
                      %ISSUE ERROR MESSAGE for HEIGHT INPUT 
                      beep%Give a beep sound 
                      errmsg50{1}=sprintf('Insuficient Input for Receiver / Station Position .\n');
                      errmsg50{2}=sprintf('Latitude & Longitude Coordinates have not been provided.\n');
                      errmsg50{3}=sprintf('Please provide Latitude & Longitude Coordinates & Try again.');
                      errordlg(errmsg50,'Coordinate(s) Input Error','modal')  
                      return
                   
                   end  %//if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),...
                   
                end %//if ~isempty(dryModel)
                
                
               %***********WET
               if ~isempty(wetModel)
                   
                  if any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                          strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),strncmpi(wetModel,'EGNOS',5),...
                          strncmpi(wetModel,'MOPS',4),strncmpi(wetModel,'Saastamoinen',12),strncmpi(wetModel,'Askne & Nordius',15),...
                          strncmpi(wetModel,'VMF gridded ZWD',15),any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),strfind(wetModel,'GTrop')])]) 
                
                     %ISSUE ERROR MESSAGE for HEIGHT INPUT 
                     beep%Give a beep sound 
                     errmsg51{1}=sprintf('Insuficient Input for Receiver / Station Position .\n');
                     errmsg51{2}=sprintf('Latitude & Longitude Coordinates have not been provided.\n');
                     errmsg51{3}=sprintf('Please provide Latitude & Longitude Coordinates & Try again.');
                     errordlg(errmsg51,'Coordinate(s) Input Error','modal')  
                     return
                   
                  end   %//if any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),...
                  
               end %//if ~isempty(wetModel)
                
             end %//if size(ReceiverPos,2)==1 

      end %//if isempty(ReceiverPos)
      
      %(3)***************CHECK HYDROSTATIC TROPOPHERIC DELAY MODEL  INPUT
      %--------------------------------------------------------------------
      %IF  HYDROSTATIC DELAY MODEL IS EMPTY([]),ASSIGN DEFAULT MODEL      
      %--------------------------------------------------------------------
      if isempty(dryModel)
         
         if ~isempty(wetModel) 
             
            if any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                    strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),...
                    strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4),strncmpi(wetModel,'Saastamoinen',12),strncmpi(wetModel,'Hopfield',8),...
                    strncmpi(wetModel,'Modified Hopfield(Goads & Goodman)',34),strncmpi(wetModel,'Askne & Nordius',15),strncmpi(wetModel,'VMF gridded ZWD',15),...
                    any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),strfind(wetModel,'GTrop')])]) 
            
                %SET DEFAULT HYDROSTATIC DELAY MODEL USING WET MODEL
                if strncmpi(wetModel,'VMF gridded ZWD',15)
                   dryModel = 'VMF gridded ZHD';
                else
                   dryModel = wetModel;
                end
                
                %PRINT SOME MESSAGE TO USER
                fprintf('Hydrostatic Correction Delay model is not provided i.e it''s Empty ([ ]).\n Using  %s as default Hydrostatic Correction model.',dryModel)     
            else
                dryModel = 'Saastamoinen';
               
                %PRINT SOME MESSAGE TO USER
                fprintf('Hydrostatic Correction Delay model is not provided i.e it''s Empty ([ ]).\n Using  %s as default Hydrostatic Correction model.',dryModel)
                 
            end
            
         else
               %SET DEFAULT HYDROSTATIC DELAY MODEL
               dryModel = 'Saastamoinen';
               
                %PRINT SOME MESSAGE TO USER
                fprintf('Hydrostatic Correction Delay model is not provided i.e it''s Empty ([ ]).\n Using  %s as default Hydrostatic Correction model.',dryModel)     
         end
         
      end
         
      %(4)***************CHECK WET TROPOPHERIC DELAY MODEL  INPUT
      %--------------------------------------------------------------------
      %IF  WET DELAY MODEL IS EMPTY([]),ASSIGN DEFAULT MODEL      
      %--------------------------------------------------------------------
      if isempty(wetModel)  
         
         if ~isempty(dryModel) 
             
            if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                    strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),...
                    strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),strncmpi(dryModel,'Saastamoinen',12),strncmpi(dryModel,'Hopfield',8),...
                    strncmpi(dryModel,'Modified Hopfield(Goads & Goodman)',34),strncmpi(dryModel,'Askne & Nordius',15),strncmpi(dryModel,'VMF gridded ZHD',15),...
                    any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),strfind(dryModel,'GTrop')])]) 
            
                %SET DEFAULT HYDROSTATIC DELAY MODEL USING WET MODEL
                if strncmpi(dryModel,'VMF gridded ZHD',15)
                   wetModel = 'VMF gridded ZWD';
                else
                    wetModel = dryModel;
                end
                
                %PRINT SOME MESSAGE TO USER
                fprintf('Wet Correction Delay model is not provided i.e it''s Empty ([ ]).\n Using  %s as default Wet Correction model.',wetModel)     
            else
                wetModel = 'Saastamoinen';
               
                %PRINT SOME MESSAGE TO USER
                fprintf('Wet Correction Delay model is not provided i.e it''s Empty ([ ]).\n Using  %s as default Wet Correction model.',wetModel)
                 
            end
            
         else
               %SET DEFAULT HYDROSTATIC DELAY MODEL
               wetModel = 'Saastamoinen';
               
                %PRINT SOME MESSAGE TO USER
                fprintf('Wet Correction Delay model is not provided i.e it''s Empty ([ ]).\n Using  %s as default Hydrostatic Correction model.',wetModel)     
         end
         
      end %//if isempty(wetModel) 
           
      %(5)**************CHECK METEOROLOGICAL PARAMETERS  INPUT
      %--------------------------------------------------------------------
      if ~isempty(metPARA)
         if size(metPARA,2)==5 %if the number of columns = 5
            T = metPARA(:,1);
            P = metPARA(:,2);
           es = metPARA(:,3);
           Tm = metPARA(:,4);
       lambda = metPARA(:,5);
             
         elseif size(metPARA,2)==4 %if the number of columns = 4
                T = metPARA(:,1);
                P = metPARA(:,2);
               es = metPARA(:,3);
               Tm = metPARA(:,4);
           lambda = []; %Assign empty([]) matrix to lambda
               
               %NOTE:'Askne & Nordius' TROPO WET DELAY MODEL REQUIRES water 
               %vapour lapse rate ( lambda ) TO PERFORM.
               if ~isempty(wetModel) 
                   
                  if strncmpi(wetModel,'Askne & Nordius',15)
                      
                     %ISSUE ERROR MESSAGE  
                     beep%Give a beep sound 
                     errmsg1{1}=sprintf('Insuficient Input for Meteorological Parameters.\n');
                     errmsg1{2}=sprintf('Water Vapour Lapse Rate ( lambda ) has not been provided.\n');
                     errmsg1{3}=sprintf('Please provide water vapour lapse rate (lambda) & Try again.');
                     errordlg(errmsg1,'Meteorological parameters Input Error','modal') 
                
                  %RETURN EMPTY([]) DELAYS
                  STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           
                  return 
                  
                  end
               
               end %//if ~isempty(wetModel) 
                   
                  
         elseif size(metPARA,2)==3 %if the number of columns = 3
                T = metPARA(:,1);
                P = metPARA(:,2);
                es = metPARA(:,3);
                Tm = []; %Assign empty([]) matrix to Tm
            lambda = []; %Assign empty([]) matrix to lambda
                
              if all([~isempty(dryModel),~isempty(wetModel)])
                
                 if any([all([~any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),strncmpi(dryModel,'VMF gridded ZHD',15),...
                                    strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),...
                                    strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),strfind(dryModel,'GTrop')])]),...
                              ~any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),strncmpi(wetModel,'VMF gridded ZWD',15),...
                                    strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),...
                                    strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4),any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),strfind(wetModel,'GTrop')])])]),...
                          any([~any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),strncmpi(dryModel,'VMF gridded ZHD',15),...
                                     strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),...
                                     strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),strfind(dryModel,'GTrop')])]),...
                               ~any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),strncmpi(wetModel,'VMF gridded ZWD',15)...
                                     strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),...
                                     strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4),any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),strfind(wetModel,'GTrop')])])])])
               
                    %NOTE:'Askne & Nordius' TROPO WET DELAY MODEL REQUIRES water 
                    %vapour lapse rate ( lambda ) TO PERFORM.
                     if strncmpi(wetModel,'Askne & Nordius',15)
                    
                         %ISSUE ERROR MESSAGE  
                         beep%Give a beep sound 
                         errmsg2{1}=sprintf('Insuficient Inputs for Meteorological Parameters.\n');
                         errmsg2{2}=sprintf('The following parameters have not been provided :.\n');
                         errmsg2{3}=sprintf('1. Mean Temperature of Water Vapor ( Tm ).\n');
                         errmsg2{4}=sprintf('2. Water Vapour Lapse Rate ( lambda ).\n');
                         errmsg2{5}=sprintf('Please provide these Meteorological parameters & Try again.');
                         errordlg(errmsg2,'Meteorological parameters Input Error','modal') 
                
                         %RETURN EMPTY([]) DELAYS
                         STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
                  
                         return
                    
                      end

                end 
                
              end %//if all([~isempty(dryModel),~isempty(dryModel)])
                         
         elseif size(metPARA,2)==2 %if the number of columns = 2
                T = metPARA(:,1);
                P = metPARA(:,2);
               es = [];%Assign empty([]) matrix to es
               Tm = []; %Assign empty([]) matrix to Tm
           lambda = []; %Assign empty([]) matrix to lambda
           
              if all([~isempty(dryModel),~isempty(wetModel)])
                
                if any([all([~any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),strncmpi(dryModel,'VMF gridded ZHD',15),...
                                   strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),...
                                   strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),strfind(dryModel,'GTrop')])]),...
                             ~any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),strncmpi(wetModel,'VMF gridded ZWD',15),...
                                   strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),...
                                   strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4),any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),strfind(wetModel,'GTrop')])])]),...
                          any([~any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),strncmpi(dryModel,'VMF gridded ZHD',15),...
                                     strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),...
                                     strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),strfind(dryModel,'GTrop')])]),...
                               ~any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),strncmpi(wetModel,'VMF gridded ZWD',15),...
                                     strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),...
                                     strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4),any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),strfind(wetModel,'GTrop')])])])])
               
                        
                       if ~strncmpi(wetModel,'Askne & Nordius',15)
                    
                          %ISSUE ERROR MESSAGE  
                          beep%Give a beep sound 
                          errmsg3{1}=sprintf('Insuficient Input for Meteorological Parameters.\n');
                          errmsg3{2}=sprintf('Water Vapour Partial Pressure (e) has not been provided.\n');
                          errmsg3{3}=sprintf('Please provide water vapour partial pressure & Try again.');
                          errordlg(errmsg3,'Meteorological parameters Input Error','modal')
                  
                          %RETURN EMPTY([]) DELAYS
                          STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
                  
                          return
                 
                       else  
                           %ISSUE ERROR MESSAGE  
                           beep%Give a beep sound 
                           errmsg4{1}=sprintf('Insuficient Inputs for Meteorological Parameters.\n');
                           errmsg4{2}=sprintf('The following parameters have not been provided :.\n');
                           errmsg4{3}=sprintf('1. Water Vapour Partial Pressure ( e ) .\n');
                           errmsg4{4}=sprintf('2. Mean Temperature of Water Vapor ( Tm ).\n');
                           errmsg4{5}=sprintf('3. Water Vapour Lapse Rate ( lambda ).\n');
                           errmsg4{6}=sprintf('Please these parameters & Try again.');
                           errordlg(errmsg4,'Meteorological parameters Input Error','modal') 
                    
                           %RETURN EMPTY([]) DELAYS
                           STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
                  
                           return
                    
                       end
                       
                end
                    
              end %//if all([~isempty(dryModel),~isempty(dryModel)])
               
                    
         elseif size(metPARA,2)==1 %if the number of columns = 1
                T = metPARA(:,1);
                P = [];%Assign empty([]) matrix to P
               es = [];%Assign empty([]) matrix to es
               Tm = []; %Assign empty([]) matrix to Tm
           lambda = []; %Assign empty([]) matrix to lambda
                
               if all([~isempty(dryModel),~isempty(wetModel)])
                   
                if any([all([~any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),strncmpi(dryModel,'VMF gridded ZHD',15),...
                                   strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),...
                                   strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),strfind(dryModel,'GTrop')])]),...
                             ~any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),strncmpi(wetModel,'VMF gridded ZWD',15),...
                                   strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),...
                                   strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4),any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),strfind(wetModel,'GTrop')])])]),...
                          any([~any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),strncmpi(dryModel,'VMF gridded ZHD',15),...
                                     strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),...
                                     strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),strfind(dryModel,'GTrop')])]),...
                               ~any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),strncmpi(wetModel,'VMF gridded ZWD',15),...
                                     strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),...
                                     strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4),any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),strfind(wetModel,'GTrop')])])])])
               
                     if ~strncmpi(wetModel,'Askne & Nordius',15)
                    
                        %ISSUE ERROR MESSAGE  
                        beep%Give a beep sound 
                        errmsg5{1}=sprintf('Insuficient Inputs for Meteorological Parameters.\n');
                        errmsg5{2}=sprintf('The following parameters have not been provided :.\n');
                        errmsg5{3}=sprintf('1. Surface Pressure ( P ) .\n');
                        errmsg5{4}=sprintf('2. Water Vapour Partial Pressure ( e ) .\n');
                        errmsg5{5}=sprintf('Please provide these meteorological parameters & Try again.');
                        errordlg(errmsg5,'Meteorological parameters Input Error','modal') 
                  
                        %RETURN EMPTY([]) DELAYS
                        STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
                  
                        return
                 
                     else  
                         %ISSUE ERROR MESSAGE  
                         beep%Give a beep sound 
                         errmsg6{1}=sprintf('Insuficient Inputs for Meteorological Parameters.\n');
                         errmsg6{2}=sprintf('The following parameters have not been provided :.\n');
                         errmsg6{3}=sprintf('1. Surface Pressure ( P ) .\n');
                         errmsg6{4}=sprintf('2. Water Vapour Partial Pressure ( e ) .\n');
                         errmsg6{5}=sprintf('2. Mean Temperature of Water Vapor ( Tm ).\n');
                         errmsg6{6}=sprintf('3. Water Vapour Lapse Rate ( lambda ).\n');
                         errmsg6{7}=sprintf('Please provide these Meteorological parameters & Try again.');
                         errordlg(errmsg6,'Meteorological parameters Input Error','modal') 
                    
                         %RETURN EMPTY([]) DELAYS
                         STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
                  
                         return
                    
                     end

                end
                
               end %//if all([~isempty(dryModel),~isempty(wetModel)])
                
         end %//if size(metPARA,2)==5 %if the number of columns = 5
                    
      else
          T = [];%Assign empty([]) matrix to T
          P = [];%Assign empty([]) matrix to P
         es = [];%Assign empty([]) matrix to es
         Tm = []; %Assign empty([]) matrix to Tm
     lambda = []; %Assign empty([]) matrix to lambda
     
          if all([~isempty(dryModel),~isempty(wetModel)])
              
           if any([all([~any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),strncmpi(dryModel,'VMF gridded ZHD',15),...
                                   strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),...
                                   strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),strfind(dryModel,'GTrop')])]),...
                             ~any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),strncmpi(wetModel,'VMF gridded ZWD',15),...
                                   strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),...
                                   strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4),any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),strfind(wetModel,'GTrop')])])]),...
                          any([~any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),strncmpi(dryModel,'VMF gridded ZHD',15),...
                                     strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),...
                                     strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),strfind(dryModel,'GTrop')])]),...
                               ~any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),strncmpi(wetModel,'VMF gridded ZWD',15),...
                                     strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),...
                                     strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4),any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),strfind(wetModel,'GTrop')])])])])
               
                                     
              if ~strncmpi(wetModel,'Askne & Nordius',15)
                    
                  %ISSUE ERROR MESSAGE  
                  beep%Give a beep sound 
                  errmsg7{1}=sprintf('Insuficient Inputs for Meteorological Parameters.\n');
                  errmsg7{2}=sprintf('The following parameters have not been provided :.\n');
                  errmsg7{3}=sprintf('1. Surface Temperature ( T ) .\n');
                  errmsg7{4}=sprintf('2. Surface Pressure ( P ) .\n');
                  errmsg7{5}=sprintf('3. Water Vapour Partial Pressure ( e ) .\n');
                  errmsg7{6}=sprintf('Please provide these meteorological parameters & Try again.');
                  errordlg(errmsg7,'Meteorological parameters Input Error','modal') 
                  
                  %RETURN EMPTY([]) DELAYS
                  STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
                  
                  return
                 
              else
                  beep%Give a beep sound 
                  errmsg8{1}=sprintf('Insuficient Inputs for Meteorological Parameters.\n');
                  errmsg8{2}=sprintf('The following parameters have not been provided :.\n');
                  errmsg8{3}=sprintf('1. Surface Temperature ( T ) .\n');
                  errmsg8{4}=sprintf('2. Surface Pressure ( P ) .\n');
                  errmsg8{5}=sprintf('3. Water Vapour Partial Pressure ( e ) .\n');
                  errmsg8{6}=sprintf('4. Mean Temperature of Water Vapor ( Tm ).\n');
                  errmsg8{7}=sprintf('5. Water Vapour Lapse Rate ( lambda ).\n');
                  errmsg8{8}=sprintf('Please provide these meteorological parameters & Try again.');
                  errordlg(errmsg8,'Meteorological parameters Input Error','modal')   
                    
                  %RETURN EMPTY([]) DELAYS
                  STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
                  
                  return
                  
              end 
                    
           end 
           
          end %//if all([~isempty(dryModel),~isempty(wetModel)])
              
      end %//if ~isempty(metPARA)
      
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
      %*****************WHEN ONLY 6 INPUTS ARE PROVIDED
      %NOTE: SOME MODELS DO NOT REQUIRE GRID FILES & RESOLUTIONS TO PERFORM
      %IF THOSE MODELS ARE NOT PROVIDED, AN ERROR MESSAGE IS ISSUED
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
      if nargin == 6 %FOR 6 INPUTS
          
         %***********CHECK DRY & WET TROPOSPHERIC DELAY MODEL TYPE
         %NOTE:IF TROPO MODEL TYPE IS ONE OF THE VMF MODELS WHICH REQUIRE 
         %GRID FILES & RESOLUTION.
         %----------------------------------------------------------------- 
         %******HYDROSTATIC
         if ~isempty(dryModel)
             
            if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),...
                    strncmpi(dryModel,'GPT2w (5° x 5°)',14), strncmpi(dryModel,'GPT3 (1° x 1°)',14),...
                    strncmpi(dryModel,'GPT3 (5° x 5°)',14)])
                    
               %**********ISSUE ERROR MESSAGE   
               if strncmpi(dryModel,'GPT2 (5° x 5°)',14)
                   
                  beep%Give a beep sound 
                  errmsg16{1}=sprintf('No input for GPT2 grid values.\n');
                  errmsg16{2}='Please provide GPT2 grid values & Try again.';
                  errordlg(errmsg16,'GPT2 Grid values Input Error','modal') 
                  
               elseif any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14)])
                   
                      beep%Give a beep sound 
                      errmsg161{1}=sprintf('No input for GPT2w grid values & grid resolution.\n');
                      errmsg161{2}='Please provide GPT2w grid values & grid resolution & Try again.';
                      errordlg(errmsg161,'GPT2w Grid values & resolution Input Error','modal') 
                
               elseif any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)])
                   
                      beep%Give a beep sound 
                      errmsg162{1}=sprintf('No input for GPT3 grid values & grid resolution.\n');
                      errmsg162{2}='Please provide GPT3 grid values & grid resolution & Try again.';
                      errordlg(errmsg162,'GPT3 Grid values & resolution Input Error','modal')
               end
                
               %RETURN EMPTY([]) DELAYS
               STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           
               return 
               
               elseif strncmpi(dryModel,'VMF gridded ZHD',15)
                
                   %ISSUE ERROR MESSAGE  
                   beep%Give a beep sound 
                   errmsg163{1}=sprintf('No input for VMF grid files .\n');
                   errmsg163{2}='Please provide VMF grid files & Try again.';
                   errordlg(errmsg163,'VMF grid files Input Error','modal') 
                   
                   %RETURN EMPTY([]) DELAYS
                   STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           
                   return 
                   
            elseif any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),strfind(dryModel,'GTrop')])
                  
                   %ISSUE ERROR MESSAGE  
                   beep%Give a beep sound 
                   errmsg164{1}=sprintf('No input for GTrop grid values or Coefficients.\n');
                   errmsg164{2}='Please provide GTrop Coefficients & Try again.';
                   errordlg(errmsg164,'GTrop Coefficients Input Error','modal')  
                   
                   %RETURN EMPTY([]) DELAYS
                   STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           
                   return 
            end  
            
         end %//if ~isempty(dryModel)
         
         %********WET
         if ~isempty(wetModel)
             
            if any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),...
                    strncmpi(wetModel,'GPT2w (5° x 5°)',14),strncmpi(wetModel,'GPT3 (1° x 1°)',14),...
                    strncmpi(wetModel,'GPT3 (5° x 5°)',14)])
              
               %**********ISSUE ERROR MESSAGE   
               if strncmpi(wetModel,'GPT2 (5° x 5°)',14)
                   
                  beep%Give a beep sound 
                  errmsg165{1}=sprintf('No input for GPT2 grid values.\n');
                  errmsg165{2}='Please provide GPT2 grid values & Try again.';
                  errordlg(errmsg165,'GPT2 Grid values Input Error','modal') 
                  
               elseif any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14)])
                   
                      beep%Give a beep sound 
                      errmsg166{1}=sprintf('No input for GPT2w grid values & grid resolution.\n');
                      errmsg166{2}='Please provide GPT2w grid values & grid resolution & Try again.';
                      errordlg(errmsg166,'GPT2w Grid values & resolution Input Error','modal') 
                
               elseif any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])
                   
                      beep%Give a beep sound 
                      errmsg167{1}=sprintf('No input for GPT3 grid values & grid resolution.\n');
                      errmsg167{2}='Please provide GPT3 grid values & grid resolution & Try again.';
                      errordlg(errmsg167,'GPT3 Grid values & resolution Input Error','modal')
               end  
                
               %RETURN EMPTY([]) DELAYS
               STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           
               return
               
            elseif strncmpi(wetModel,'VMF gridded ZWD',15)
                
                   %ISSUE ERROR MESSAGE  
                   beep%Give a beep sound 
                   errmsg168{1}=sprintf('No input for VMF grid files .\n');
                   errmsg168{2}='Please provide VMF grid files & Try again.';
                   errordlg(errmsg168,'VMF grid files Input Error','modal') 
                   
                   %RETURN EMPTY([]) DELAYS
                   STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           
                   return 
                   
            elseif any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),strfind(wetModel,'GTrop')])
                  
                   %ISSUE ERROR MESSAGE  
                   beep%Give a beep sound 
                   errmsg169{1}=sprintf('No input for GTrop grid values or Coefficients.\n');
                   errmsg169{2}='Please provide GTrop Coefficients & Try again.';
                   errordlg(errmsg169,'GTrop Coefficients Input Error','modal')  
                   
                   %RETURN EMPTY([]) DELAYS
                   STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           
                   return 
             
            end
            
         end %//if ~isempty(wetModel)
          
     %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*
      %*****************WHEN ALL 10 INPUTS ARE PROVIDED
      %NOTE: SOME MODELS REQUIRE GRID FILES & RESOLUTIONS TO PERFORM
      %IF THOSE MODELS ARE PROVIDED, CHECK THE INPUTS FORMATS AND ISSUE AN
      %ERROR MESSAGE IF ANY IS FOUND EMPTY([]) 
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*  
     %****************CHECH IF ALL INPUTS ARE PROVIDED
      elseif any(nargin == [14,13,12]) %FOR 16 INPUTS
             
            %1.CHECK GRID FILE FOR HYDROSTATIC TROPO MODELS
            if any([all([isempty(grid_dryModel),isempty(gridRES_dryModel)]),...
                    any([isempty(grid_dryModel),isempty(gridRES_dryModel)])])
               
               %IF PROVIDED MODEL(s) IS ANY OF THE GPT DRY MODELS,ISSUE MSG
               if ~isempty(dryModel)
                   
                  if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),...
                          strncmpi(dryModel,'GPT2w (5° x 5°)',14),strncmpi(dryModel,'GPT3 (1° x 1°)',14),...
                          strncmpi(dryModel,'GPT3 (5° x 5°)',14)])
                   
                     if isempty(grid_dryModel)& isempty(gridRES_dryModel)

                        %ISSUE ERROR MESSAGE
                         if strncmpi(dryModel,'GPT2 (5° x 5°)',14)
                   
                            beep%Give a beep sound 
                            errmsg20{1}=sprintf('No input for GPT2 grid values & grid resolution to compute ...\n');
                            errmsg20{2}=sprintf('Hydrostatic Tropospheric Delay, i.e. inputs are empty ( [ ] ).\n');
                            errmsg20{2}='Please provide GPT2 grid values & grid resolution & Try again.';
                            errordlg(errmsg20,'GPT2 grid values & resolution Input Error','modal') 
                         
                         elseif any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14)])
                         
                                beep%Give a beep sound 
                                errmsg21{1}=sprintf('No input for GPT2w grid values & grid resolution to compute ...\n');
                                errmsg21{2}=sprintf('Hydrostatic Tropospheric Delay, i.e. inputs are empty ( [ ] ).\n');
                                errmsg21{2}='Please provide GPT2w grid values & grid resolution & Try again.';
                                errordlg(errmsg21,'GPT2w grid values & resolution Input Error','modal') 
                         
                         elseif any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)])
                         
                                beep%Give a beep sound 
                                errmsg22{1}=sprintf('No input for GPT3 grid values & grid resolution to compute ...\n');
                                errmsg22{2}=sprintf('Hydrostatic Tropospheric Delay, i.e. inputs are empty ( [ ] ).\n');
                                errmsg22{2}='Please provide GPT3 grid values & grid resolution & Try again.';
                                errordlg(errmsg22,'GPT3 grid values & resolution Input Error','modal') 
                                
                         end 
                     
                         %RETURN EMPTY([]) DELAYS
                         STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           
                         return 
                     
                  elseif isempty(grid_dryModel)& ~isempty(gridRES_dryModel)
                           
                         %**********ISSUE ERROR MESSAGE   
                         if strncmpi(dryModel,'GPT2 (5° x 5°)',14)
                   
                            beep%Give a beep sound 
                            errmsg23{1}=sprintf('No input for GPT2 grid values to compute Hydrostatic ...\n');
                            errmsg23{2}=sprintf('Tropospheric Delay, i.e. input is empty ( [ ] ).\n');
                            errmsg23{2}='Please provide GPT2 grid values & Try again.';
                            errordlg(errmsg23,'GPT2 grid values Input Error','modal') 
                  
                         elseif any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14)])
                   
                                beep%Give a beep sound 
                                errmsg24{1}=sprintf('No input for GPT2w grid values to compute Hydrostatic ...\n');
                                errmsg24{2}=sprintf('Tropospheric Delay, i.e. input is empty ( [ ] ).\n');
                                errmsg24{2}='Please provide GPT2w grid values & Try again.';
                                errordlg(errmsg24,'GPT2w grid values Input Error','modal') 
                
                         elseif any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)])
                   
                                beep%Give a beep sound 
                                errmsg25{1}=sprintf('No input for GPT3 grid values to compute Hydrostatic ...\n');
                                errmsg25{2}=sprintf('Tropospheric Delay, i.e. input is empty ( [ ] ).\n');
                                errmsg25{2}='Please provide GPT3 grid values & Try again.';
                                errordlg(errmsg25,'GPT3 grid values Input Error','modal')
                         end  
                     
                         %RETURN EMPTY([]) DELAYS
                         STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           
                         return 
                         
                  elseif ~isempty(grid_dryModel)& isempty(gridRES_dryModel)
                      
                         %ISSUE ERROR MESSAGE  
                         beep%Give a beep sound 
                         errmsg26{1}=sprintf('No input for GPT grid resolution to compute Hydrostatic ...\n');
                         errmsg26{2}=sprintf('Tropospheric Delay, i.e. input is empty ( [ ] ).\n');
                         errmsg26{3}='Please provide GPT grid resolution & Try again.';
                         errordlg(errmsg26,'GPT grid resolution Input Error','modal') 
                     
                         %RETURN EMPTY([]) DELAYS
                         STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           
                         return 
                         
                     end  %//if isempty(grid_dryModel)& isempty(gridRES_dryModel)
                  
                  end  %//if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),...
               
               end  %//if ~isempty(dryModel)
               
            end  %//if any([all([isempty(grid_dryModel),isempty(gridRES_dryModel)]),...
            
            
           %2.CHECK GRID FILE FOR WET TROPO MODELS
            if any([all([isempty(grid_wetModel),isempty(gridRES_wetModel)]),...
                    any([isempty(grid_wetModel),isempty(gridRES_wetModel)])])
               
              if ~isempty(wetModel)
                  
                 %IF PROVIDED MODEL(s) IS ANY OF THE GPT WET MODELS,ISSUE MSG
                 if any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),...
                         strncmpi(wetModel,'GPT2w (5° x 5°)',14),strncmpi(wetModel,'GPT3 (1° x 1°)',14),...
                         strncmpi(wetModel,'GPT3 (5° x 5°)',14)])
                   
                    if isempty(grid_wetModel)& isempty(gridRES_wetModel)

                       %ISSUE ERROR MESSAGE
                       if strncmpi(wetModel,'GPT2 (5° x 5°)',14)
                   
                          beep%Give a beep sound 
                          errmsg27{1}=sprintf('No input for GPT2 grid values & grid resolution to compute ...\n');
                          errmsg27{2}=sprintf('Wet Tropospheric Delay, i.e. inputs are empty ( [ ] ).\n');
                          errmsg27{2}='Please provide GPT2 grid values & grid resolution & Try again.';
                          errordlg(errmsg27,'GPT2 grid values & resolution Input Error','modal') 
                         
                       elseif any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14)])
                         
                              beep%Give a beep sound 
                              errmsg28{1}=sprintf('No input for GPT2w grid values & grid resolution to compute ...\n');
                              errmsg28{2}=sprintf('Wet Tropospheric Delay, i.e. inputs are empty ( [ ] ).\n');
                              errmsg28{2}='Please provide GPT2w grid values & grid resolution & Try again.';
                              errordlg(errmsg28,'GPT2w grid values & resolution Input Error','modal') 
                         
                         elseif any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])
                         
                                beep%Give a beep sound 
                                errmsg29{1}=sprintf('No input for GPT3 grid values & grid resolution to compute ...\n');
                                errmsg29{2}=sprintf('Wet Tropospheric Delay, i.e. inputs are empty ( [ ] ).\n');
                                errmsg29{2}='Please provide GPT3 grid values & grid resolution & Try again.';
                                errordlg(errmsg29,'GPT3 grid values & resolution Input Error','modal') 
                                
                       end  
                     
                       %RETURN EMPTY([]) DELAYS
                       STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           
                       return 
                     
                  elseif isempty(grid_wetModel)& ~isempty(gridRES_wetModel)
                              
                         %**********ISSUE ERROR MESSAGE   
                         if strncmpi(wetModel,'GPT2 (5° x 5°)',14)
                   
                            beep%Give a beep sound 
                            errmsg30{1}=sprintf('No input for GPT2 grid values to compute Wet ...\n');
                            errmsg30{2}=sprintf('Tropospheric Delay, i.e. input is empty ( [ ] ).\n');
                            errmsg30{2}='Please provide GPT2 grid values & Try again.';
                            errordlg(errmsg30,'GPT2 grid values Input Error','modal') 
                  
                         elseif any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14)])
                   
                                beep%Give a beep sound 
                                errmsg31{1}=sprintf('No input for GPT2w grid values to compute Wet ...\n');
                                errmsg31{2}=sprintf('Tropospheric Delay, i.e. input is empty ( [ ] ).\n');
                                errmsg31{2}='Please provide GPT2w grid values & Try again.';
                                errordlg(errmsg31,'GPT2w grid values Input Error','modal') 
                
                         elseif any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])
                   
                                beep%Give a beep sound 
                                errmsg32{1}=sprintf('No input for GPT3 grid values to compute Wet ...\n');
                                errmsg32{2}=sprintf('Tropospheric Delay, i.e. input is empty ( [ ] ).\n');
                                errmsg32{2}='Please provide GPT3 grid values & Try again.';
                                errordlg(errmsg32,'GPT3 grid values Input Error','modal')
                         end  
                     
                         %RETURN EMPTY([]) DELAYS
                         STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           
                         return 
                         
                  elseif ~isempty(grid_wetModel)& isempty(gridRES_wetModel)
                      
                         %ISSUE ERROR MESSAGE  
                         beep%Give a beep sound 
                         errmsg33{1}=sprintf('No input for GPT grid resolution to compute Wet ...\n');
                         errmsg33{2}=sprintf('Tropospheric Delay, i.e. input is empty ( [ ] ).\n');
                         errmsg33{3}='Please provide GPT grid resolution & Try again.';
                         errordlg(errmsg33,'GPT grid resolution Input Error','modal') 
                     
                         %RETURN EMPTY([]) DELAYS
                         STD =[];SHD = [];SWD = [];ZTD = [];ZHD = [];ZWD = [];
           
                         return 
                         
                    end  %//if isempty(grid_wetModel)& isempty(gridRES_wetModel)
                  
                 end  %//if any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),...
               
              end  %//if ~isempty(wetModel)
               
            end  %//if any([all([isempty(grid_wetModel),isempty(gridRES_wetModel)]),... 
           
            %****CHECK TIME VARIATION INPUTS, IF EMPTY([]) ASSIGN DEFAULTs
          if isempty( Timevar_dry)
              
             %GET TIME VARIATION STATUS FROM goGPS GUI
             Timevar_dry  = getappdata(0,'Timevar_dry');%time variation 
             
             if isempty( Timevar_dry)
                Timevar_dry = 0;%with time variation (annual and semiannual terms)
             end
             
          end
            
          if isempty( Timevar_wet) 
             
              %GET TIME VARIATION STATUS FROM goGPS GUI
              Timevar_wet  = getappdata(0,'Timevar_wet');%time variation 
              
              if isempty( Timevar_dry)
                 Timevar_wet = 0;%with time variation (annual and semiannual terms)
              end  
          end
         
          
        if nargin == 13
            
           VMFgrids = [];%Assign empty([]) matrix to VMFgrids
           
        elseif nargin == 12
            
               GTropCoeff = [];%Assign empty([]) matrix to GTropCoeff
               VMFgrids = []; %Assign empty([]) matrix to VMFgrids
            
        end
      
           
      end  %// if nargin == 6  
                  
          %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-          
          
          %CHECK IF SATELLITE ELEVATION ANGLE & HEIGHT OF RECEIVER ARE...
          %                                               ENTERED INSTEAD
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
                Satpos=SatPos;%Satellite position in XYZ                   
             end   %\\if size(SatPos,2)==1
             
           else
               satpos_empty=1; %flag to indicate that satpos is empty([])
               setappdata(0,'satpos_empty',satpos_empty) 
               
           end %\\if (~isempty(SatPos)| ~isnan(SatPos))
           
           %2.RECEIVER POSITION [ReceiverPos]
           if (~isempty(ReceiverPos) | ~isnan(ReceiverPos))
             
             %***ROUND TO THE NEAREST THOUSAND 
             t1=roundn(ReceiverPos,4) ;
               
             if size(ReceiverPos,2)==1 
               if all(t1(:,1))==0                    
               h=ReceiverPos;%Assigning ReceiverPos to h(Elipsoidal Height)
               end  
               
            else 
                Rpos=ReceiverPos;%Receiver position(XYZ/LAT LONG h) 
                        
            end   %//if size(ReceiverPos,2)==1
                            
           else
               receiverpos_empty=1; %flag to indicate that receiverpos is empty([])
               setappdata(0,'receiverpos_empty',receiverpos_empty)
               
           end% if (~isempty(ReceiverPos) | ~isnan(ReceiverPos))               
                                                                 
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
               UTCtime = [];
               time_empty=1;
               setappdata(0,'time_empty',time_empty) 
                             
          end %if any(~isempty(Time(:,:)))| any(~isnan(Time(:,:)))
        
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
             %ISSUE ERROR MESSAGE for insufficient INPUTs
              beep%Give a beep sound 
              errmsg{1}=sprintf('Insuficient Data Input / Wrong Data Input format .');
              errmsg{2}='';
              errmsg{3}='Please Check file / Data format & Try Again.';
              errordlg(errmsg,'Tropospheric Delay Input Error','modal')  
              
              %Return empty ([]) outputs
              STD=[]; SHD=[]; SWD=[]; ZTD=[]; ZHD=[]; ZWD=[];
              
              return                          
end %switch nargin
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-                                   

%CONVERT TEMPERATURE(T) IN Celcius(°C) to Kelvin(K) 
 T = T + 273.15; %Temperature in kelvin(K)
            
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-   
   
%************TROPOHERIC DELAY MODELING USING VARIOUS MODELS
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-    
%ASSIGNMENT
grid_dry = grid_dryModel;
grid_wet = grid_wetModel;
gridRES_dry = gridRES_dryModel;
gridRES_wet = gridRES_wetModel;

Timevar = [Timevar_dry Timevar_wet];

%(1)*********COMPUTE ZENITH TROPOSPHERIC DELAY
if (strncmpi(dryModel,'no model',8) & strncmpi(wetModel,'no model',8))~=1
    
%Call the "ZenithTropDelay.m" subroutine
[ZTD, ZHD, ZWD] =ZenithTropDelay(UTCtime,T,P,es,Tm,lambda,dryModel,wetModel,...
                                 grid_dry,gridRES_dry,grid_wet,gridRES_wet,...
                                 Timevar_dry,Timevar_wet,GTropCoeff,VMFgrids) ;   

                                                                                                               
%(2)*********COMPUTE SLANT TROPOSPHERIC DELAY
if any([~exist('satpos_empty','var'),exist('satEL_deg','var')])
    
   %(2.1)****COMPUTE MAPPING FUNCTIONs(MF)
   
   %Call the "MF.m" function    
   [MFh, MFw] = MF(UTCtime,latD,lonD,h,satEL_deg,dryModel,wetModel,T,P,es,...
                   grid_dry,gridRES_dry,grid_wet,gridRES_wet,Timevar,VMFgrids);                         
                            
   
   %Call the "SlantTropDelay.m" subroutine
   [STD,SHD,SWD] = SlantTropDelay(ZHD,ZWD,MFh,MFw);
else
     STD=zeros(size(ZHD));
     SHD=zeros(size(ZHD));
     SWD=zeros(size(ZHD));     
end

else
    %Create zero matrix
     ZTD=0;
     ZHD=0;
     ZWD=0;    
     STD=0;
     SHD=0;
     SWD=0;        
end
%***********************END OF TropModels.m ***********************    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-                                             
                                             
%A.******SUBROUTINE TO COMPUTE ZENITH TROPOSPHERIC DELAY
function [ZTD, ZHD, ZWD] =ZenithTropDelay(UTCtime,T,P,es,Tm,lambda,dryModel,wetModel,...
                                          grid_dry,gridRES_dry,grid_wet,gridRES_wet,...
                                          Timevar_dry,Timevar_wet,GTropCoeff,VMFgrids)

%DESCRIPTION:
%            "ZenithTropDelay" Computes Zenith Tropospheric Delay Using ...* 
%             The Various models Specify by the User.                      *
%******SYNTAX:
%             [ZTD, ZHD, ZWD] = ZenithTropDelay(UTCtime,T,P,es,Tm,lambda,  *
%                                                        dryModel,wetModel)*
%******INPUT:                                                              *
%      UTCtime =  UTC time in [Year,Month,Day,Hour,Minute,Seconds]         * 
%            P = atmospheric pressure [mbar/hPa]                           *
%            T = atmospheric temperature in degrees Kelvin [k]             *
%           es = partial pressure of water vapor [mbar/hPa]                *
%           Tm = mean temperature of the water vapor                       *
%     dryModel = Hydrostatic/Dry Tropospheric Delay Model                  *
%     wetModel = Non-Hydrostatic/Wet Tropospheric Delay Model              *
%       lambda = water vapour `lapse rate'(dimensionless)]                 *
%****OUTPUT:
%           ZHD = Zenith Hydrostatic Delay                                 *
%           ZWD = Zenith Wet Delay                                         *
%           ZTD = Zenith Total Delay                                       *

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%*****RETRIEVE STORED DATA
latD = getappdata(0,'latD'); %Get station latitude Coords [degrees]
lonD = getappdata(0,'lonD'); %Get station latitude Coords [degrees]
latR = getappdata(0,'latR'); %Get station latitude Coords [radians]
lonR = getappdata(0,'lonR'); %Get station latitude Coords [radians]
   h = getappdata(0,'ht'); %Get station height [meters]
     
time_empty=getappdata(0,'time_empty');%Flag to indicate time is empty
receiverpos_empty=getappdata(0,'receiverpos_empty');

%************************COMPUTE ZENITH TROPOSPHERIC DELAYs

%*************COMPUTE ORTHOMETRIC HEIGHT(H)
%--------------------------------------------------------------------------  
%NOTE:
%SAASTAMOINEN,DAVIS ET AL,MARINI, Askne & Nordius,UNB3m,EGNOS, & MOPS MODELS 
%REQUIRE +VE ORTHOMETRIC HEIGHT---------------------------------------------

%COMPUTE CORRECTION FACTOR FOR LOCAL GRAVITY ACCELERATION 

if (~isempty(latR) && ~isempty(h))
    
  %COMPUTE UNDULATION & ORTHOMETRIC HEIGHT 
  [~,horth] = GETundulation(latD,lonD,h);  
  
  dgref=1-0.00266.*cos(2.*latR)-0.00028e-3.*horth;

else
    dgref=1;
    
end %//if (~isempty(latRAD) && ~isempty(h))

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%**************COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)***********************
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%*********FIND OUT THE TYPE OF MODEL PROVIDED

%--------------------------------------------------------------------------
%1.*************CONVENTIONAL DRY TROPOSPHERIC MODELs
%--------------------------------------------------------------------------
  %*******FOR SAASTAMOINEN DRY MODEL
if any([strncmpi(dryModel,'Saastamoinen',12),strncmpi(dryModel,'Saas',4),...        
        strncmpi(dryModel,'Saas(Refined)',14),strncmpi(dryModel,'Saastamoinen(Refined)',21)])
                  
   %*****COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
   if size(P,2) > 1
    
      for i = 1:size(P,2)
      
          ZHD(:,i)=(0.002277 .* P(:,i))./dgref;  %Zenith Hydrostatic or dry delay [m]    
      end  
  
   else  
       ZHD =(0.002277 .* P)./dgref;   %Zenith Hydrostatic or dry delay [m]    
   end 
  
 %*******FOR HOPFIELD DRY MODEL  
elseif any([strncmpi(dryModel,'Hopfield',8),strncmpi(dryModel,'Hop',3)])
       
      %***USING SMITH & WEINTRAUB REFRACTIVITY MODEL (1953)
       k1 = 77.6 ;%Refractivity constant(K/mbar or K/hpa) 
      
       if any([size(P,2) > 1, size(T,2) > 1])
          
          for i = 1 : size(P,2) 
              
              N_dry0 = k1 .* (P(:,i) ./ T(:,i)) ; %Dry refractivity along signal path @ MSL      
        
              %COMPUTE MEAN DRY TROPOSPHERE HEIGHT(Hd)
              Hd = 40136 + 148.72 .* (T(:,i) - 273.15);%[meters]
      
              %*****COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
              ZHD(:,i) = (1.0e-6 / 5) .* N_dry0 .* Hd;%[m] Zenith hydrostatic delay
       
          end
          
       else
           %***USING SMITH & WEINTRAUB REFRACTIVITY MODEL (1953)
           N_dry0 = k1 .* (P ./ T) ; %Dry refractivity along signal path @ MSL      
        
           %COMPUTE MEAN DRY TROPOSPHERE HEIGHT(Hd)
           Hd = 40136 + 148.72 .* (T - 273.15);%[meters]
      
           %*****COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
           ZHD =(1.0e-6 / 5) .* N_dry0 .* Hd;%[m] Zenith hydrostatic delay
           
       end %//if any([size(P,2) > 1, size(T,2) > 1])
  
      %*****FOR DAVIS ET AL DRY MODEL
elseif any([strncmpi(dryModel,'Davis et al',11),strncmpi(dryModel,'Davis',5),strncmpi(dryModel,'Davis etal',10)])
    
      %*****COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
      if size(P,2) > 1
    
        for i = 1 : size(P,2)
      
            ZHD(:,i) = (0.0022768 .* P(:,i))./dgref; %Zenith Hydrostatic or dry delay [m]    
        end   
  
      else   
          ZHD = (0.0022768 .* P)./dgref;  %Zenith Hydrostatic or dry delay [m]    
      end  
      
      %*******FOR BLACK DRY MODEL
elseif  strncmpi(dryModel,'Black',5)
    
      %*****COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
      if any([size(P,2) > 1,size(T,2) > 1])
          
         for i = 1 : size(P,2) 
             
             try
                ZHD(:,i) = (2.343*1.0e-3) .* (P(:,i)) .* ((T(:,i)-4.12)./T(:,i));%[m] 
             catch   
                  %COMPUTE HEIGHT DRY TROPOSPHERE(Hd) 
                  Hd = 148.98.*(T(:,i)-4.12);%[m] Height of Dry Troposphere above the earth
                 
                  ZHD(:,i) = 77.6*(1.0e-6/5).*(P(:,i)./T(:,i)).*Hd;%[m]
             end 
             
         end
         
      else
           %*****COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
           try
              ZHD = (2.343*1.0e-3) .* (P) .* ((T-4.12)./T);%[m] 
           catch  
                %COMPUTE HEIGHT DRY TROPOSPHERE(Hd) 
                Hd = 148.98.*(T-4.12);%[m] Height of Dry Troposphere above the earth
                
               ZHD = 77.6*(1.0e-6/5).*(P./T).*Hd;%[m]
           end 
           
      end %//if any([size(P,2) > 1, size(T,2) > 1])
      
      
      %FOR Modified Hopfield Model by Goad & Goodman(1974
elseif  any([strncmpi(dryModel,'Modified Hopfield(Goads & Goodman)',34),...
             strncmpi(dryModel,'Goad & Goodman',14),strncmpi(dryModel,'GG',2),...
             strncmpi(dryModel,'Modified Hop',12),strncmpi(dryModel,'Modified Hopfield',17)])
    
      %**COMPUTE DRY TROPO REFRACTIVITY(N) @ MSL
      %***USE ESSEN & FROOME REFRACTIVITY MODEL
      k1 = 77.624; %Refractivity constant(K/mbar or K/hpa)  
      
      if any([size(P,2) > 1,size(T,2) > 1])
          
        for i = 1 : size(P,2) 
            
            N_dry0 = k1 .* P(:,i) ./ T(:,i) ; %Dry refractivity along signal path       
            
            %******COMPUTE TROPOPAUSE HEIGHT FOR DRY DELAY
            try
               %Using Saastamoinen range Correction model,compute....
               %DRY "TOP" OF TROPOSPHERE
               hDry = (5*(0.002277).*P(:,i))./(N_dry0.*1.0e-6);%[meters]
            
            catch   
                 %Using Hopfield range Correction model,compute....
                 %DRY "TOP" OF TROPOSPHERE
                 hDry = 40136 + 148.72 .* (T(:,i) - 273.15);% [meters]
           
            end  %\\try

           %*****COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
           ZHD(:,i) = (1.0e-6 / 5) .* N_dry0 .* hDry;%[m]
           
        end
        
      else
          %******COMPUTE TROPOPAUSE HEIGHT FOR DRY DELAY
          N_dry0 = k1 .* P ./ T ; %Dry refractivity along signal path       
            
          %******COMPUTE TROPOPAUSE HEIGHT FOR DRY DELAY
          try
            %Using Saastamoinen range Correction model,compute....
            %DRY "TOP" OF TROPOSPHERE
            hDry = (5*(0.002277).*P)./(N_dry0.*1.0e-6);%[meters]
            
          catch   
               %Using Hopfield range Correction model,compute....
               %DRY "TOP" OF TROPOSPHERE
               hDry = 40136 + 148.72 .* (T - 273.16);% [meters]
           
          end   %\\try

          %*****COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
          ZHD = (1.0e-6 / 5) .* N_dry0 .* hDry;%[m] 
          
      end %//if any([size(P,2) > 1, size(T,2) > 1])
      
      
  %FOR Aske and Nordius DRY MODEL
elseif  any([strncmpi(dryModel,'Askne & Nordius',15),strncmpi(dryModel,'Askne and Nordius',17)])
      
      %*****DEFINE CONSTANTs
      k1 = 77.604; 				   % K/hPa      
      Md = 28.9644; %Molar Mass of Dry Air  in [g/mol]     
    dMtr = Md*10^-3;%molar mass of dry air in kg/mol     
       R = 8.3143;%universal gas constant in J/K/mol          
      Rd = R/dMtr ; %specific gas constant for dry consituents 
    
      %ACCELERATION GRAVITY AT THE ATMOSPHERIC COLUMN
      gm = 9.784 .* dgref;%[m/s^2]
      
      if size(P,2) > 1
    
         for i = 1:size(P,2)
             
             %*****COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
             ZHD(:,i) = (1.0e-6.*k1.*Rd.*P(:,i))./gm;%[m] 
         end 
         
      else  
          %*****COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
          ZHD = (1.0e-6.*k1.*Rd.*P)./gm;%[m]  
          
      end %//if size(P,2) > 1

%--------------------------------------------------------------------------
%******************'BLIND' DRY TROPOSPHERIC MODELs 
%--------------------------------------------------------------------------

%(1)*************GLOBAL PRESSURE & TEMPERATURE (GPT) MODELS
elseif any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),...
            strncmpi(dryModel,'GPT2w (5° x 5°)',14),strncmpi(dryModel,'GPT3 (1° x 1°)',14),...
            strncmpi(dryModel,'GPT3 (5° x 5°)',14)])
       
       if all([isempty(receiverpos_empty),isempty(time_empty)])
           
          [ZHD,ZWD_gpt,~,T_GPT] = ZTD_GPT(UTCtime,latR,lonR,h,grid_dry,grid_wet,gridRES_dry,gridRES_wet,dryModel,wetModel,Timevar_dry,Timevar_wet); 
          
       end
       
       
  %********GET & SAVE VMF COEFFIENTS     
   if strncmpi(dryModel,'GPT2 (5° x 5°)',14)
       
       setappdata(0,'ah_vmf1',getappdata(0,'ah_gpt2_5'))
       setappdata(0,'aw_gpt2_vmf1',getappdata(0,'aw_gpt2_5'))
       setappdata(0,'Drymodel','GPT2')
       
   elseif strncmpi(dryModel,'GPT2w (1° x 1°)',14)
       
          setappdata(0,'ah_vmf1',getappdata(0,'ah_gpt2w_1'))
          setappdata(0,'aw_gpt2w_vmf1',getappdata(0,'aw_gpt2w_1'))
          setappdata(0,'Drymodel','GPT2w_1')
          
   elseif strncmpi(dryModel,'GPT2w (5° x 5°)',14)
       
          setappdata(0,'ah_vmf1',getappdata(0,'ah_gpt2w_5'))
          setappdata(0,'aw_gpt2w_vmf1',getappdata(0,'aw_gpt2w_5'))
          setappdata(0,'Drymodel','GPT2w_5')
          
   elseif strncmpi(dryModel,'GPT3 (1° x 1°)',14)
          
          setappdata(0,'ah_vmf3',getappdata(0,'ah_gpt3_1'))
          setappdata(0,'aw_gpt3_vmf3',getappdata(0,'aw_gpt3_1'))
          setappdata(0,'Drymodel','GPT3_1') 
       
       
   elseif strncmpi(dryModel,'GPT3 (5° x 5°)',14)
          
          setappdata(0,'ah_vmf3',getappdata(0,'ah_gpt3_5'))
          setappdata(0,'aw_gpt3_vmf3',getappdata(0,'aw_gpt3_5'))
          setappdata(0,'Drymodel','GPT3_5')
          
   end
          
 %(2)****************VIENNA MAPPING FUNCTION(VMF1/VMF3) GRID MODELS 
elseif any([strncmpi(dryModel,'VMF gridded ZHD',15),strfind(dryModel,'VMF gridded ZHD')])
       
       %Call the "SearchReadVMFgrid.m" function(Internal) 
       [ZHD,ZWD_vmf,~,ah,aw,VMF_model] = searchReadVMFgrids(UTCtime,latD,lonD,h,VMFgrids);
       
       %SAVE VMF COEFFICIENTS
       if strcmpi(VMF_model,'VMF1')
          setappdata(0,'ah_vmfg1',ah)
          setappdata(0,'aw_vmfg1',aw)
          
       elseif strcmpi(VMF_model,'VMF3')
              setappdata(0,'ah_vmfg3',ah)
              setappdata(0,'aw_vmfg3',aw)
       end
       
      setappdata(0,'DryModel',VMF_model) 
          
%(3)****************UNB3m,EGNOS & MOPS MODELS 
elseif any([strncmpi(dryModel,'UNB3m',5),strncmpi(dryModel,'EGNOS',5),...
            strncmpi(dryModel,'MOPS',4)])
        
      if all([isempty(receiverpos_empty),isempty(time_empty)])
          
         [ZHD,ZWD_uem,~,T_uem]=ZTD_UNB3m_EGNOS_MOPS(UTCtime,latD,horth,dryModel,wetModel) ;   
         
      end
      
%(4)********************************GTrop MODEL
elseif any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),strfind(dryModel,'GTrop')])
    
       [ZHD,ZWD_GTrop] = getZTD_GTrop(UTCtime,latD,lonD,h,GTropCoeff);     
       
end %//if any([strncmpi(dryModel,'Saastamoinen',12),...

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%*****************COMPUTE ZENITH WET DELAY(ZHD)***********************
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%*********FIND OUT THE TYPE OF MODEL PROVIDED

%--------------------------------------------------------------------------
%********************CONVENTIONAL WET TROPOSPHERIC MODELs
%--------------------------------------------------------------------------

%(1)*********FOR SAASTAMOINEN WET MODEL
if any([strncmpi(wetModel,'Saastamoinen',12),strncmpi(wetModel,'Saas',4)])
    
  %*****COMPUTE ZENITH WET DELAY(ZWD)
  if any([size(T,2) > 1,size(es,2) > 1])
          
     for i = 1 : size(T,2) 
         
         ZWD(:,i) = 0.002277.*((1255 ./ T(:,i)) + 0.05).* es(:,i);% Zenith wet delay[m] 
     end
     
  else
      %*****COMPUTE ZENITH WET DELAY(ZWD)
      ZWD = 0.002277.*((1255 ./ T) + 0.05).* es;% Zenith wet delay[m]
      
  end
  
  %FOR HOPFIELD WET MODEL  
elseif  any([strncmpi(wetModel,'Hop',3),strncmpi(wetModel,'Hopfield',8)])
       
       %***USING SMITH & WEINTRAUB REFRACTIVITY MODEL (1953)       
       k3=3.73e5;%Refractivity constant(K^2/mbar or K^2/hpa)
       
       %COMPUTE MEAN WET TROPOSPHERE HEIGHT(Hw)
       Hw = 11000;%[meters]
       
       if any([size(T,2) > 1,size(es,2) > 1])
          
          for i = 1 : size(T,2) 
         
              N_wet0 = k3 .* (es(:,i) ./ T(:,i).^2);%Wet refractivity along signal path @ MSL 
       
             %*****COMPUTE ZENITH WET DELAY(ZWD)
             ZWD(:,i) = (1.0e-6 / 5) .* N_wet0 .* Hw;%[m]
          end
       
       else
           %COMPUTE REFRACTIVITY
           N_wet0 = k3 .* (es ./ T.^2);%Wet refractivity along signal path @ MSL 
       
          %*****COMPUTE ZENITH WET DELAY(ZWD)
          ZWD = (1.0e-6 / 5) .* N_wet0 .* Hw;%[m]
       
       end 
       
      %FOR BLACK WET MODEL
elseif strcmpi(wetModel,'Black')
          
      if ~isempty(latR) %Chech if latRAD is not empty        
        
        %*****INITIALIZING ZWD OUTPUT VARIABLE
        ZWD=zeros(size(latR,1),1);%Assign zeros of nxm to ZWD
      
        %******CHECK WHICH REGION Latitude FALLs
        for i=1:size(latD,1)   
            
           %Regions include:
           %1.Tropics/Low-Latitudes. lat range:(b/n (23.5°)N and  (23.5°)S 
           %2.Temperate zones/Mid-latitudes. lat range:From 23.5N to 66.5N and b/n 23.5S and 66.5S
           %3.Polar/High Latitudes. lat range:From 66.5N to the North Pole we have the Arctic, 
           %                                 and from 66.5S to the South Pole, the Antarctic.

           %*******GET REGIONAL EMPIRICAL CONSTANTs
           if abs(latD(i))<23.5  %for Tropics/Low-Latitude regions 
               
              ZWD(i,1)=0.28; %[m] 
        
           elseif  abs(latD(i)) >= 23.5 && abs(latD(i)) <= 66.5 %for Mid-latitude region
                 
                   ZWD(i,1)=0.20 ;%[m] 
      
           elseif  abs(latD(i))>66.5 %for Polar regions
               
                   ZWD(i,1)=0.05; %[m]   
           else          
                ZWD(i,1)=0.28; %[m]
    
           end     %if ~isempty(latRAD)
        
        end   %\\for i=1:size(lat,1)
        
      else 
           %COMPUTE HEIGHT OF WET TROPOSPHERE      
           Hw = 13000;%[meters] Height of wet Troposphere above the earth
           
           if any([size(T,2) > 1,size(es,2) > 1])
          
              for i = 1 : size(T,2) 
                  %*****COMPUTE ZENITH WET DELAY(ZWD)
                  ZWD(:,i) = 3.73e5*(1e-6/5).*(es(:,i)./T(:,i).^2).*Hw;%[m]
              end
              
           else
               %*****COMPUTE ZENITH WET DELAY(ZWD)
               ZWD = 3.73e5*(1e-6/5).*(es./T.^2).*Hw;%[m]
               
           end
           
      end  %\\if ~isempty(latR)     
       
     %FOR Modified Hopfield WET Model by Goad & Goodman(1974
elseif any([strncmpi(wetModel,'Modified Hopfield(Goads & Goodman)',34),...
            strncmpi(wetModel,'Goad & Goodman',14),strncmpi(wetModel,'GG',2),...
            strncmpi(wetModel,'Modified Hop',12),strncmpi(wetModel,'Modified Hopfield',17)])
    
      %**COMPUTE WET TROPO REFRACTIVITY(N) @ MSL
      %***USE ESSEN & FROOME REFRACTIVITY MODEL      
      k2 = -12.92; %Refractivity constant(K/mbar or K/hpa)
      k3=371900;%Refractivity constant(K^2/mbar or K^2/hpa) 
      
      if any([size(T,2) > 1,size(es,2) > 1])
          
         for i = 1 : size(T,2) 
            
             N_wet0 = k2.*(es(:,i)./ T(:,i))+k3.*( es(:,i) ./ T(:,i).^2);%Wet refractivity along signal path  
    
             %******COMPUTE TROPOPAUSE HEIGHT FOR WET DELAY
             try
                %Using Saastamoinen range Correction model,compute....         
                %*****WET "TOP" OF TROPOSPHERE
                hWet = (5*(0.002277)./(N_wet0.*1.0e-6)).*((1255./T(:,i))+0.05).*es(:,i);%[meters]
   
             catch   
                  %Using Hopfield range Correction model,compute....          
                  %*****WET "TOP" OF TROPOSPHERE
                  hWet = 11000;% [meters]
    
             end  %\\try

             %*****COMPUTE ZENITH WET DELAY(ZWD)
             ZWD(:,i) = (1.0e-6 / 5) .* N_wet0 .* hWet;%[m]
             
         end
         
      else
           %USE ESSEN & FROOME REFRACTIVITY MODEL COMPUTE WET TROPO REFRACTIVITY(N)
           N_wet0 = k2.*(es./ T)+k3.*( es ./ T.^2);%Wet refractivity along signal path  
    
           %******COMPUTE TROPOPAUSE HEIGHT FOR WET DELAY
           try
              %Using Saastamoinen range Correction model,compute....         
              %*****WET "TOP" OF TROPOSPHERE
              hWet = (5*(0.002277)./(N_wet0.*1.0e-6)).*((1255./T)+0.05).*es;%[meters]
   
           catch   
                %Using Hopfield range Correction model,compute....          
                %*****WET "TOP" OF TROPOSPHERE
                hWet = 11000;% [meters]
    
           end  %\\try

           %*****COMPUTE ZENITH WET DELAY(ZWD)
           ZWD = (1.0e-6 / 5) .* N_wet0 .* hWet;%[m]
           
      end
      
   %FOR CHAO WET MODEL
elseif strncmpi(wetModel,'Chao',4)
    
       if any([size(T,2) > 1,size(es,2) > 1])
          
          for i = 1 : size(T,2)
              
              %COMPUTE alpha(Temperature lapse rate)     
              alpha = abs(-5.93-0.0359.*(T(:,i)-273.15));%NB: T is in degrees celsius

              %CONVERT es in mbar to N/m^2
              es = es(:,i).*1.0e2; %es converted to Nm-2

              %*****COMPUTE ZENITH WET DELAY(ZWD)
              ZWD(:,i) = (1.63.*(es(:,i).^1.23./T(:,i).^2))+(2.05.*alpha.*(es(:,i).^1.46./T(:,i).^3));%[m]
          end
      
       else
           %COMPUTE alpha(Temperature lapse rate)     
           alpha = abs(-5.93-0.0359.*(T-273.15));%NB: T is in degrees celsius
          
           %CONVERT es in mbar to N/m^2
           es = es.*1.0e2; %es converted to Nm-2

           %*****COMPUTE ZENITH WET DELAY(ZWD)
           ZWD = (1.63.*(es.^1.23./T.^2))+(2.05.*alpha.*(es.^1.46./T.^3));%[m]
      
       end
      
      %FOR IFADIS WET MODEL
elseif strncmpi(wetModel,'Ifadis',6) 
    
       if any([size(P,2) > 1,size(T,2) > 1,size(es,2) > 1])
          
          for i = 1 : size(P,2)
              
              %*****COMPUTE ZENITH WET DELAY(ZWD)
              ZWD(:,i) = 0.00554-(0.880e-4.*(P(:,i)-1000))+(0.272e-4.*es(:,i))+(2.771.*(es(:,i)./T(:,i)));%[m]
          end
          
       else
           %*****COMPUTE ZENITH WET DELAY(ZWD)
          ZWD = 0.00554-(0.880e-4.*(P-1000))+(0.272e-4.*es)+(2.771.*(es./T));%[m]
          
       end
      
      %FOR CALLAHAN WET MODEL
elseif  strncmpi(wetModel,'Callahan',8) 
    
      if any([size(T,2) > 1,size(es,2) > 1])
          
          for i = 1 : size(T,2)
              
              %*****COMPUTE ZENITH WET DELAY(ZWD)      
             ZWD = 1035.*(es(:,i)./T(:,i).^2);%[m]
          end
          
      else 
          %*****COMPUTE ZENITH WET DELAY(ZWD)      
          ZWD = 1035.*(es./T.^2);%[m]
          
      end
      
      %FOR BERMAN 70 WET MODEL
elseif  strncmpi(wetModel,'Berman 70',9)    
      
      %**Define Hygrometric constants
      A = 17.1485;
      B = 4684.1;
      C = 38.45;
      
      if any([size(T,2) > 1,size(es,2) > 1])
          
         for i = 1 : size(T,2)
             
             %COMPUTE alpha(Temperature lapse rate)     
             alpha = abs(-5.93-0.0359.*(T(:,i)-273.15));%[K/Km] NB: T is in degrees celsius

             %*****COMPUTE ZENITH WET DELAY(ZWD)
             ZWD(:,i) = (373./(alpha.*(B-(A.*C)))).*(((1-(C./T(:,i))).^2).*es(:,i));%[m]
         end
         
      else
          %COMPUTE alpha(Temperature lapse rate)     
          alpha = abs(-5.93-0.0359.*(T-273.15));%[K/Km] NB: T is in degrees celsius

          %*****COMPUTE ZENITH WET DELAY(ZWD)
          ZWD = (373./(alpha.*(B-(A.*C)))).*(((1-(C./T)).^2).*es);%[m]
          
      end
     %FOR BERMAN 74 WET MODEL
elseif strncmpi(wetModel,'Berman 74',9) 
      
      %**Define Constant Term
       K = 0.3224;
       
       if any([size(T,2) > 1,size(es,2) > 1])
          
          for i = 1 : size(T,2)

              %*****COMPUTE ZENITH WET DELAY(ZWD)
              ZWD(:,i) = 10.946.*K.*(es(:,i)./T(:,i));%[m]
          end
          
       else
           %*****COMPUTE ZENITH WET DELAY(ZWD)
           ZWD = 10.946.*K.*(es./T);%[m]
        
       end 

     %FOR BERMAN TMOD WET MODEL
elseif strncmpi(wetModel,'Berman TMOD',11) 
      
      %**Define Constant Term
      K = 0.3281;

      if any([size(T,2) > 1,size(es,2) > 1])
          
          for i = 1 : size(T,2)
              
              %*****COMPUTE ZENITH WET DELAY(ZWD)
              ZWD(:,i) = 10.946.*K.*(es(:,i)./T(:,i));%[m] 
          end
          
      else 
          %*****COMPUTE ZENITH WET DELAY(ZWD)
          ZWD = 10.946.*K.*(es./T);%[m]  
      end
      
      %FOR Askne and Nordius WET MODEL
elseif  any([strncmpi(wetModel,'Askne & Nordius',15),strncmpi(wetModel,'Askne and Nordius',17)])
      
      %*****DEFINE CONSTANTs
      k1 = 77.604; 				   % K/hPa
      k2 = 64.79; 				   % K/hPa
      Mw = 18.0152;%Molar Mass of Water in [g/mol]
      Md = 28.9644; %Molar Mass of Dry Air  in [g/mol]
     k2p = k2 - k1*Mw/Md;          % K/hPa  [called k2 prime]
     k3  = 377600; 				   % KK/hPa

    %ACCELERATION GRAVITY AT THE ATMOSPHERIC COLUMN  
    gm = 9.784 .* dgref;%[ m/s^2 ]
        
     dMtr = 28.965*10^-3;%molar mass of dry air in kg/mol     
     R = 8.3143;%universal gas constant in J/K/mol          
     Rd = R/dMtr ; %specific gas constant for dry consituents 

     if any([size(Tm,2) > 1,size(es,2) > 1])
    
        for i = 1 : size(Tm,2)
      
            ZWD(:,i) = ((1e-6.*(k2p + k3./Tm(:,i)).*Rd)./(gm.*(lambda(:,i) + 1))).*es(:,i);%[m]
        end 
  
     else
         %*****COMPUTE ZENITH WET DELAY(ZWD)
         ZWD = ((1e-6.*(k2p + k3./Tm).*Rd)./(gm.*(lambda + 1))).*es;%[m]
     end 
     

     %MENDES WET MODEL
elseif any([strncmpi(wetModel,'Mendes',6),strncmpi(wetModel,'Mendes & Langley',16)])
    
      %*****COMPUTE ZENITH WET DELAY(ZWD)
      if size(es,2) > 1
    
        for i = 1 : size(es,2)
            
            ZWD(:,i) = 0.122+0.00943.*es(:,i);%[m]
        end
        
      else
          %*****COMPUTE ZENITH WET DELAY(ZWD)
          ZWD =0.122+0.00943.*es;%[m]
      end
      
%--------------------------------------------------------------------------
%******************'BLIND' DRY TROPOSPHERIC MODELs 
%--------------------------------------------------------------------------

%(1)*************GLOBAL PRESSURE & TEMPERATURE (GPT) MODELS
elseif any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),...
            strncmpi(wetModel,'GPT2w (5° x 5°)',14),strncmpi(wetModel,'GPT3 (1° x 1°)',14),...
            strncmpi(wetModel,'GPT3 (5° x 5°)',14)])
       
        if exist('ZWD_gpt','var')
            
             
            ZWD = ZWD_gpt;
              
        else 
             if all([isempty(receiverpos_empty),isempty(time_empty)])
           
                [~,ZWD,~,T_GPT] = ZTD_GPT(UTCtime,latR,lonR,h,grid_dry,grid_wet,gridRES_dry,gridRES_wet,dryModel,wetModel,Timevar_dry,Timevar_wet); 
       
             end  
                
        end %//if exist('ZWD_gpt','var')
        
        
        %********GET & SAVE VMF COEFFIENTS     
        if strncmpi(wetModel,'GPT2 (5° x 5°)',14)
       
           setappdata(0,'aw_vmf1',getappdata(0,'aw_gpt2_5'))
           setappdata(0,'ah_gpt2_vmf1',getappdata(0,'ah_gpt2_5'))
           setappdata(0,'Wetmodel','GPT2')
       
        elseif strncmpi(wetModel,'GPT2w (1° x 1°)',14)
       
              setappdata(0,'aw_vmf1',getappdata(0,'aw_gpt2w_1'))
              setappdata(0,'ah_gpt2w_vmf1',getappdata(0,'ah_gpt2w_1'))
              setappdata(0,'Wetmodel','GPT2w_1')
          
        elseif strncmpi(wetModel,'GPT2w (5° x 5°)',14)
       
               setappdata(0,'aw_vmf1',getappdata(0,'aw_gpt2w_5'))
               setappdata(0,'ah_gpt2w_vmf1',getappdata(0,'ah_gpt2w_5'))
               setappdata(0,'Wetmodel','GPT2w_5')
          
        elseif strncmpi(wetModel,'GPT3 (1° x 1°)',14)
          
               setappdata(0,'aw_vmf3',getappdata(0,'aw_gpt3_1'))
               setappdata(0,'ah_gpt3_vmf3',getappdata(0,'ah_gpt3_1'))
               setappdata(0,'Wetmodel','GPT3_1') 
              
        elseif strncmpi(wetModel,'GPT3 (5° x 5°)',14)
          
               setappdata(0,'aw_vmf3',getappdata(0,'aw_gpt3_5'))
               setappdata(0,'ah_gpt3_vmf3',getappdata(0,'ah_gpt3_5'))
               setappdata(0,'Wetmodel','GPT3_5')
          
        end 
           
                
%(2)****************VIENNA MAPPING FUNCTION(VMF1/VMF3) GRID MODELS 
elseif  any([strncmpi(wetModel,'VMF gridded ZWD',15),strfind(wetModel,'VMF gridded ZWD')])
    
        if exist('ZWD_vmf','var')
            
           ZWD = ZWD_vmf;
           
           
           %GET SAVED VMF COEFFICIENTS & VMF MODEL TYPE
           VMF_model = getappdata(0,'DryModel');
           
           %SAVE VMF COEFFICIENTS
           if strcmpi(VMF_model,'VMF1')
              ah = getappdata(0,'ah_vmfg1');
              aw = getappdata(0,'aw_vmfg1');
          
           elseif strcmpi(VMF_model,'VMF3')
                  ah = getappdata(0,'ah_vmfg3');
                  aw = getappdata(0,'aw_vmfg3');
           end 
    
           
        else 
             %Call the "SearchReadVMFgrid.m" function(internal) 
            [~,ZWD,~,ah,aw,VMF_model] = searchReadVMFgrids(UTCtime,latD,lonD,h,VMFgrids);

               
        end
        
     %SAVE VMF COEFFICIENTS
     if strcmpi(VMF_model,'VMF1')
        setappdata(0,'ah_vmfg1',ah)
        setappdata(0,'aw_vmfg1',aw)
          
     elseif  strcmpi(VMF_model,'VMF3')
             setappdata(0,'ah_vmfg3',ah)
             setappdata(0,'aw_vmfg3',aw)
     end 
       
     setappdata(0,'WetModel',VMF_model)    
          
%(3)****************UNB3m,EGNOS & MOPS MODELS 
elseif any([strncmpi(wetModel,'UNB3m',5),strncmpi(wetModel,'EGNOS',5),...
            strncmpi(wetModel,'MOPS',4)])
        
        if exist('ZWD_uem','var')
            
           ZWD = ZWD_uem;
              
        else 
             if all([isempty(receiverpos_empty),isempty(time_empty)])
          
                [~,ZWD,~,T_uem] = ZTD_UNB3m_EGNOS_MOPS(UTCtime,latD,horth,dryModel,wetModel) ;
             end 
             
        end %//if exist('ZWD_uem','var') 
        
%(4)********************************GTrop MODEL
elseif any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),strfind(wetModel,'GTrop')])
    
       if exist('ZWD_GTrop','var')
            
          ZWD = ZWD_GTrop;
              
       else  
           if all([isempty(receiverpos_empty),isempty(time_empty)])
    
              [~,ZWD] = getZTD_GTrop(UTCtime,latD,lonD,h,GTropCoeff); 
           end
           
       end %//if exist('ZWD_GTrop','var')
           
        
end %\\if any([strcmpi(wetModel,'Saas'),strcmpi(wetModel,'Saastamoinen')...
            
%****FINALLY, COMPUTE ZENITH TOTAL DELAY(ZTD)
if all([exist('ZHD','var'),exist('ZWD','var')])
    
   if any([size(ZHD,2) > 1,size(ZWD,2) > 1]) 
    
      for k = 1: size(ZHD,2)
       
          ZTD(:,k) = ZHD(:,k) + ZWD(:,k);%#ok<*AGROW> %Zenith Total Delay [m]
      end 
   
   else 
       ZTD = ZHD + ZWD;%Zenith Total Delay [m]
   end 
  
elseif all([exist('ZHD','var'),~exist('ZWD','var')])
    
      ZWD = zeros(size(ZHD)); 
      ZTD = zeros(size(ZHD));
      
elseif  all([~exist('ZHD','var'),exist('ZWD','var')])
      ZHD = zeros(size(ZWD)) ;
      ZTD = zeros(size(ZWD));
end

%SAVE TEMPERATURE OUTPUTS BY GPT,UNB3m,EGNOS & MOPS MODELS
%1.GPT MODELS(GPT2,GPT2w,GPT3)
if exist('T_GPT','var')  
  setappdata(0,'T_GPT',T_GPT)
end

%2.UNB3m,EGNOS & MOPS MODELS
if exist('T_uem','var')
   %CONVERT T_uem in Kelvin TO degree Celsius[°C]
   %NOTE:UNB3m,EGNOS,MOPS Models T OUTPUT IS IN degree Kelvin
   T_uem = T_uem - 273.15; 
   setappdata(0,'T_uem',T_uem) 
end

%***********************END OF ZenithTropDelay.m ***********************    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%B.SUBROUTINE TO COMPUTE SLANT TROPOSPHERIC DELAY (SATELLITE IN VIEW)
function [STD,SHD,SWD] =SlantTropDelay(ZHD,ZWD,MFh,MFw)
                                                                                              
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-                                             
%DESCRIPTION:
%            "SlantTropDelay" Computes Slant Tropospheric Delay Using ...  * 
%             The Various Tropospheric and Mapping Function models Specify *
%             by the User.
%USAGE:                                                                    *
%      [STD,SHD,SWD] =SlantTropDelay(ZHD,ZWD,MFh,MFw)                      *

%******INPUT:                                                              *                                                              *
%           ZHD:Zenith Hydrostatic Delay computed from various tropo models*
%           ZWD:Zenith Wet Delay computed from various tropo models        *   
%           MFh:Hydrostatic Mapping Function computed from Various MF model*   
%           MFw:Wet Mapping Function computed from Various MF models       * 
%****OUTPUT:                                                               *
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

%*********GET MET. PARAMETERS & ZENITH TROPO DELAYs FROM GPT MODELS
%         ----------------------------------------------------------
%C.SUBROUTINE TO COMPUTE MET. PARAMETERS & ZENITH TROPO DELAYs FROM GPT MODELS
function [ZHD,ZWD,ZTD,T] = ZTD_GPT(UTCtime,lat,lon,hell,grid_Dry,grid_Wet,grid_res_Dry,...
                                 grid_res_Wet,Drymodel,Wetmodel,Timevar_Dry,Timevar_Wet)

%********COMPUTE METEOROLOGICAL PARAMETERS USING GPT2w MODEL
%1.****NOW COMPUTE MET PARAMETERs 
%GET SIZE OF USER INPUT TIME & POSITION
nrow_time = size(UTCtime,1);
nrow_pos = size(lat,1);

%1.******UTCdate
Yr  = UTCtime(:,1);%get Year
Mn  = UTCtime(:,2);%get Month
Day = UTCtime(:,3);%get Day

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
       
      %*****INITIALIZE OUTPUTs 
      ZHD = zeros(nrow_time,1);
      [ZWD,ZTD,T] = deal(ZHD);%create copy of ZHD IN ZWD,ZTD,T
     
      for i = 1:nrow_time      
          
          %******CHECK IF BOTH DRY & WET TROPO MODELS ARE THE GPT TYPE
          if any([strncmpi(Drymodel,Wetmodel,14),strfind(Drymodel,Wetmodel)])
              
              %******CHECK GPT MODEL TYPE
             if any([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),strncmpi(Wetmodel,'GPT2 (5° x 5°)',14),...
                     strfind(Drymodel,'GPT2 (5° x 5°)'),strfind(Wetmodel,'GPT2 (5° x 5°)')])
            
                %Call the "GPT2_5_x_5.m" Function
                [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)] = GPT2_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry); 
             
          elseif any([strncmpi(Drymodel,'GPT2w (1° x 1°)',14),strncmpi(Drymodel,'GPT2w (5° x 5°)',14),...
                      strncmpi(Wetmodel,'GPT2w (1° x 1°)',14),strncmpi(Wetmodel,'GPT2w (5° x 5°)',14),...
                      strfind(Drymodel,'GPT2w (1° x 1°)'),strfind(Drymodel,'GPT2w (5° x 5°)'),...
                      strfind(Wetmodel,'GPT2w (1° x 1°)'),strfind(Wetmodel,'GPT2w (5° x 5°)')])
                  
                 if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                    %Call the "GPT2w_1_x_1.m" Function
                    [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)] = GPT2w_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);
                    
                 elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                        %Call the "GPT2w_5_x_5.m" Function
                        [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)] = GPT2w_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);
                         
                 end  
          
          elseif any([strncmpi(Drymodel,'GPT3 (1° x 1°)',14),strncmpi(Drymodel,'GPT3 (5° x 5°)',14),...
                      strncmpi(Wetmodel,'GPT3 (1° x 1°)',14),strncmpi(Wetmodel,'GPT3 (5° x 5°)',14),...
                      strfind(Drymodel,'GPT3 (1° x 1°)'),strfind(Drymodel,'GPT3 (5° x 5°)'),...
                      strfind(Wetmodel,'GPT3 (1° x 1°)'),strfind(Wetmodel,'GPT3 (5° x 5°)')])
                  
                 if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
              
                    %Call the the 1 degree grid version; "GPT3_1_x_1.m" 
                    [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)] = GPT3_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry); 
                     
                 elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
              
                        %Call the the 5 degree grid version;"GPT3_5_x_5.m" 
                        [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)] = GPT3_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);        
                 end
                 
             end
             
          elseif all([any([strncmpi(Drymodel,'GPT',3),strfind(Drymodel,'GPT')]),any([~strncmpi(Wetmodel,'GPT',3),~strfind(Wetmodel,'GPT')])])
                 
                 if any([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),strfind(Drymodel,'GPT2 (5° x 5°)')])
                 
                    %Call the "GPT2_5_x_5.m" Function
                    [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)] = GPT2_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry); 
             
                 
                 elseif any([strncmpi(Drymodel,'GPT2w (1° x 1°)',14),strncmpi(Drymodel,'GPT2w (5° x 5°)',14),...
                             strfind(Drymodel,'GPT2w (1° x 1°)'),strfind(Drymodel,'GPT2w (5° x 5°)')]) 
                         
                        if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                           %Call the "GPT2w_1_x_1.m" Function
                           [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)] = GPT2w_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);
                    
                        elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                               %Call the "GPT2w_5_x_5.m" Function
                              [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)] = GPT2w_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);
                         
                        end
                        
                 elseif any([strncmpi(Drymodel,'GPT3 (1° x 1°)',14),strncmpi(Drymodel,'GPT3 (5° x 5°)',14),...
                             strfind(Drymodel,'GPT3 (1° x 1°)'),strfind(Drymodel,'GPT3 (5° x 5°)')]) 
                         
                        if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                           %Call the "GPT2w_1_x_1.m" Function
                           [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)] = GPT3_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);
                    
                        elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                               %Call the "GPT2w_5_x_5.m" Function
                              [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)] = GPT3_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);
                         
                        end
                        
                 end %//if any([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),strfind(Drymodel,'GPT2 (5° x 5°)',14)])
            
          elseif all([any([~strncmpi(Drymodel,'GPT',3),~strfind(Drymodel,'GPT')]),any([strncmpi(Wetmodel,'GPT',3),strfind(Wetmodel,'GPT')])])
                 
                 if any([strncmpi(Wetmodel,'GPT2 (5° x 5°)',14),strfind(Wetmodel,'GPT2 (5° x 5°)')])
                 
                    %Call the "GPT2_5_x_5.m" Function
                    [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)] = GPT2_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Wet,Timevar_Wet); 
             
                 
                 elseif any([strncmpi(Wetmodel,'GPT2w (1° x 1°)',14),strncmpi(Wetmodel,'GPT2w (5° x 5°)',14),...
                             strfind(Wetmodel,'GPT2w (1° x 1°)'),strfind(Wetmodel,'GPT2w (5° x 5°)')]) 
                         
                        if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                           %Call the "GPT2w_1_x_1.m" Function
                           [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)] = GPT2w_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Wet,Timevar_Wet);
                    
                        elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                               %Call the "GPT2w_5_x_5.m" Function
                              [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)] = GPT2w_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Wet,Timevar_Wet);
                         
                        end
                        
                 elseif any([strncmpi(Wetmodel,'GPT3 (1° x 1°)',14),strncmpi(Wetmodel,'GPT3 (5° x 5°)',14),...
                             strfind(Wetmodel,'GPT3 (1° x 1°)'),strfind(Wetmodel,'GPT3 (5° x 5°)')]) 
                         
                        if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                           %Call the "GPT2w_1_x_1.m" Function
                           [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)] = GPT3_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Wet,Timevar_Wet);
                    
                        elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                               %Call the "GPT2w_5_x_5.m" Function
                              [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)] = GPT3_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Wet,Timevar_Wet);
                         
                        end
                        
                 end %//if any([strncmpi(Wetmodel,'GPT2 (5° x 5°)',14),strfind(Wetmodel,'GPT2 (5° x 5°)',14)])
            
          else  
              if all([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),any([strncmpi(Wetmodel,'GPT2w (1° x 1°)',14),...
                      strncmpi(Wetmodel,'GPT2w (5° x 5°)',14)])])
                 
                %*****COMPUTE ZHD
                %Call the "GPT2_5_x_5.m" Function
                [ZHD(i,1),~,~,T(i,1)] = GPT2_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);   
             
                %*****COMPUTE ZWD
                if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                    %Call the "GPT2w_1_x_1.m" Function
                    [~,ZWD(i,1),~,T(i,1)] = GPT2w_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Wet,Timevar_Wet);
                    
                 elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                        %Call the "GPT2w_5_x_5.m" Function
                        [~,ZWD(i,1),~,T(i,1)] = GPT2w_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Wet,Timevar_Wet);       
                end
                 
                %COMPUTE ZTD
                ZTD(i,1) = ZHD(i,1) + ZWD(i,1);
                
              elseif all([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),any([strncmpi(Wetmodel,'GPT3 (1° x 1°)',14),...
                          strncmpi(Wetmodel,'GPT3 (5° x 5°)',14)])])
                     
                     %*****COMPUTE ZHD
                     %Call the "GPT2_5_x_5.m" Function
                     [ZHD(i,1),~,~,T(i,1)] = GPT2_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);   
             
                     %*****COMPUTE ZWD
                     if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [~,ZWD(i,1),~,T(i,1)] = GPT3_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Wet,Timevar_Wet);
                    
                     elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [~,ZWD(i,1),~,T(i,1)] = GPT3_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Wet,Timevar_Wet);       
                     end 
                 
                     %****COMPUTE ZTD
                     ZTD(i,1) = ZHD(i,1) + ZWD(i,1);
                     
              elseif all([any([strncmpi(Drymodel,'GPT2w (1° x 1°)',14),strncmpi(Drymodel,'GPT2w (5° x 5°)',14)]),...
                               strncmpi(Wetmodel,'GPT2 (5° x 5°)',14)]) 
                     
                     %*****COMPUTE ZHD
                     if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [ZHD(i,1),~,~,T(i,1)] = GPT2w_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);
                    
                     elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [ZHD(i,1),~,~,T(i,1)] = GPT2w_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);       
                     end
                     
                     %*****COMPUTE ZWD
                     %Call the "GPT2_5_x_5.m" Function
                     [~,ZWD(i,1),~,T(i,1)] = GPT2_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Wet,Timevar_Wet);   
             
                     %****COMPUTE ZTD
                     ZTD(i,1) = ZHD(i,1) + ZWD(i,1);
                     
              elseif all([any([strncmpi(Drymodel,'GPT3 (1° x 1°)',14),strncmpi(Drymodel,'GPT3 (5° x 5°)',14)]),...
                               strncmpi(Wetmodel,'GPT2 (5° x 5°)',14)])
                           
                     %*****COMPUTE ZHD
                     if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [ZHD(i,1),~,~,T(i,1)] = GPT3_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);
                    
                     elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [ZHD(i,1),~,~,T(i,1)] = GPT3_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);       
                     end
                     
                     %*****COMPUTE ZWD
                     %Call the "GPT2_5_x_5.m" Function
                     [~,ZWD(i,1),~,T(i,1)] = GPT2_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Wet,Timevar_Wet);   
             
                     %****COMPUTE ZTD
                     ZTD(i,1) = ZHD(i,1) + ZWD(i,1); 
                     
              elseif all([any([strncmpi(Drymodel,'GPT2w (1° x 1°)',14),strncmpi(Drymodel,'GPT2w (5° x 5°)',14)]),...
                          any([strncmpi(Wetmodel,'GPT3 (1° x 1°)',14),strncmpi(Wetmodel,'GPT3 (5° x 5°)',14)])]) 
                     
                     %*****COMPUTE ZHD
                     if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [ZHD(i,1),~,~,T(i,1)] = GPT2w_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);
                    
                     elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [ZHD(i,1),~,~,T(i,1)] = GPT2w_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);       
                     end
                     
                     %*****COMPUTE ZWD       
                     if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [~,ZWD(i,1),~,T(i,1)] = GPT3_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Wet,Timevar_Wet);
                    
                     elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [~,ZWD(i,1),~,T(i,1)] = GPT3_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Wet,Timevar_Wet);       
                     end 
                 
                     %****COMPUTE ZTD
                     ZTD(i,1) = ZHD(i,1) + ZWD(i,1); 
                     
                     
              elseif all([any([strncmpi(Drymodel,'GPT3 (1° x 1°)',14),strncmpi(Drymodel,'GPT3 (5° x 5°)',14)]),...
                          any([strncmpi(Wetmodel,'GPT2w (1° x 1°)',14),strncmpi(Wetmodel,'GPT2w (5° x 5°)',14)])])      
                     
                     %*****COMPUTE ZHD
                     if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [ZHD(i,1),~,~,T(i,1)] = GPT3_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);
                    
                     elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [ZHD(i,1),~,~,T(i,1)] = GPT3_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Dry,Timevar_Dry);       
                     end
                     
                     %*****COMPUTE ZWD
                     if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [~,ZWD(i,1),~,T(i,1)] = GPT2w_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Wet,Timevar_Wet);
                    
                     elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [~,ZWD(i,1),~,T(i,1)] = GPT2w_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid_Wet,Timevar_Wet);       
                     end 
                 
                     %****COMPUTE ZTD
                     ZTD(i,1) = ZHD(i,1) + ZWD(i,1);
                                          
              end %//if any([strncmpi(Drymodel,Wetmodel,14),strfind(Drymodel,Wetmodel)])
          
          end %//if any([strncmpi(Drymodel,Wetmodel,14),strfind(Drymodel,Wetmodel)])
      
      end %//for i = 1:nrow_time
      
   else  
       %*****INITIALIZE OUTPUTs 
       ZHD = zeros(nrow_pos,nrow_time);
       [ZWD,ZTD,T] = deal(ZHD);%create copy of ZHD in ZWD,ZTD,T
    
     for i = 1:nrow_time %LOOP OVER TIME
         
        for j = 1:nrow_pos %LOOP OVER POSITIONS
            
            %******CHECK IF BOTH DRY & WET TROPO MODELS ARE THE GPT TYPE
            if any([strncmpi(Drymodel,Wetmodel,14),strfind(Drymodel,Wetmodel)])
              
               %******CHECK GPT MODEL TYPE
               if any([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),strncmpi(Wetmodel,'GPT2 (5° x 5°)',14),...
                       strfind(Drymodel,'GPT2 (5° x 5°)'),strfind(Wetmodel,'GPT2 (5° x 5°)')])
                
                  %Call the "GPT2_5_x_5.m" Function
                  [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
               
               elseif  any([strncmpi(Drymodel,'GPT2w (1° x 1°)',14),strncmpi(Drymodel,'GPT2w (5° x 5°)',14),...
                            strncmpi(Wetmodel,'GPT2w (1° x 1°)',14),strncmpi(Wetmodel,'GPT2w (5° x 5°)',14),...
                            strfind(Drymodel,'GPT2w (1° x 1°)'),strfind(Drymodel,'GPT2w (5° x 5°)'),...
                            strfind(Wetmodel,'GPT2w (1° x 1°)'),strfind(Wetmodel,'GPT2w (5° x 5°)')])
                  
                       if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                       
                          %Call the the 1 degree grid version; "GPT2w_1_x_1.m"  
                          [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
           
                       elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                       
                              %Call the the 5 degree grid version; "GPT2w_5_x_5.m"  
                              [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);    
                       end  
                   
               elseif any([strncmpi(Drymodel,'GPT3 (1° x 1°)',14),strncmpi(Drymodel,'GPT3 (5° x 5°)',14),...
                           strncmpi(Wetmodel,'GPT3 (1° x 1°)',14),strncmpi(Wetmodel,'GPT3 (5° x 5°)',14),...
                           strfind(Drymodel,'GPT3 (1° x 1°)'),strfind(Drymodel,'GPT3 (5° x 5°)'),...
                           strfind(Wetmodel,'GPT3 (1° x 1°)'),strfind(Wetmodel,'GPT3 (5° x 5°)')])
                  
                       %****CHECK GRID RESOLUTION
                       if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                       
                          %Call the "GPT2w_1_x_1.m" Function 
                          [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
           
                       elseif  grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                       
                               %Call the "GPT2w_5_x_5.m" Function 
                               [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);                  
                       end
                       
               end %//if any([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),strncmpi(Wetmodel,'GPT2 (5° x 5°)',14),...
                   %          strfind(Drymodel,'GPT2 (5° x 5°)',14),strfind(Wetmodel,'GPT2 (5° x 5°)',14)])
                       
            elseif all([any([strncmpi(Drymodel,'GPT',3),strfind(Drymodel,'GPT')]),any([~strncmpi(Wetmodel,'GPT',3),~strfind(Wetmodel,'GPT')])])
                 
                   if any([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),strfind(Drymodel,'GPT2 (5° x 5°)')])
                 
                      %Call the "GPT2_5_x_5.m" Function
                      [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry); 
             
                 
                   elseif any([strncmpi(Drymodel,'GPT2w (1° x 1°)',14),strncmpi(Drymodel,'GPT2w (5° x 5°)',14),...
                               strfind(Drymodel,'GPT2w (1° x 1°)'),strfind(Drymodel,'GPT2w (5° x 5°)')]) 
                         
                          if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                             %Call the "GPT2w_1_x_1.m" Function
                             [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                    
                          elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                                 %Call the "GPT2w_5_x_5.m" Function
                                 [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                         
                          end 
                        
                 elseif any([strncmpi(Drymodel,'GPT3 (1° x 1°)',14),strncmpi(Drymodel,'GPT3 (5° x 5°)',14),...
                             strfind(Drymodel,'GPT3 (1° x 1°)'),strfind(Drymodel,'GPT3 (5° x 5°)')]) 
                         
                        if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                           %Call the "GPT2w_1_x_1.m" Function
                           [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                    
                        elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                               %Call the "GPT2w_5_x_5.m" Function
                              [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                         
                        end
                        
                   end %//if any([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),strfind(Drymodel,'GPT2 (5° x 5°)',14)])
            
            elseif all([any([~strncmpi(Drymodel,'GPT',3),~strfind(Drymodel,'GPT')]),any([strncmpi(Wetmodel,'GPT',3),strfind(Wetmodel,'GPT')])])
                 
                   if any([strncmpi(Wetmodel,'GPT2 (5° x 5°)',14),strfind(Wetmodel,'GPT2 (5° x 5°)')])
                 
                      %Call the "GPT2_5_x_5.m" Function
                      [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet); 
             
                 
                   elseif any([strncmpi(Wetmodel,'GPT2w (1° x 1°)',14),strncmpi(Wetmodel,'GPT2w (5° x 5°)',14),...
                               strfind(Wetmodel,'GPT2w (1° x 1°)'),strfind(Wetmodel,'GPT2w (5° x 5°)')]) 
                         
                        if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                           %Call the "GPT2w_1_x_1.m" Function
                           [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                    
                        elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                               %Call the "GPT2w_5_x_5.m" Function
                              [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                         
                        end
                        
                 elseif any([strncmpi(Wetmodel,'GPT3 (1° x 1°)',14),strncmpi(Wetmodel,'GPT3 (5° x 5°)',14),...
                             strfind(Wetmodel,'GPT3 (1° x 1°)'),strfind(Wetmodel,'GPT3 (5° x 5°)')]) 
                         
                        if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                           %Call the "GPT2w_1_x_1.m" Function
                           [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                    
                        elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                               %Call the "GPT2w_5_x_5.m" Function
                              [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                         
                        end
                        
                   end %//if any([strncmpi(Wetmodel,'GPT2 (5° x 5°)',14),strfind(Wetmodel,'GPT2 (5° x 5°)',14)])          
                               
            else   
                 if all([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),any([strncmpi(Wetmodel,'GPT2w (1° x 1°)',14),...
                         strncmpi(Wetmodel,'GPT2w (5° x 5°)',14)])])
                 
                    %*****COMPUTE ZHD
                    %Call the "GPT2_5_x_5.m" Function
                    [ZHD(j,i),~,~,T(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);   
             
                %*****COMPUTE ZWD
                if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                    %Call the "GPT2w_1_x_1.m" Function
                    [~,ZWD(j,i),~,T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                    
                 elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                        %Call the "GPT2w_5_x_5.m" Function
                        [~,ZWD(j,i),~,T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);       
                end
                 
                %COMPUTE ZTD
                ZTD(j,i) = ZHD(j,i) + ZWD(j,i);
                
              elseif all([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),any([strncmpi(Wetmodel,'GPT3 (1° x 1°)',14),...
                          strncmpi(Wetmodel,'GPT3 (5° x 5°)',14)])])
                     
                     %*****COMPUTE ZHD
                     %Call the "GPT2_5_x_5.m" Function
                     [ZHD(j,i),~,~,T(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);   
             
                     %*****COMPUTE ZWD
                     if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [~,ZWD(j,i),~,T(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                    
                     elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [~,ZWD(j,i),~,T(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);       
                     end 
                 
                     %****COMPUTE ZTD
                     ZTD(j,i) = ZHD(j,i) + ZWD(j,i);
                     
              elseif all([any([strncmpi(Drymodel,'GPT2w (1° x 1°)',14),strncmpi(Drymodel,'GPT2w (5° x 5°)',14)]),...
                               strncmpi(Wetmodel,'GPT2 (5° x 5°)',14)]) 
                     
                     %*****COMPUTE ZHD
                     if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [ZHD(j,i),~,~,T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                    
                     elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [ZHD(j,i),~,~,T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);       
                     end
                     
                     %*****COMPUTE ZWD
                     %Call the "GPT2_5_x_5.m" Function
                     [~,ZWD(j,i),~,T(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);   
             
                     %****COMPUTE ZTD
                     ZTD(j,i) = ZHD(j,i) + ZWD(j,i);
                     
              elseif all([any([strncmpi(Drymodel,'GPT3 (1° x 1°)',14),strncmpi(Drymodel,'GPT3 (5° x 5°)',14)]),...
                               strncmpi(Wetmodel,'GPT2 (5° x 5°)',14)])
                           
                     %*****COMPUTE ZHD
                     if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [ZHD(j,i),~,~,T(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                    
                     elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [ZHD(j,i),~,~,T(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);       
                     end
                     
                     %*****COMPUTE ZWD
                     %Call the "GPT2_5_x_5.m" Function
                     [~,ZWD(j,i),~,T(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);   
             
                     %****COMPUTE ZTD
                     ZTD(j,i) = ZHD(j,i) + ZWD(j,i); 
                     
              elseif all([any([strncmpi(Drymodel,'GPT2w (1° x 1°)',14),strncmpi(Drymodel,'GPT2w (5° x 5°)',14)]),...
                          any([strncmpi(Wetmodel,'GPT3 (1° x 1°)',14),strncmpi(Wetmodel,'GPT3 (5° x 5°)',14)])]) 
                     
                     %*****COMPUTE ZHD
                     if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [ZHD(j,i),~,~,T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                    
                     elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [ZHD(j,i),~,~,T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);       
                     end
                     
                     %*****COMPUTE ZWD       
                     if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [~,ZWD(j,i),~,T(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                    
                     elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [~,ZWD(j,i),~,T(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);       
                     end 
                 
                     %****COMPUTE ZTD
                     ZTD(j,i) = ZHD(j,i) + ZWD(j,i); 
                     
                     
              elseif all([any([strncmpi(Drymodel,'GPT3 (1° x 1°)',14),strncmpi(Drymodel,'GPT3 (5° x 5°)',14)]),...
                          any([strncmpi(Wetmodel,'GPT2w (1° x 1°)',14),strncmpi(Wetmodel,'GPT2w (5° x 5°)',14)])])      
                     
                     %*****COMPUTE ZHD
                     if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [ZHD(j,i),~,~,T(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                    
                     elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [ZHD(j,i),~,~,T(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);       
                     end
                     
                     %*****COMPUTE ZWD
                     if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [~,ZWD(j,i),~,T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                    
                     elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [~,ZWD(j,i),~,T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);       
                     end 
                 
                     %****COMPUTE ZTD
                     ZTD(j,i) = ZHD(j,i) + ZWD(j,i);
                                          
                 end  %//if any([strncmpi(Drymodel,Wetmodel,14),strfind(Drymodel,Wetmodel)])
          
            end %//if any([strncmpi(Drymodel,Wetmodel,14),strfind(Drymodel,Wetmodel)])
                                 
        end %//for j = 1:nrow_pos
                
     end %//for i = 1:nrow_time
     
   end %//if isequal(Identical,1)
    
else %IF USER POSITION & OBSERVATION TIME ARE UNEQUAL
    %*****INITIALIZE OUTPUTs 
    ZHD = zeros(nrow_pos,nrow_time);
    [ZWD,ZTD,T] = deal(ZHD);%create copy of ZHD in ZWD,ZTD,T
    
    for i = 1:nrow_time %LOOP OVER TIME
        
        for j = 1:nrow_pos %LOOP OVER POSITIONS
            
           %******CHECK IF BOTH DRY & WET TROPO MODELS ARE THE GPT TYPE
            if any([strncmpi(Drymodel,Wetmodel,14),strfind(Drymodel,Wetmodel)])
              
               %******CHECK GPT MODEL TYPE
               if any([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),strncmpi(Wetmodel,'GPT2 (5° x 5°)',14),...
                       strfind(Drymodel,'GPT2 (5° x 5°)'),strfind(Wetmodel,'GPT2 (5° x 5°)')])
                
                  %Call the "GPT2_5_x_5.m" Function
                  [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
               
               elseif  any([strncmpi(Drymodel,'GPT2w (1° x 1°)',14),strncmpi(Drymodel,'GPT2w (5° x 5°)',14),...
                            strncmpi(Wetmodel,'GPT2w (1° x 1°)',14),strncmpi(Wetmodel,'GPT2w (5° x 5°)',14),...
                            strfind(Drymodel,'GPT2w (1° x 1°)'),strfind(Drymodel,'GPT2w (5° x 5°)'),...
                            strfind(Wetmodel,'GPT2w (1° x 1°)'),strfind(Wetmodel,'GPT2w (5° x 5°)')])
                  
                       if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                       
                          %Call the the 1 degree grid version; "GPT2w_1_x_1.m"  
                          [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
           
                       elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                       
                              %Call the the 5 degree grid version; "GPT2w_5_x_5.m"  
                              [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);    
                       end  
                   
               elseif any([strncmpi(Drymodel,'GPT3 (1° x 1°)',14),strncmpi(Drymodel,'GPT3 (5° x 5°)',14),...
                           strncmpi(Wetmodel,'GPT3 (1° x 1°)',14),strncmpi(Wetmodel,'GPT3 (5° x 5°)',14),...
                           strfind(Drymodel,'GPT3 (1° x 1°)'),strfind(Drymodel,'GPT3 (5° x 5°)'),...
                           strfind(Wetmodel,'GPT3 (1° x 1°)'),strfind(Wetmodel,'GPT3 (5° x 5°)')])
                  
                       %****CHECK GRID RESOLUTION
                       if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                       
                          %Call the "GPT2w_1_x_1.m" Function 
                          [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
           
                       elseif  grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                       
                               %Call the "GPT2w_5_x_5.m" Function 
                               [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);                  
                       end
                       
               end %//if any([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),strncmpi(Wetmodel,'GPT2 (5° x 5°)',14),...
                   %          strfind(Drymodel,'GPT2 (5° x 5°)',14),strfind(Wetmodel,'GPT2 (5° x 5°)',14)])
                       
            elseif all([any([strncmpi(Drymodel,'GPT',3),strfind(Drymodel,'GPT')]),any([~strncmpi(Wetmodel,'GPT',3),~strfind(Wetmodel,'GPT')])])
                 
                   if any([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),strfind(Drymodel,'GPT2 (5° x 5°)')])
                 
                      %Call the "GPT2_5_x_5.m" Function
                      [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry); 
             
                 
                   elseif any([strncmpi(Drymodel,'GPT2w (1° x 1°)',14),strncmpi(Drymodel,'GPT2w (5° x 5°)',14),...
                               strfind(Drymodel,'GPT2w (1° x 1°)'),strfind(Drymodel,'GPT2w (5° x 5°)')]) 
                         
                          if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                             %Call the "GPT2w_1_x_1.m" Function
                             [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                    
                          elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                                 %Call the "GPT2w_5_x_5.m" Function
                                 [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                         
                          end 
                        
                 elseif any([strncmpi(Drymodel,'GPT3 (1° x 1°)',14),strncmpi(Drymodel,'GPT3 (5° x 5°)',14),...
                             strfind(Drymodel,'GPT3 (1° x 1°)'),strfind(Drymodel,'GPT3 (5° x 5°)')]) 
                         
                        if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                           %Call the "GPT2w_1_x_1.m" Function
                           [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                    
                        elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                               %Call the "GPT2w_5_x_5.m" Function
                              [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                         
                        end
                        
                   end %//if any([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),strfind(Drymodel,'GPT2 (5° x 5°)',14)])
            
            elseif all([any([~strncmpi(Drymodel,'GPT',3),~strfind(Drymodel,'GPT')]),any([strncmpi(Wetmodel,'GPT',3),strfind(Wetmodel,'GPT')])])
                 
                   if any([strncmpi(Wetmodel,'GPT2 (5° x 5°)',14),strfind(Wetmodel,'GPT2 (5° x 5°)')])
                 
                      %Call the "GPT2_5_x_5.m" Function
                      [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet); 
             
                 
                   elseif any([strncmpi(Wetmodel,'GPT2w (1° x 1°)',14),strncmpi(Wetmodel,'GPT2w (5° x 5°)',14),...
                               strfind(Wetmodel,'GPT2w (1° x 1°)'),strfind(Wetmodel,'GPT2w (5° x 5°)')]) 
                         
                        if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                           %Call the "GPT2w_1_x_1.m" Function
                           [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                    
                        elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                               %Call the "GPT2w_5_x_5.m" Function
                              [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                         
                        end
                        
                 elseif any([strncmpi(Wetmodel,'GPT3 (1° x 1°)',14),strncmpi(Wetmodel,'GPT3 (5° x 5°)',14),...
                             strfind(Wetmodel,'GPT3 (1° x 1°)'),strfind(Wetmodel,'GPT3 (5° x 5°)')]) 
                         
                        if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                           %Call the "GPT2w_1_x_1.m" Function
                           [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                    
                        elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                               %Call the "GPT2w_5_x_5.m" Function
                              [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                         
                        end
                        
                   end %//if any([strncmpi(Wetmodel,'GPT2 (5° x 5°)',14),strfind(Wetmodel,'GPT2 (5° x 5°)',14)])          
                               
            else   
                 if all([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),any([strncmpi(Wetmodel,'GPT2w (1° x 1°)',14),...
                         strncmpi(Wetmodel,'GPT2w (5° x 5°)',14)])])
                 
                    %*****COMPUTE ZHD
                    %Call the "GPT2_5_x_5.m" Function
                    [ZHD(j,i),~,~,T(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);   
             
                %*****COMPUTE ZWD
                if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                    %Call the "GPT2w_1_x_1.m" Function
                    [~,ZWD(j,i),~,T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                    
                 elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                        %Call the "GPT2w_5_x_5.m" Function
                        [~,ZWD(j,i),~,T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);       
                end
                 
                %COMPUTE ZTD
                ZTD(j,i) = ZHD(j,i) + ZWD(j,i);
                
              elseif all([strncmpi(Drymodel,'GPT2 (5° x 5°)',14),any([strncmpi(Wetmodel,'GPT3 (1° x 1°)',14),...
                          strncmpi(Wetmodel,'GPT3 (5° x 5°)',14)])])
                     
                     %*****COMPUTE ZHD
                     %Call the "GPT2_5_x_5.m" Function
                     [ZHD(j,i),~,~,T(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);   
             
                     %*****COMPUTE ZWD
                     if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [~,ZWD(j,i),~,T(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                    
                     elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [~,ZWD(j,i),~,T(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);       
                     end 
                 
                     %****COMPUTE ZTD
                     ZTD(j,i) = ZHD(j,i) + ZWD(j,i);
                     
              elseif all([any([strncmpi(Drymodel,'GPT2w (1° x 1°)',14),strncmpi(Drymodel,'GPT2w (5° x 5°)',14)]),...
                               strncmpi(Wetmodel,'GPT2 (5° x 5°)',14)]) 
                     
                     %*****COMPUTE ZHD
                     if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [ZHD(j,i),~,~,T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                    
                     elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [ZHD(j,i),~,~,T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);       
                     end
                     
                     %*****COMPUTE ZWD
                     %Call the "GPT2_5_x_5.m" Function
                     [~,ZWD(j,i),~,T(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);   
             
                     %****COMPUTE ZTD
                     ZTD(j,i) = ZHD(j,i) + ZWD(j,i);
                     
              elseif all([any([strncmpi(Drymodel,'GPT3 (1° x 1°)',14),strncmpi(Drymodel,'GPT3 (5° x 5°)',14)]),...
                               strncmpi(Wetmodel,'GPT2 (5° x 5°)',14)])
                           
                     %*****COMPUTE ZHD
                     if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [ZHD(j,i),~,~,T(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                    
                     elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [ZHD(j,i),~,~,T(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);       
                     end
                     
                     %*****COMPUTE ZWD
                     %Call the "GPT2_5_x_5.m" Function
                     [~,ZWD(j,i),~,T(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);   
             
                     %****COMPUTE ZTD
                     ZTD(j,i) = ZHD(j,i) + ZWD(j,i); 
                     
              elseif all([any([strncmpi(Drymodel,'GPT2w (1° x 1°)',14),strncmpi(Drymodel,'GPT2w (5° x 5°)',14)]),...
                          any([strncmpi(Wetmodel,'GPT3 (1° x 1°)',14),strncmpi(Wetmodel,'GPT3 (5° x 5°)',14)])]) 
                     
                     %*****COMPUTE ZHD
                     if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [ZHD(j,i),~,~,T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                    
                     elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [ZHD(j,i),~,~,T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);       
                     end
                     
                     %*****COMPUTE ZWD       
                     if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [~,ZWD(j,i),~,T(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                    
                     elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [~,ZWD(j,i),~,T(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);       
                     end 
                 
                     %****COMPUTE ZTD
                     ZTD(j,i) = ZHD(j,i) + ZWD(j,i); 
                     
                     
              elseif all([any([strncmpi(Drymodel,'GPT3 (1° x 1°)',14),strncmpi(Drymodel,'GPT3 (5° x 5°)',14)]),...
                          any([strncmpi(Wetmodel,'GPT2w (1° x 1°)',14),strncmpi(Wetmodel,'GPT2w (5° x 5°)',14)])])      
                     
                     %*****COMPUTE ZHD
                     if grid_res_Dry == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [ZHD(j,i),~,~,T(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);
                    
                     elseif grid_res_Dry == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [ZHD(j,i),~,~,T(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Dry,Timevar_Dry);       
                     end
                     
                     %*****COMPUTE ZWD
                     if grid_res_Wet == 1 %if grid resolution is 1 x 1 degree
                     
                        %Call the "GPT2w_1_x_1.m" Function
                        [~,ZWD(j,i),~,T(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);
                    
                     elseif grid_res_Wet == 5 %if grid resolution is 5 x 5 degree
                     
                            %Call the "GPT2w_5_x_5.m" Function
                            [~,ZWD(j,i),~,T(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid_Wet,Timevar_Wet);       
                     end 
                 
                     %****COMPUTE ZTD
                     ZTD(j,i) = ZHD(j,i) + ZWD(j,i);
                                          
                 end   %//if any([strncmpi(Drymodel,Wetmodel,14),strfind(Drymodel,Wetmodel)])
          
            end  %//if any([strncmpi(Drymodel,Wetmodel,14),strfind(Drymodel,Wetmodel)])
                                 
        end  %//for j = 1:nrow_pos
                
    end  %//for i = 1:nrow_time
    
end  %//if isequal(nrow_time,nrow_pos)


%***********************END OF ZTD_GPT.m ***********************    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%C_1.SUBROUTINE TO COMPUTE MET. PARAMETERS & ZENITH TROPO DELAYs FROM GPT2w 5°x5° GRID MODEL
function [ZHD,ZWD,ZTD,T] = GPT2w_5_x_5(UTCtime,lat,lon,hell,grid,Timevar)

%**************************************************************************
%DESCRIPTION:
%            This subroutine determines pressure, temperature, temperature *
%            lapse rate, mean temperature of the water vapor, water vapor  * 
%            pressure, hydrostatic and wet mapping function coefficients   * 
%            ah and aw, water vapour decrease factor and geoid undulation  * 
%            for specific sites near the Earth surface. It is based on a   *
%            5° x 5° external grid file ('gpt2_5w.grd') with mean          *
%            values as well as sine and cosine amplitudes for the annual   * 
%            and semiannual variation of the coefficients.
%USAGE:
%      [ZHD,ZWD,ZTD] = GPT2w_5(UTCtime,lat,lon,hell,Timevar)               *
%INPUTs:                                                                   *
%1.   UTCtime:  UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *                                                                
%2.       lat:  Station geodetic latitude in radians [-pi/2:+pi/2] (n x 1) *
%3.       lon:  Station geodetic longitude in radians [-pi:pi] or [0:2pi]  *
%4.      hell:  ellipsoidal height in [meters] (n x 1)                     *
%5    Timevar:  case 1: no time variation but static quantities            *
%               case 0: with time variation (annual and semiannual terms)  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% OUTPUTs:
%1.     ZHD: Zenith Hydrostaic Tropospheric Delay in meters                *
%2.     ZWD: Zenith Wet Tropospheric Delay  in meters                      *
%3.    ZTD: Zenith Total Tropospheric Delay in meters                     *
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

%SAVE MAPPING FUNCTION COEFFICIENTS
setappdata(0,'ah_gpt2w_5',ah)
setappdata(0,'aw_gpt2w_5',aw)

%--------------------------------------------------------------------------
%(2)***COMPUTE ZENITH TROPOSPHERIC DELAYS USING EXTRACTED MET PARAMETERS
%--------------------------------------------------------------------------

%**************COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
%              =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%--------------------------------------------------------------------------
%NOTE:
%      DAVIS ET AL DRY MODEL REQUIRES +VE ORTOMETRIC HEIGHT TO PERFORM
%--------------------------------------------------------------------------
%COMPUTE CORRECTION FACTOR FOR LOCAL GRAVITY ACCELERATION    
dgref = 1-0.00266.*cos(2.*lat)-0.00028e-3.*(hell-undu);

 gm   = 9.784 .* dgref;%Acceleration Gravity at the Atmospheric column in m/s^2

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

%*****COMPUTE ZENITH WET DELAY(ZWD)    
ZWD = ((1e-6.*(k2p + k3./Tm).*Rd)./(gm.*(la + 1))).*e;%[m]
                
%****FINALLY, COMPUTE ZENITH TOTAL DELAY(ZTD)
ZTD = ZHD + ZWD;%Zenith Total Delay [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GPT2w_5 x 5.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%C_2.SUBROUTINE TO COMPUTE MET. PARAMETERS & ZENITH TROPO DELAYs FROM GPT2w 1°x1° GRID MODEL
function [ZHD,ZWD,ZTD,T] = GPT2w_1_x_1(UTCtime,lat,lon,hell,grid,Timevar)
       
%**************************************************************************
%DESCRIPTION:
%            This subroutine determines pressure, temperature, temperature *
%            lapse rate, mean temperature of the water vapor, water vapor  * 
%            pressure, hydrostatic and wet mapping function coefficients   * 
%            ah and aw, water vapour decrease factor and geoid undulation  * 
%            for specific sites near the Earth surface. It is based on a   *
%            1° x 1° external grid file ('gpt2_1w.grd') with mean          *
%            values as well as sine and cosine amplitudes for the annual   * 
%            and semiannual variation of the coefficients.
%USAGE:
%      [ZHD,ZWD,ZTD] = GPT2w_1_x_1(UTCtime,lat,lon,hell,Timevar)           *
%INPUTs:                                                                   *
%1.   UTCtime:  UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *                                                                
%2.       lat:  Station geodetic latitude in radians [-pi/2:+pi/2] (n x 1) *
%3.       lon:  Station geodetic longitude in radians [-pi:pi] or [0:2pi]  *
%4.      hell:  ellipsoidal height in [meters] (n x 1)                     *
%5    Timevar:  case 1: no time variation but static quantities            *
%               case 0: with time variation (annual and semiannual terms)  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% OUTPUTs:
%1.     ZHD: Zenith Hydrostaic Tropospheric Delay in meters                *
%2.     ZWD: Zenith Wet Tropospheric Delay  in meters                      *
%3.    ZTD: Zenith Total Tropospheric Delay in meters                     *
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

%SAVE MAPPING FUNCTION COEFFICIENTS
setappdata(0,'ah_gpt2w_1',ah)
setappdata(0,'aw_gpt2w_1',aw)

%--------------------------------------------------------------------------
%(2)***COMPUTE ZENITH TROPOSPHERIC DELAYS USING EXTRACTED MET PARAMETERS
%--------------------------------------------------------------------------
%NOTE:
%      DAVIS ET AL DRY MODEL REQUIRES +VE ORTOMETRIC HEIGHT TO PERFORM
%--------------------------------------------------------------------------

%COMPUTE CORRECTION FACTOR FOR LOCAL GRAVITY ACCELERATION    
dgref = 1-0.00266.*cos(2.*lat)-0.00028e-3.*(hell-undu);

 gm   = 9.784 .* dgref;%Acceleration Gravity at the Atmospheric column in m/s^2

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

%*****COMPUTE ZENITH WET DELAY(ZWD)    
ZWD = ((1e-6.*(k2p + k3./Tm).*Rd)./(gm.*(la + 1))).*e;%[m]
                
%****FINALLY, COMPUTE ZENITH TOTAL DELAY(ZTD)
ZTD = ZHD + ZWD;%Zenith Total Delay [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GPT2w_1_x_1.m  %%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%C_3.SUBROUTINE TO COMPUTE MET. PARAMETERS & ZENITH TROPO DELAYs FROM GPT3 5°x5° GRID MODEL
function [ZHD,ZWD,ZTD,T] = GPT3_5_x_5(UTCtime,lat,lon,h_ell,grid,Timevar)
         
%**************************************************************************
%DESCRIPTION:
%            This subroutine determines pressure, temperature, temperature *
%            lapse rate, mean temperature of the water vapor, water vapor  * 
%            pressure, hydrostatic and wet mapping function coefficients   * 
%            ah and aw, water vapour decrease factor, geoid undulation, and* 
%            empirical tropospheric gradients for specific sites near the  *
%            Earth surface. It is based on a  5° x 5° external grid        *
%            file ('gpt3_5.grd') with mean  values as well as sine and     * 
%            cosine amplitudes for the annual and semiannual variation of  *
%            the coefficients.
%USAGE:
%      [ZHD,ZWD,ZTD] = GPT3_5_x_5(UTCtime,lat,lon,hell,grid,Timevar)       *
%INPUTs:                                                                   *
%1.   UTCtime:  UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *                                                                
%2.       lat:  Station geodetic latitude in radians [-pi/2:+pi/2] (n x 1) *
%3.       lon:  Station geodetic longitude in radians [-pi:pi] or [0:2pi]  *
%4.      hell:  ellipsoidal height in [meters] (n x 1)                     *
%5.      grid: grid values in cells extracted from the grid file           *
%6    Timevar:  case 1: no time variation but static quantities            *
%               case 0: with time variation (annual and semiannual terms)  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% OUTPUTs:
%1.     ZHD: Zenith Hydrostaic Tropospheric Delay in meters                *
%2.     ZWD: Zenith Wet Tropospheric Delay  in meters                      *
%3.     ZTD: Zenith Total Tropospheric Delay in meters                     *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%NOTE:                                                                     *     
%    The hydrostatic mapping function coefficients have to be used with the*
%     height dependent Vienna Mapping Function 1  because the              *
%     coefficients refer to zero height.                                   *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% REFERENCE:                                                               *
%           Landskron,D., and Boehm, J.(2018).VMF3/GPT3: re?ned discrete   * 
%           and empirical troposphere mapping functions.J Geod Vol(92),    *
%           pp.349–360

%Original codes by Daniel Landskron (2016)                                 *  
%==========================================================================+
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

%****MODIFIED JULIAN DATE(mjd) & DAY OF YEAR(doy)
[~, mjd ,doy]=utc2JulianDay_DoY(Yr,Mn,Day,H,MIN,SECs);

%CHECK IF SPECIFIED YEAR IS LEAP YEAR
leapYear = ((mod(Yr,4) == 0 & mod(Yr,100) ~= 0) | mod(Yr,400) == 0);

if leapYear == 1 && Mn > 2
    doy = doy + 1;
end

doy = doy + mjd-floor(mjd);   % add decimal places

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
    cosfy = cos(doy/365.25*2*pi);   % coefficient for A1
    coshy = cos(doy/365.25*4*pi);   % coefficient for B1
    sinfy = sin(doy/365.25*2*pi);   % coefficient for A2
    sinhy = sin(doy/365.25*4*pi);   % coefficient for B2
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
Gn_h_grid = grid{11};   % hydrostatic north gradient
Ge_h_grid = grid{12};   % hydrostatic east gradient
Gn_w_grid = grid{13};   % wet north gradient
Ge_w_grid = grid{14};   % wet east gradient

%DETERMINE NUMBER OF STATIONS
nstat = length(lat);%Number of stations

%INITIALIZATION OF NEW VECTORS
p    = zeros([nstat , 1]);
T    = zeros([nstat , 1]);
dT   = zeros([nstat , 1]);
Tm   = zeros([nstat , 1]);
e    = zeros([nstat , 1]);
ah   = zeros([nstat , 1]);
aw   = zeros([nstat , 1]);
la   = zeros([nstat , 1]);
undu = zeros([nstat , 1]);
Gn_h = zeros([nstat , 1]);
Ge_h = zeros([nstat , 1]);
Gn_w = zeros([nstat , 1]);
Ge_w = zeros([nstat , 1]);

%LOOP OVER STATIONs
for k = 1:nstat
    
    %ONLY +VE LONGITUDE IN [degrees]
    if lon(k) < 0
        plon = (lon(k) + 2*pi)*180/pi;
    else
        plon = lon(k)*180/pi;
    end
    % transform to polar distance in degrees
    ppod = (-lat(k) + pi/2)*180/pi; 

    % find the index (line in the grid file) of the nearest point
    ipod = floor((ppod+5)/5); 
    ilon = floor((plon+5)/5);
    
    % normalized (to one) differences, can be positive or negative
    diffpod = (ppod - (ipod*5 - 2.5))/5;
    difflon = (plon - (ilon*5 - 2.5))/5;
    if ipod == 37
        ipod = 36;
    end
    if ilon == 73
		ilon = 1;
    end
    if ilon == 0
        ilon = 72;
    end

    % get the number of the corresponding line
    indx(1) = (ipod - 1)*72 + ilon;
    
    % near the poles: nearest neighbour interpolation, otherwise: bilinear
    bilinear = 0;
    if ppod > 2.5 && ppod < 177.5 
           bilinear = 1;          
    end          
    
    % case of nearest neighbourhood
    if bilinear == 0

        ix = indx(1);
        
        % transforming ellipsoidal height to orthometric height
        undu(k) = u_grid(ix);
        hgt = h_ell(k)-undu(k);
            
        % pressure, temperature at the height of the grid
        T0 = T_grid(ix,1) + T_grid(ix,2)*cosfy + T_grid(ix,3)*sinfy + T_grid(ix,4)*coshy + T_grid(ix,5)*sinhy;
        p0 = p_grid(ix,1) + p_grid(ix,2)*cosfy + p_grid(ix,3)*sinfy + p_grid(ix,4)*coshy + p_grid(ix,5)*sinhy;
         
        % specific humidity
        Q = Q_grid(ix,1) + Q_grid(ix,2)*cosfy + Q_grid(ix,3)*sinfy + Q_grid(ix,4)*coshy + Q_grid(ix,5)*sinhy;
            
        % lapse rate of the temperature
        dT(k) = dT_grid(ix,1) + dT_grid(ix,2)*cosfy + dT_grid(ix,3)*sinfy + dT_grid(ix,4)*coshy + dT_grid(ix,5)*sinhy; 

        % station height - grid height
        redh = hgt - Hs_grid(ix);

        % temperature at station height in Celsius
        T(k) = T0 + dT(k)*redh - 273.15;
        
        % temperature lapse rate in degrees / km
        dT(k) = dT(k)*1000;

        % virtual temperature in Kelvin
        Tv = T0*(1+0.6077*Q);
        
        c = gm*dMtr/(Rg*Tv);
        
        % pressure in hPa
        p(k) = (p0*exp(-c*redh))/100;
            
        % hydrostatic and wet coefficients ah and aw 
        ah(k) = ah_grid(ix,1) + ah_grid(ix,2)*cosfy + ah_grid(ix,3)*sinfy + ah_grid(ix,4)*coshy + ah_grid(ix,5)*sinhy;
        aw(k) = aw_grid(ix,1) + aw_grid(ix,2)*cosfy + aw_grid(ix,3)*sinfy + aw_grid(ix,4)*coshy + aw_grid(ix,5)*sinhy;
		
		% water vapour decrease factor la
        la(k) = la_grid(ix,1) + ...
                la_grid(ix,2)*cosfy + la_grid(ix,3)*sinfy + ...
                la_grid(ix,4)*coshy + la_grid(ix,5)*sinhy; 
		
		% mean temperature Tm
        Tm(k) = Tm_grid(ix,1) + ...
                Tm_grid(ix,2)*cosfy + Tm_grid(ix,3)*sinfy + ...
                Tm_grid(ix,4)*coshy + Tm_grid(ix,5)*sinhy;
            
        % north and east gradients (total, hydrostatic and wet)
        Gn_h(k) = Gn_h_grid(ix,1) + Gn_h_grid(ix,2)*cosfy + Gn_h_grid(ix,3)*sinfy + Gn_h_grid(ix,4)*coshy + Gn_h_grid(ix,5)*sinhy;
        Ge_h(k) = Ge_h_grid(ix,1) + Ge_h_grid(ix,2)*cosfy + Ge_h_grid(ix,3)*sinfy + Ge_h_grid(ix,4)*coshy + Ge_h_grid(ix,5)*sinhy;
        Gn_w(k) = Gn_w_grid(ix,1) + Gn_w_grid(ix,2)*cosfy + Gn_w_grid(ix,3)*sinfy + Gn_w_grid(ix,4)*coshy + Gn_w_grid(ix,5)*sinhy;
        Ge_w(k) = Ge_w_grid(ix,1) + Ge_w_grid(ix,2)*cosfy + Ge_w_grid(ix,3)*sinfy + Ge_w_grid(ix,4)*coshy + Ge_w_grid(ix,5)*sinhy;
		
		% water vapor pressure in hPa
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
        
        % get the number of the line
        indx(2) = (ipod1 - 1)*72 + ilon;  % along same longitude
        indx(3) = (ipod  - 1)*72 + ilon1; % along same polar distance
        indx(4) = (ipod1 - 1)*72 + ilon1; % diagonal
                
        % transforming ellipsoidal height to orthometric height: Hortho = -N + Hell
        undul = u_grid(indx);
        hgt = h_ell(k)-undul;
        
        % pressure, temperature at the height of the grid
        T0 = T_grid(indx,1) + T_grid(indx,2)*cosfy + T_grid(indx,3)*sinfy + T_grid(indx,4)*coshy + T_grid(indx,5)*sinhy;
        p0 = p_grid(indx,1) + p_grid(indx,2)*cosfy + p_grid(indx,3)*sinfy + p_grid(indx,4)*coshy + p_grid(indx,5)*sinhy;
        
        % humidity
        Ql = Q_grid(indx,1) + Q_grid(indx,2)*cosfy + Q_grid(indx,3)*sinfy + Q_grid(indx,4)*coshy + Q_grid(indx,5)*sinhy;
        
        % reduction = stationheight - gridheight
        Hs1 = Hs_grid(indx);
        redh = hgt - Hs1;
        
        % lapse rate of the temperature in degree / m
        dTl = dT_grid(indx,1) + dT_grid(indx,2)*cosfy + dT_grid(indx,3)*sinfy + dT_grid(indx,4)*coshy + dT_grid(indx,5)*sinhy;
        
        % temperature reduction to station height
        Tl = T0 + dTl.*redh - 273.15;
        
        % virtual temperature
        Tv = T0.*(1+0.6077*Ql);
        c = gm*dMtr./(Rg*Tv);
        
        % pressure in hPa
        pl = (p0.*exp(-c.*redh))/100;
        
        % hydrostatic and wet coefficients ah and aw
        ahl = ah_grid(indx,1) + ah_grid(indx,2)*cosfy + ah_grid(indx,3)*sinfy + ah_grid(indx,4)*coshy + ah_grid(indx,5)*sinhy;
        awl = aw_grid(indx,1) + aw_grid(indx,2)*cosfy + aw_grid(indx,3)*sinfy + aw_grid(indx,4)*coshy + aw_grid(indx,5)*sinhy;
        
        % water vapour decrease factor la
        lal = la_grid(indx,1) + la_grid(indx,2)*cosfy + la_grid(indx,3)*sinfy + la_grid(indx,4)*coshy + la_grid(indx,5)*sinhy;
        
        % mean temperature of the water vapor Tm
        Tml = Tm_grid(indx,1) + Tm_grid(indx,2)*cosfy + Tm_grid(indx,3)*sinfy + Tm_grid(indx,4)*coshy + Tm_grid(indx,5)*sinhy;
        
        % north and east gradients (total, hydrostatic and wet)
        Gn_hl = Gn_h_grid(indx,1) + Gn_h_grid(indx,2)*cosfy + Gn_h_grid(indx,3)*sinfy + Gn_h_grid(indx,4)*coshy + Gn_h_grid(indx,5)*sinhy;
        Ge_hl = Ge_h_grid(indx,1) + Ge_h_grid(indx,2)*cosfy + Ge_h_grid(indx,3)*sinfy + Ge_h_grid(indx,4)*coshy + Ge_h_grid(indx,5)*sinhy;
        Gn_wl = Gn_w_grid(indx,1) + Gn_w_grid(indx,2)*cosfy + Gn_w_grid(indx,3)*sinfy + Gn_w_grid(indx,4)*coshy + Gn_w_grid(indx,5)*sinhy;
        Ge_wl = Ge_w_grid(indx,1) + Ge_w_grid(indx,2)*cosfy + Ge_w_grid(indx,3)*sinfy + Ge_w_grid(indx,4)*coshy + Ge_w_grid(indx,5)*sinhy;
        
        % water vapor pressure in hPa
        e0 = Ql.*p0./(0.622+0.378*Ql)/100; % on the grid
        el = e0.*(100.*pl./p0).^(lal+1);  % on the station height - (14) Askne and Nordius, 1987
			
            
        dnpod1 = abs(diffpod); % distance nearer point
        dnpod2 = 1 - dnpod1;   % distance to distant point
        dnlon1 = abs(difflon);
        dnlon2 = 1 - dnlon1;
        
        % pressure
        R1 = dnpod2*pl(1)+dnpod1*pl(2);
        R2 = dnpod2*pl(3)+dnpod1*pl(4);
        p(k) = dnlon2*R1+dnlon1*R2;
            
        % temperature
        R1 = dnpod2*Tl(1)+dnpod1*Tl(2);
        R2 = dnpod2*Tl(3)+dnpod1*Tl(4);
        T(k) = dnlon2*R1+dnlon1*R2;
        
        % temperature in degree per km
        R1 = dnpod2*dTl(1)+dnpod1*dTl(2);
        R2 = dnpod2*dTl(3)+dnpod1*dTl(4);
        dT(k) = (dnlon2*R1+dnlon1*R2)*1000;
            
        % water vapor pressure in hPa
		R1 = dnpod2*el(1)+dnpod1*el(2);
        R2 = dnpod2*el(3)+dnpod1*el(4);
        e(k) = dnlon2*R1+dnlon1*R2;
            
        % ah and aw
        R1 = dnpod2*ahl(1)+dnpod1*ahl(2);
        R2 = dnpod2*ahl(3)+dnpod1*ahl(4);
        ah(k) = dnlon2*R1+dnlon1*R2;
        R1 = dnpod2*awl(1)+dnpod1*awl(2);
        R2 = dnpod2*awl(3)+dnpod1*awl(4);
        aw(k) = dnlon2*R1+dnlon1*R2;
        
        % undulation
        R1 = dnpod2*undul(1)+dnpod1*undul(2);
        R2 = dnpod2*undul(3)+dnpod1*undul(4);
        undu(k) = dnlon2*R1+dnlon1*R2;
		
		% water vapor decrease factor la
        R1 = dnpod2*lal(1)+dnpod1*lal(2);
        R2 = dnpod2*lal(3)+dnpod1*lal(4);
        la(k) = dnlon2*R1+dnlon1*R2;
        
        % gradients
        R1 = dnpod2*Gn_hl(1)+dnpod1*Gn_hl(2);
        R2 = dnpod2*Gn_hl(3)+dnpod1*Gn_hl(4);
        Gn_h(k) = (dnlon2*R1 + dnlon1*R2);
        R1 = dnpod2*Ge_hl(1)+dnpod1*Ge_hl(2);
        R2 = dnpod2*Ge_hl(3)+dnpod1*Ge_hl(4);
        Ge_h(k) = (dnlon2*R1 + dnlon1*R2);
        R1 = dnpod2*Gn_wl(1)+dnpod1*Gn_wl(2);
        R2 = dnpod2*Gn_wl(3)+dnpod1*Gn_wl(4);
        Gn_w(k) = (dnlon2*R1 + dnlon1*R2);
        R1 = dnpod2*Ge_wl(1)+dnpod1*Ge_wl(2);
        R2 = dnpod2*Ge_wl(3)+dnpod1*Ge_wl(4);
        Ge_w(k) = (dnlon2*R1 + dnlon1*R2);
		
		% mean temperature of the water vapor Tm
        R1 = dnpod2*Tml(1)+dnpod1*Tml(2);
        R2 = dnpod2*Tml(3)+dnpod1*Tml(4);
        Tm(k) = dnlon2*R1+dnlon1*R2;
                                        
    end %//if bilinear == 0
    
end%//for k = 1:nstat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF MET PARAMETERS EXTRACTION

%SAVE MAPPING FUNCTION COEFFICIENTS
setappdata(0,'ah_gpt3_5',ah)
setappdata(0,'aw_gpt3_5',aw)
%--------------------------------------------------------------------------
%(2)***COMPUTE ZENITH TROPOSPHERIC DELAYS USING EXTRACTED MET PARAMETERS
%--------------------------------------------------------------------------
%NOTE:
%      DAVIS ET AL DRY MODEL REQUIRES +VE ORTOMETRIC HEIGHT TO PERFORM
%--------------------------------------------------------------------------

%COMPUTE CORRECTION FACTOR FOR LOCAL GRAVITY ACCELERATION 
dgref = 1-0.00266.*cos(2.*lat)-0.00028e-3.*(h_ell-undu);
   gm = 9.784 * dgref;%Acceleration Gravity at the Atmospheric column in m/s^2
   
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

%*****COMPUTE ZENITH WET DELAY(ZWD)    
ZWD = ((1e-6.*(k2p + k3./Tm).*Rd)./(gm.*(la + 1))).*e;%[m]
                
%****FINALLY, COMPUTE ZENITH TOTAL DELAY(ZTD)
ZTD = ZHD + ZWD;%Zenith Total Delay [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GPT3_5 x 5.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%C_4.SUBROUTINE TO COMPUTE MET. PARAMETERS & ZENITH TROPO DELAYs FROM GPT3 1°x1° GRID MODEL
function [ZHD,ZWD,ZTD,T] = GPT3_1_x_1(UTCtime,lat,lon,h_ell,grid,Timevar)
          
%**************************************************************************
%DESCRIPTION:
%            This subroutine determines pressure, temperature, temperature *
%            lapse rate, mean temperature of the water vapor, water vapor  * 
%            pressure, hydrostatic and wet mapping function coefficients   * 
%            ah and aw, water vapour decrease factor, geoid undulation, and* 
%            empirical tropospheric gradients for specific sites near the  *
%            Earth surface. It is based on a  1° x 1° external grid        *
%            file ('gpt3_1.grd') with mean  values as well as sine and     * 
%            cosine amplitudes for the annual and semiannual variation of  *
%            the coefficients.
%USAGE:
%      [ZHD,ZWD,ZTD] = GPT3_1_x_1(UTCtime,lat,lon,hell,grid,Timevar)       *
%INPUTs:                                                                   *
%1.   UTCtime:  UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *                                                                
%2.       lat:  Station geodetic latitude in radians [-pi/2:+pi/2] (n x 1) *
%3.       lon:  Station geodetic longitude in radians [-pi:pi] or [0:2pi]  *
%4.      hell:  ellipsoidal height in [meters] (n x 1)                     *
%5.      grid: grid values in cells extracted from the grid file           *
%6.   Timevar:  case 1: no time variation but static quantities            *
%               case 0: with time variation (annual and semiannual terms)  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% OUTPUTs:
%1.     ZHD: Zenith Hydrostaic Tropospheric Delay in meters                *
%2.     ZWD: Zenith Wet Tropospheric Delay  in meters                      *
%3.    ZTD: Zenith Total Tropospheric Delay in meters                     *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%NOTE:                                                                     *     
%    The hydrostatic mapping function coefficients have to be used with the*
%     height dependent Vienna Mapping Function 1  because the              *
%     coefficients refer to zero height.                                   *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% REFERENCE:                                                               *
%           Landskron,D., and Boehm, J.(2018).VMF3/GPT3: re?ned discrete   * 
%           and empirical troposphere mapping functions.J Geod Vol(92),    *
%           pp.349–360

%Original codes by Daniel Landskron (2016)                                 *  
%==========================================================================+
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

%****COMPUTE MODIFIED JULIAN DATE(mjd) & DAY OF YEAR(doy)
[~, mjd ,doy]=utc2JulianDay_DoY(Yr,Mn,Day,H,MIN,SECs);

%CHECK IF SPECIFIED YEAR IS LEAP YEAR
leapYear = ((mod(Yr,4) == 0 & mod(Yr,100) ~= 0) | mod(Yr,400) == 0);

if leapYear == 1 && Mn > 2
    doy = doy + 1;
end

doy = doy + mjd-floor(mjd);   % add decimal places

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
    cosfy = cos(doy/365.25*2*pi);   % coefficient for A1
    coshy = cos(doy/365.25*4*pi);   % coefficient for B1
    sinfy = sin(doy/365.25*2*pi);   % coefficient for A2
    sinhy = sin(doy/365.25*4*pi);   % coefficient for B2
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
Gn_h_grid = grid{11};   % hydrostatic north gradient
Ge_h_grid = grid{12};   % hydrostatic east gradient
Gn_w_grid = grid{13};   % wet north gradient
Ge_w_grid = grid{14};   % wet east gradient

%DETERMINE NUMBER OF STATIONS
nstat = length(lat);%Number of stations

%INITIALIZATION OF NEW VECTORS
p    = zeros([nstat , 1]);
T    = zeros([nstat , 1]);
dT   = zeros([nstat , 1]);
Tm   = zeros([nstat , 1]);
e    = zeros([nstat , 1]);
ah   = zeros([nstat , 1]);
aw   = zeros([nstat , 1]);
la   = zeros([nstat , 1]);
undu = zeros([nstat , 1]);
Gn_h = zeros([nstat , 1]);
Ge_h = zeros([nstat , 1]);
Gn_w = zeros([nstat , 1]);
Ge_w = zeros([nstat , 1]);

%****LOOP OVER STATIONs
for k = 1:nstat
    
    %Only +VE longitude in [degrees]
    if lon(k) < 0
        plon = (lon(k) + 2*pi)*180/pi;
    else
        plon = lon(k)*180/pi;
    end
    % transform to polar distance in degrees
    ppod = (-lat(k) + pi/2)*180/pi; 

    % find the index (line in the grid file) of the nearest point
    % changed for the 1 degree grid
    ipod = floor(ppod+1); 
    ilon = floor(plon+1);
    
    % normalized (to one) differences, can be positive or negative
    % changed for the 1 degree grid
    diffpod = (ppod - (ipod - 0.5));
    difflon = (plon - (ilon - 0.5));
	% changed for the 1 degree grid
    if ipod == 181
        ipod = 180;
    end
    if ilon == 361
		ilon = 1;
    end
    if ilon == 0
        ilon = 360;
    end

    % get the number of the corresponding line
    % changed for the 1 degree grid
    indx(1) = (ipod - 1)*360 + ilon;
    
    % near the poles: nearest neighbour interpolation, otherwise: bilinear
    % with the 1 degree grid the limits are lower and upper
    bilinear = 0;
    if ppod > 0.5 && ppod < 179.5 
           bilinear = 1;          
    end           
    
    % case of nearest neighbourhood
    if bilinear == 0

        ix = indx(1);
        
        % transforming ellipsoidal height to orthometric height
        undu(k) = u_grid(ix);
        hgt = h_ell(k)-undu(k);
            
        % pressure, temperature at the height of the grid
        T0 = T_grid(ix,1) + T_grid(ix,2)*cosfy + T_grid(ix,3)*sinfy + T_grid(ix,4)*coshy + T_grid(ix,5)*sinhy;
        p0 = p_grid(ix,1) + p_grid(ix,2)*cosfy + p_grid(ix,3)*sinfy + p_grid(ix,4)*coshy + p_grid(ix,5)*sinhy;
         
        % specific humidity
        Q = Q_grid(ix,1) + Q_grid(ix,2)*cosfy + Q_grid(ix,3)*sinfy + Q_grid(ix,4)*coshy + Q_grid(ix,5)*sinhy;
            
        % lapse rate of the temperature
        dT(k) = dT_grid(ix,1) + dT_grid(ix,2)*cosfy + dT_grid(ix,3)*sinfy + dT_grid(ix,4)*coshy + dT_grid(ix,5)*sinhy; 

        % station height - grid height
        redh = hgt - Hs_grid(ix);

        % temperature at station height in Celsius
        T(k) = T0 + dT(k)*redh - 273.15;
        
        % temperature lapse rate in degrees / km
        dT(k) = dT(k)*1000;

        % virtual temperature in Kelvin
        Tv = T0*(1+0.6077*Q);
        
        c = gm*dMtr/(Rg*Tv);
        
        % pressure in hPa
        p(k) = (p0*exp(-c*redh))/100;
            
        % hydrostatic and wet coefficients ah and aw 
        ah(k) = ah_grid(ix,1) + ah_grid(ix,2)*cosfy + ah_grid(ix,3)*sinfy + ah_grid(ix,4)*coshy + ah_grid(ix,5)*sinhy;
        aw(k) = aw_grid(ix,1) + aw_grid(ix,2)*cosfy + aw_grid(ix,3)*sinfy + aw_grid(ix,4)*coshy + aw_grid(ix,5)*sinhy;
		
		% water vapour decrease factor la
        la(k) = la_grid(ix,1) + ...
                la_grid(ix,2)*cosfy + la_grid(ix,3)*sinfy + ...
                la_grid(ix,4)*coshy + la_grid(ix,5)*sinhy; 
		
		% mean temperature Tm
        Tm(k) = Tm_grid(ix,1) + ...
                Tm_grid(ix,2)*cosfy + Tm_grid(ix,3)*sinfy + ...
                Tm_grid(ix,4)*coshy + Tm_grid(ix,5)*sinhy;
            
        % north and east gradients (total, hydrostatic and wet)
        Gn_h(k) = Gn_h_grid(ix,1) + Gn_h_grid(ix,2)*cosfy + Gn_h_grid(ix,3)*sinfy + Gn_h_grid(ix,4)*coshy + Gn_h_grid(ix,5)*sinhy;
        Ge_h(k) = Ge_h_grid(ix,1) + Ge_h_grid(ix,2)*cosfy + Ge_h_grid(ix,3)*sinfy + Ge_h_grid(ix,4)*coshy + Ge_h_grid(ix,5)*sinhy;
        Gn_w(k) = Gn_w_grid(ix,1) + Gn_w_grid(ix,2)*cosfy + Gn_w_grid(ix,3)*sinfy + Gn_w_grid(ix,4)*coshy + Gn_w_grid(ix,5)*sinhy;
        Ge_w(k) = Ge_w_grid(ix,1) + Ge_w_grid(ix,2)*cosfy + Ge_w_grid(ix,3)*sinfy + Ge_w_grid(ix,4)*coshy + Ge_w_grid(ix,5)*sinhy;
		
		% water vapor pressure in hPa
		e0 = Q*p0/(0.622+0.378*Q)/100; % on the grid
		e(k) = e0*(100*p(k)/p0)^(la(k)+1);   % on the station height - (14) Askne and Nordius, 1987
		
     else % bilinear interpolation
        
        ipod1 = ipod + sign(diffpod);
        ilon1 = ilon + sign(difflon);
		% changed for the 1 degree grid
        if ilon1 == 361
            ilon1 = 1;
        end
        if ilon1 == 0
            ilon1 = 360;
        end
        
        % get the number of the line
		% changed for the 1 degree grid
        indx(2) = (ipod1 - 1)*360 + ilon;  % along same longitude
        indx(3) = (ipod  - 1)*360 + ilon1; % along same polar distance
        indx(4) = (ipod1 - 1)*360 + ilon1; % diagonal
                
        % transforming ellipsoidal height to orthometric height: Hortho = -N + Hell
        undul = u_grid(indx);
        hgt = h_ell(k)-undul;
        
        % pressure, temperature at the height of the grid
        T0 = T_grid(indx,1) + T_grid(indx,2)*cosfy + T_grid(indx,3)*sinfy + T_grid(indx,4)*coshy + T_grid(indx,5)*sinhy;
        p0 = p_grid(indx,1) + p_grid(indx,2)*cosfy + p_grid(indx,3)*sinfy + p_grid(indx,4)*coshy + p_grid(indx,5)*sinhy;
        
        % humidity
        Ql = Q_grid(indx,1) + Q_grid(indx,2)*cosfy + Q_grid(indx,3)*sinfy + Q_grid(indx,4)*coshy + Q_grid(indx,5)*sinhy;
        
        % reduction = stationheight - gridheight
        Hs1 = Hs_grid(indx);
        redh = hgt - Hs1;
        
        % lapse rate of the temperature in degree / m
        dTl = dT_grid(indx,1) + dT_grid(indx,2)*cosfy + dT_grid(indx,3)*sinfy + dT_grid(indx,4)*coshy + dT_grid(indx,5)*sinhy;
        
        % temperature reduction to station height
        Tl = T0 + dTl.*redh - 273.15;
        
        % virtual temperature
        Tv = T0.*(1+0.6077*Ql);
        c = gm*dMtr./(Rg*Tv);
        
        % pressure in hPa
        pl = (p0.*exp(-c.*redh))/100;
        
        % hydrostatic and wet coefficients ah and aw
        ahl = ah_grid(indx,1) + ah_grid(indx,2)*cosfy + ah_grid(indx,3)*sinfy + ah_grid(indx,4)*coshy + ah_grid(indx,5)*sinhy;
        awl = aw_grid(indx,1) + aw_grid(indx,2)*cosfy + aw_grid(indx,3)*sinfy + aw_grid(indx,4)*coshy + aw_grid(indx,5)*sinhy;
        
        % water vapour decrease factor la
        lal = la_grid(indx,1) + la_grid(indx,2)*cosfy + la_grid(indx,3)*sinfy + la_grid(indx,4)*coshy + la_grid(indx,5)*sinhy;
        
        % mean temperature of the water vapor Tm
        Tml = Tm_grid(indx,1) + Tm_grid(indx,2)*cosfy + Tm_grid(indx,3)*sinfy + Tm_grid(indx,4)*coshy + Tm_grid(indx,5)*sinhy;
        
        % north and east gradients (total, hydrostatic and wet)
        Gn_hl = Gn_h_grid(indx,1) + Gn_h_grid(indx,2)*cosfy + Gn_h_grid(indx,3)*sinfy + Gn_h_grid(indx,4)*coshy + Gn_h_grid(indx,5)*sinhy;
        Ge_hl = Ge_h_grid(indx,1) + Ge_h_grid(indx,2)*cosfy + Ge_h_grid(indx,3)*sinfy + Ge_h_grid(indx,4)*coshy + Ge_h_grid(indx,5)*sinhy;
        Gn_wl = Gn_w_grid(indx,1) + Gn_w_grid(indx,2)*cosfy + Gn_w_grid(indx,3)*sinfy + Gn_w_grid(indx,4)*coshy + Gn_w_grid(indx,5)*sinhy;
        Ge_wl = Ge_w_grid(indx,1) + Ge_w_grid(indx,2)*cosfy + Ge_w_grid(indx,3)*sinfy + Ge_w_grid(indx,4)*coshy + Ge_w_grid(indx,5)*sinhy;
        
        % water vapor pressure in hPa
        e0 = Ql.*p0./(0.622+0.378*Ql)/100; % on the grid
        el = e0.*(100.*pl./p0).^(lal+1);  % on the station height - (14) Askne and Nordius, 1987
			
            
        dnpod1 = abs(diffpod); % distance nearer point
        dnpod2 = 1 - dnpod1;   % distance to distant point
        dnlon1 = abs(difflon);
        dnlon2 = 1 - dnlon1;
        
        % pressure
        R1 = dnpod2*pl(1)+dnpod1*pl(2);
        R2 = dnpod2*pl(3)+dnpod1*pl(4);
        p(k) = dnlon2*R1+dnlon1*R2;
            
        % temperature
        R1 = dnpod2*Tl(1)+dnpod1*Tl(2);
        R2 = dnpod2*Tl(3)+dnpod1*Tl(4);
        T(k) = dnlon2*R1+dnlon1*R2;
        
        % temperature in degree per km
        R1 = dnpod2*dTl(1)+dnpod1*dTl(2);
        R2 = dnpod2*dTl(3)+dnpod1*dTl(4);
        dT(k) = (dnlon2*R1+dnlon1*R2)*1000;
            
        % water vapor pressure in hPa
		R1 = dnpod2*el(1)+dnpod1*el(2);
        R2 = dnpod2*el(3)+dnpod1*el(4);
        e(k) = dnlon2*R1+dnlon1*R2;
            
        % ah and aw
        R1 = dnpod2*ahl(1)+dnpod1*ahl(2);
        R2 = dnpod2*ahl(3)+dnpod1*ahl(4);
        ah(k) = dnlon2*R1+dnlon1*R2;
        R1 = dnpod2*awl(1)+dnpod1*awl(2);
        R2 = dnpod2*awl(3)+dnpod1*awl(4);
        aw(k) = dnlon2*R1+dnlon1*R2;
        
        % undulation
        R1 = dnpod2*undul(1)+dnpod1*undul(2);
        R2 = dnpod2*undul(3)+dnpod1*undul(4);
        undu(k) = dnlon2*R1+dnlon1*R2;
		
		% water vapor decrease factor la
        R1 = dnpod2*lal(1)+dnpod1*lal(2);
        R2 = dnpod2*lal(3)+dnpod1*lal(4);
        la(k) = dnlon2*R1+dnlon1*R2;
        
        % gradients
        R1 = dnpod2*Gn_hl(1)+dnpod1*Gn_hl(2);
        R2 = dnpod2*Gn_hl(3)+dnpod1*Gn_hl(4);
        Gn_h(k) = (dnlon2*R1 + dnlon1*R2);
        R1 = dnpod2*Ge_hl(1)+dnpod1*Ge_hl(2);
        R2 = dnpod2*Ge_hl(3)+dnpod1*Ge_hl(4);
        Ge_h(k) = (dnlon2*R1 + dnlon1*R2);
        R1 = dnpod2*Gn_wl(1)+dnpod1*Gn_wl(2);
        R2 = dnpod2*Gn_wl(3)+dnpod1*Gn_wl(4);
        Gn_w(k) = (dnlon2*R1 + dnlon1*R2);
        R1 = dnpod2*Ge_wl(1)+dnpod1*Ge_wl(2);
        R2 = dnpod2*Ge_wl(3)+dnpod1*Ge_wl(4);
        Ge_w(k) = (dnlon2*R1 + dnlon1*R2);
		
		% mean temperature of the water vapor Tm
        R1 = dnpod2*Tml(1)+dnpod1*Tml(2);
        R2 = dnpod2*Tml(3)+dnpod1*Tml(4);
        Tm(k) = dnlon2*R1+dnlon1*R2;
                    
    end %//if bilinear == 0
    
end%//for k = 1:nstat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF MET PARAMETERS EXTRACTION

%SAVE MAPPING FUNCTION COEFFICIENTS
setappdata(0,'ah_gpt3_1',ah)
setappdata(0,'aw_gpt3_1',aw)

%--------------------------------------------------------------------------
%(2)***COMPUTE ZENITH TROPOSPHERIC DELAYS USING EXTRACTED MET PARAMETERS
%--------------------------------------------------------------------------
%NOTE:
%      DAVIS ET AL DRY MODEL REQUIRES +VE ORTOMETRIC HEIGHT TO PERFORM
%--------------------------------------------------------------------------

%COMPUTE CORRECTION FACTOR FOR LOCAL GRAVITY ACCELERATION 
dgref = 1-0.00266.*cos(2.*lat)-0.00028e-3.*(h_ell-undu);
   gm = 9.784 * dgref;%Acceleration Gravity at the Atmospheric column in m/s^2
   
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

%*****COMPUTE ZENITH WET DELAY(ZWD)    
ZWD = ((1e-6.*(k2p + k3./Tm).*Rd)./(gm.*(la + 1))).*e;%[m]
                
%****FINALLY, COMPUTE ZENITH TOTAL DELAY(ZTD)
ZTD = ZHD + ZWD;%Zenith Total Delay [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GPT3_1_x_1.m  %%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%C_5.SUBROUTINE TO COMPUTE MET. PARAMETERS & ZENITH TROPO DELAYs FROM GPT2 5°x5° GRID MODEL
function [ZHD,ZWD,ZTD,T] = GPT2_5_x_5(UTCtime,lat,lon,hell,grid,Timevar)
        
%**************************************************************************
%DESCRIPTION:                                                              *
%            This subroutine determines pressure, temperature, temperature *
%            lapse rate,water vapor pressure, hydrostatic and wet mapping  * 
%            function coefficients ah and aw, water vapour decrease        * 
%            factor and geoid undulation for specific sites near the Earth * 
%            surface. It is based on a 5° x 5° external grid file          *
%            ('gpt2_5.grd') with mean  values as well as sine and cosine   *
%             amplitudes for the annual  and semiannual variation of the   * 
%             coefficients.                                                *
%USAGE:                                                                    *
%      [ZHD,ZWD,ZTD] = GPT2_5_x_5(UTCtime,lat,lon,hell,grid,Timevar)       *    
%INPUTs:                                                                   *
%1.   UTCtime:  UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *                                                                
%2.       lat:  Station geodetic latitude in radians [-pi/2:+pi/2] (n x 1) *
%3.       lon:  Station geodetic longitude in radians [-pi:pi] or [0:2pi]  *
%4.      hell:  ellipsoidal height in [meters] (n x 1)                     *
%5.      grid:  grid values in cells extracted from the grid file          *
%6    Timevar:  case 1: no time variation but static quantities            *
%               case 0: with time variation (annual and semiannual terms)  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% OUTPUTs:                                                                 *
%1.     ZHD: Zenith Hydrostaic Tropospheric Delay in meters                *
%2.     ZWD: Zenith Wet Tropospheric Delay  in meters                      *
%3.     ZTD: Zenith Total Tropospheric Delay in meters                     *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%NOTE:                                                                     *     
%    The hydrostatic mapping function coefficients have to be used with the*
%     height dependent Vienna Mapping Function 1  because the              *
%     coefficients refer to zero height.                                   *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%%REFERENCE:                                                               *
%          Lagler, K., Schindelegger, M., Böhm, J., Krásná, H., Nilsson, T.*  
%         (2013),GPT2: Empirical slant delay model for radio space geodetic*  
%          techniques,Geophys. Res. Lett., Vol. 40, pp. 1069–1073, DOI:    * 
%          10.1002/grl.50288.                                              *
%Original codes by Böhm et al 2014                                         *
%Modified by: Samuel Osah,Msc Geomatic Engineering ,2016                   *    
%             Email: osahsamuel@yahoo.ca                                   *
%             Tel:+233 (0)246137410/+233 (0)509438484                      *     
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

%DETERMINE NUMBER OF STATIONS
nstat = length(lat);%Number of stations

%INITIALIZE NEW VECTORs
p =  zeros([nstat, 1]);
T =  zeros([nstat, 1]);
dT = zeros([nstat, 1]);
e =  zeros([nstat, 1]);
ah = zeros([nstat, 1]);
aw = zeros([nstat, 1]);
undu = zeros([nstat, 1]);

%LOOP OVER STATIONs
for k = 1:nstat
    
    %Only +VE longitude in [degrees]    
    if lon(k) < 0
        plon = (lon(k) + 2*pi)*180/pi;
    else
        plon = lon(k)*180/pi;
    end
    %Transform to polar distance in degrees
    ppod = (-lat(k) + pi/2)*180/pi; 

    % find the index (line in the grid file) of the nearest point
    ipod = floor((ppod+5)/5); 
    ilon = floor((plon+5)/5);
    
    %Normalized (to one) differences, can be positive or negative
    diffpod = (ppod - (ipod*5 - 2.5))/5;
    difflon = (plon - (ilon*5 - 2.5))/5;
    %added by HCY
    if ipod == 37
        ipod = 36;
    end

    %Get the number of the corresponding line
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

        %Temperature at station height in [Celsius]
        T(k) = T0 + dT(k)*redh - 273.15;
        
        %Temperature lapse rate in [degrees / km]
        dT(k) = dT(k)*1000;

        %Virtual temperature in Kelvin
        Tv = T0*(1+0.6077*Q);        
        c = gm*dMtr/(Rg*Tv);
        
        %Pressure in [hPa]
        p(k) = (p0*exp(-c*redh))/100;
        
        %Water vapour pressure in [hPa]
        e(k) = (Q*p(k))/(0.622+0.378*Q);
            
        %Hydrostatic coefficient ah 
        ah(k) = ah_grid(ix,1) + ...
                ah_grid(ix,2)*cosfy + ah_grid(ix,3)*sinfy+ ...
                ah_grid(ix,4)*coshy + ah_grid(ix,5)*sinhy;
            
        %Wet coefficient aw
        aw(k) = aw_grid(ix,1) + ...
                aw_grid(ix,2)*cosfy + aw_grid(ix,3)*sinfy + ...
                aw_grid(ix,4)*coshy + aw_grid(ix,5)*sinhy;           
                    
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
            %Hortho = -N + Hell
            undul(l) = u_grid(indx(l)); %#ok<*AGROW>
            hgt = hell(k)-undul(l);
        
            %Pressure, temperature at the heigtht of the grid
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
        
        %Temperature in [degree per km]
        R1 = dnpod2*dTl(1)+dnpod1*dTl(2);
        R2 = dnpod2*dTl(3)+dnpod1*dTl(4);
        dT(k) = (dnlon2*R1+dnlon1*R2)*1000;
            
        %Humidity
        R1 = dnpod2*Ql(1)+dnpod1*Ql(2);
        R2 = dnpod2*Ql(3)+dnpod1*Ql(4);
        Q = dnlon2*R1+dnlon1*R2;
        e(k) = (Q*p(k))/(0.622+0.378*Q);
            
        %Hydrostatic
        R1 = dnpod2*ahl(1)+dnpod1*ahl(2);
        R2 = dnpod2*ahl(3)+dnpod1*ahl(4);
        ah(k) = dnlon2*R1+dnlon1*R2;
           
        %Wet
        R1 = dnpod2*awl(1)+dnpod1*awl(2);
        R2 = dnpod2*awl(3)+dnpod1*awl(4);
        aw(k) = dnlon2*R1+dnlon1*R2;
        
        %Undulation (N)
        R1 = dnpod2*undul(1)+dnpod1*undul(2);
        R2 = dnpod2*undul(3)+dnpod1*undul(4);
        undu(k) = dnlon2*R1+dnlon1*R2;
                    
    end %//if bilinear == 0
    
end%//for k = 1:nstat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF MET PARAMETERS EXTRACTION

%SAVE MAPPING FUNCTION COEFFICIENTS
setappdata(0,'ah_gpt2_5',ah)
setappdata(0,'aw_gpt2_5',aw)

%--------------------------------------------------------------------------
%(2)***COMPUTE ZENITH TROPOSPHERIC DELAYS USING EXTRACTED MET PARAMETERS
%--------------------------------------------------------------------------

%**************COMPUTE ZENITH HYDROSTATIC DELAY(ZHD)
%              =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%--------------------------------------------------------------------------
%NOTE:
%     USE OF SAASTAMOINEN MODEL REQUIRES +VE ORTOMETRIC HEIGHT TO PERFORM
%--------------------------------------------------------------------------
%COMPUTE CORRECTION FACTOR FOR LOCAL GRAVITY ACCELERATION 
dgref = 1-0.00266.*cos(2.*lat)-0.00028e-3.*(hell-undu);
 
%CONVERT TEMPERATURE(T) IN Celcius(°C) to Kelvin(K) 
 Temp = T + 273.15; %Temperature in kelvin(K)
 
%USE OF SAASTAMOINEN DRY MODEL   
ZHD =(0.002277 .* p)./dgref; %Zenith Hydrostatic or dry delay [m]
            
%**************COMPUTE ZENITH WET DELAY(ZHD)
%              =-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%*****COMPUTE ZENITH WET DELAY(ZWD) USING OF SAASTAMOINEN MODEL 
ZWD = (0.002277.*((1255 ./ Temp) + 0.05).* e)./dgref;% Zenith wet delay[m]
                
%****FINALLY, COMPUTE ZENITH TOTAL DELAY(ZTD)
ZTD = ZHD + ZWD;%Zenith Total Delay [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GPT2_5_x_5.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%D.SUBROUTINE TO COMPUTE MET. PARAMETERs & ZENITH TROPO DELAYs-UNB3m,EGNOS
function [ZHD,ZWD,ZTD,T]=ZTD_UNB3m_EGNOS_MOPS(UTCtime,lat,hgt,Drymodel,Wetmodel)

%*********INTERPOLATE METEOROLOGICAL PARAMETERs & COMPUTE ZHD,ZWD,ZTD
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
      
     %*******INITIALIZE OUTPUTs
     ZHD = zeros(nrow_time,1);
     [ZHD,ZWD,ZTD,T] = deal(ZHD);%create copy of ZHD in ZWD,ZTD,T
   
     for i=1:nrow_time
      
         if any([strncmpi(Drymodel,Wetmodel,5),strfind(Drymodel,Wetmodel)])
             
            %****Compute ZHD,ZWD ZTD
            [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(i),hgt(i),Drymodel); 
            
         elseif all([any([strncmpi(Drymodel,'UNB3m',5),strncmpi(Drymodel,'EGNOS',5),strncmpi(Drymodel,'MOPS',4)]),any([~strncmpi(Wetmodel,'UNB3m',5),~strncmpi(Wetmodel,'EGNOS',5),~strncmpi(Wetmodel,'MOPS',4)])])
                
                %****Compute ZHD,ZWD ZTD
                [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(i),hgt(i),Drymodel); 
            
         elseif all([any([~strncmpi(Drymodel,'UNB3m',5),~strncmpi(Drymodel,'EGNOS',5),~strncmpi(Drymodel,'MOPS',4)]),any([strncmpi(Wetmodel,'UNB3m',5),strncmpi(Wetmodel,'EGNOS',5),strncmpi(Wetmodel,'MOPS',4)])])
                
                %****Compute ZHD,ZWD ZTD
                [ZHD(i,1),ZWD(i,1),ZTD(i,1),T(i,1)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(i),hgt(i),Wetmodel); 
                             
         else    
              if all([strncmpi(Drymodel,'UNB3m',5),strncmpi(Wetmodel,'EGNOS',5)])      
                
                 %****Compute ZHD USING UNB3m MODEL 
                [ZHD(i,1),~,~,T(i,1)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(i),hgt(i),Drymodel); 
                
                %****Compute ZWD USING EGNOS MODEL
                [~,ZWD(i,1),~,T(i,1)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(i),hgt(i),Wetmodel);
                
                %*****COMPUTE ZTD
                ZTD(i,1) = ZHD(i,1) + ZWD(i,1);
                
              elseif all([strncmpi(Drymodel,'UNB3m',5),strncmpi(Wetmodel,'MOPS',4)])      
                
                     %****Compute ZHD USING UNB3m MODEL
                     [ZHD(i,1),~,~,T(i,1)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(i),hgt(i),Drymodel); 
                
                     %****Compute ZWD USING MOPS MODEL
                     [~,ZWD(i,1),~,T(i,1)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(i),hgt(i),Wetmodel);
                
                     %****COMPUTE ZTD
                     ZTD(i,1) = ZHD(i,1) + ZWD(i,1);
                     
                     
              elseif all([strncmpi(Drymodel,'EGNOS',5),strncmpi(Wetmodel,'UNB3m',5)])      
                
                     %****Compute ZHD USING EGNOS MODEL
                     [ZHD(i,1),~,~,T(i,1)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(i),hgt(i),Drymodel); 
                
                    %****Compute ZWD USING UNB3m MODEL
                    [~,ZWD(i,1),~,T(i,1)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(i),hgt(i),Wetmodel);
                
                    %*****COMPUTE ZTD
                    ZTD(i,1) = ZHD(i,1) + ZWD(i,1);
                
              elseif all([strncmpi(Drymodel,'MOPS',4),strncmpi(Wetmodel,'UNB3m',5)])      
                
                     %****Compute ZHD USING MOPS MODEL
                     [ZHD(i,1),~,~,T(i,1)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(i),hgt(i),Drymodel); 
                
                     %****Compute ZWD USING UNB3m MODEL
                     [~,ZWD(i,1),~,T(i,1)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(i),hgt(i),Wetmodel);
                
                     %****COMPUTE ZTD
                     ZTD(i,1) = ZHD(i,1) + ZWD(i,1);
                     
              end %//if all([strncmpi(Drymodel,'UNB3m',5),strncmpi(Wetmodel,'EGNOS',5)]) 
                                   
         end %//if any([strncmpi(Drymodel,Wetmodel),strfind(Drymodel,Wetmodel)])
         
     end %//for i=1:nrow_time
     
  else
      %*****INITIALIZE OUTPUTs 
      ZHD = zeros(nrow_pos,nrow_time);
      [ZWD,ZTD,T] = deal(ZHD);%create copy of ZHD in ZWD,ZTD,T
      
      for i=1:nrow_time
          
          for j=1:nrow_pos
            
              if any([strncmpi(Drymodel,Wetmodel,5),strfind(Drymodel,Wetmodel)])
             
                 %****Compute ZHD,ZWD ZTD
                 [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Drymodel); 
            
              elseif all([any([strncmpi(Drymodel,'UNB3m',5),strncmpi(Drymodel,'EGNOS',5),strncmpi(Drymodel,'MOPS',4)]),any([~strncmpi(Wetmodel,'UNB3m',5),~strncmpi(Wetmodel,'EGNOS',5),~strncmpi(Wetmodel,'MOPS',4)])])
                
                     %****Compute ZHD,ZWD ZTD
                    [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Drymodel); 
            
         
              elseif all([any([~strncmpi(Drymodel,'UNB3m',5),~strncmpi(Drymodel,'EGNOS',5),~strncmpi(Drymodel,'MOPS',4)]),any([strncmpi(Wetmodel,'UNB3m',5),strncmpi(Wetmodel,'EGNOS',5),strncmpi(Wetmodel,'MOPS',4)])])
                
                     %****Compute ZHD,ZWD ZTD
                     [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Wetmodel); 
                             
              else     
                   if all([strncmpi(Drymodel,'UNB3m',5),strncmpi(Wetmodel,'EGNOS',5)])      
                
                      %****Compute ZHD USING UNB3m MODEL 
                      [ZHD(j,i),~,~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Drymodel); 
                
                      %****Compute ZWD USING EGNOS MODEL
                      [~,ZWD(j,i),~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Wetmodel);
                
                      %*****COMPUTE ZTD
                      ZTD(j,i) = ZHD(j,i) + ZWD(j,i);
                
                   elseif all([strncmpi(Drymodel,'UNB3m',5),strncmpi(Wetmodel,'MOPS',4)])      
                
                          %****Compute ZHD USING UNB3m MODEL
                          [ZHD(j,i),~,~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Drymodel); 
                
                          %****Compute ZWD USING MOPS MODEL
                          [~,ZWD(j,i),~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Wetmodel);
                
                          %****COMPUTE ZTD
                          ZTD(j,i) = ZHD(j,i) + ZWD(j,i);
                     
                   elseif all([strncmpi(Drymodel,'EGNOS',5),strncmpi(Wetmodel,'UNB3m',5)])      
                
                          %****Compute ZHD USING EGNOS MODEL
                          [ZHD(j,i),~,~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Drymodel); 
                
                          %****Compute ZWD USING UNB3m MODEL
                          [~,ZWD(j,i),~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Wetmodel);
                
                          %*****COMPUTE ZTD
                          ZTD(j,i) = ZHD(j,i) + ZWD(j,i);
                
                   elseif all([strncmpi(Drymodel,'MOPS',4),strncmpi(Wetmodel,'UNB3m',5)])      
                
                          %****Compute ZHD USING MOPS MODEL
                          [ZHD(j,i),~,~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Drymodel); 
                
                          %****Compute ZWD USING UNB3m MODEL
                          [~,ZWD(j,i),~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Wetmodel);
                
                          %****COMPUTE ZTD
                          ZTD(j,i) = ZHD(i,1) + ZWD(j,i);
                     
                   end  %//if all([strncmpi(Drymodel,'UNB3m',5),strncmpi(Wetmodel,'EGNOS',5)]) 
                                   
              end  %//if any([strncmpi(Drymodel,Wetmodel),strfind(Drymodel,Wetmodel)])
              
          end %//for j=1:nrow_pos 
          
      end %//for i=1:nrow_time
      
  end %// if isequal(Identical,1)
        
else
    %*****INITIALIZE OUTPUTs 
    ZHD = zeros(nrow_pos,nrow_time);
   [ZWD,ZTD,T] = deal(ZHD);%create copy of ZHD in ZWD,ZTD,T
    
    for i=1:nrow_time
        
        for j=1:nrow_pos
            
            if any([strncmpi(Drymodel,Wetmodel,5),strfind(Drymodel,Wetmodel)])
             
                 %****Compute ZHD,ZWD ZTD
                 [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Drymodel); 
            
              elseif all([any([strncmpi(Drymodel,'UNB3m',5),strncmpi(Drymodel,'EGNOS',5),strncmpi(Drymodel,'MOPS',4)]),any([~strncmpi(Wetmodel,'UNB3m',5),~strncmpi(Wetmodel,'EGNOS',5),~strncmpi(Wetmodel,'MOPS',4)])])
                
                     %****Compute ZHD,ZWD ZTD
                    [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Drymodel); 
            
         
              elseif all([any([~strncmpi(Drymodel,'UNB3m',5),~strncmpi(Drymodel,'EGNOS',5),~strncmpi(Drymodel,'MOPS',4)]),any([strncmpi(Wetmodel,'UNB3m',5),strncmpi(Wetmodel,'EGNOS',5),strncmpi(Wetmodel,'MOPS',4)])])
                
                     %****Compute ZHD,ZWD ZTD
                     [ZHD(j,i),ZWD(j,i),ZTD(j,i),T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Wetmodel); 
                             
              else     
                   if all([strncmpi(Drymodel,'UNB3m',5),strncmpi(Wetmodel,'EGNOS',5)])      
                
                      %****Compute ZHD USING UNB3m MODEL 
                      [ZHD(j,i),~,~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Drymodel); 
                
                      %****Compute ZWD USING EGNOS MODEL
                      [~,ZWD(j,i),~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Wetmodel);
                
                      %*****COMPUTE ZTD
                      ZTD(j,i) = ZHD(j,i) + ZWD(j,i);
                
                   elseif all([strncmpi(Drymodel,'UNB3m',5),strncmpi(Wetmodel,'MOPS',4)])      
                
                          %****Compute ZHD USING UNB3m MODEL
                          [ZHD(j,i),~,~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Drymodel); 
                
                          %****Compute ZWD USING MOPS MODEL
                          [~,ZWD(j,i),~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Wetmodel);
                
                          %****COMPUTE ZTD
                          ZTD(j,i) = ZHD(j,i) + ZWD(j,i);
                     
                   elseif all([strncmpi(Drymodel,'EGNOS',5),strncmpi(Wetmodel,'UNB3m',5)])      
                
                          %****Compute ZHD USING EGNOS MODEL
                          [ZHD(j,i),~,~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Drymodel); 
                
                          %****Compute ZWD USING UNB3m MODEL
                          [~,ZWD(j,i),~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Wetmodel);
                
                          %*****COMPUTE ZTD
                          ZTD(j,i) = ZHD(j,i) + ZWD(j,i);
                
                   elseif all([strncmpi(Drymodel,'MOPS',4),strncmpi(Wetmodel,'UNB3m',5)])      
                
                          %****Compute ZHD USING MOPS MODEL
                          [ZHD(j,i),~,~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Drymodel); 
                
                          %****Compute ZWD USING UNB3m MODEL
                          [~,ZWD(j,i),~,T(j,i)]=getZTD_UNB3m_EGNOS_MOPS(UTCtime(i,:),lat(j),hgt(j),Wetmodel);
                
                          %****COMPUTE ZTD
                          ZTD(j,i) = ZHD(i,1) + ZWD(j,i);
                     
                   end  %//if all([strncmpi(Drymodel,'UNB3m',5),strncmpi(Wetmodel,'EGNOS',5)]) 
                                   
            end   %//if any([strncmpi(Drymodel,Wetmodel),strfind(Drymodel,Wetmodel)])
              
        end  %//for j=1:nrow_pos 
          
    end  %//for i=1:nrow_time
      
end %//if isequal(nrow_time,nrow_pos)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZTD_UNB3m_EGNOS_MOPS.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%D_1.SUBROUTINE TO COMPUTE MET. PARAMETERs & ZENITH TROPO DELAYs-UNB3m,EGNOS
function [ZHD,ZWD,ZTD,T]=getZTD_UNB3m_EGNOS_MOPS(UTCtime,lat,hgt,model)
%**************************************************************************
%DESCRIPTION:                                                              * 
%           "getZTD_UNB3m_EGNOS_MOPS"computes surface Meteorological parameters*
%           temperature,pressure,relative humidity or water vapor pressure,*
%           mean temperature of water vapor and water vapor lapse rates    *
%           (beta & lambda) for a given latitude, height and day of year.  * 
%           The Estimated Meteorological parameters are then use to compute*
%           the Zenith Hydrostatic, Wet and Total Tropospheric Delays      *
%           (ZHD,ZWD,ZTD).                                                 *
%USAGE:                                                                    *
%      [ZHD,ZWD,ZTD]=getMETparaZTD_table(UTCtime,lat,hgt,model)            *
%                                                                          *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%****INPUTs:                                                               *
%1.       UTCtime:.........UTC time in [Year,Month,Day,Hour,Minute,Seconds]*                         
%2.           lat:.........Station geodetic latitude in [degrees]          *
%3.           hgt:.........Station Orthometric height in [meters]          *
%4.         model:.........Model type in string eg:'UNB3m','EGNOS','MOPS'  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%***OUTPUTs:                                                               *                 * 
%1.       ZHD : Zenith Hydrostaic Tropospheric Delay in meters             *
%2.       ZWD : Zenith Wet Tropospheric Delay  in meters                   *
%3.       ZTD : Zenith Total Tropospheric Delay in meters                  *
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
%**************************************************************************
%GET UTC TIME COMPONENTs
%1.******UTCdate
  Yr = UTCtime(:,1);%get Year
  Mn = UTCtime(:,2);%get Month
 Day = UTCtime(:,3);%get Day
%2.UTCtime
   H = UTCtime(:,4);%get Hour
 MIN = UTCtime(:,5);%get Minute
SECs = UTCtime(:,6);%get Seconds

%*******************Initialize UNB3m LOOK-UP TABLE
%(1)AVG:AVERAGE METEOROLOGICAL PARAMETERS
AVG1=[ 15.0  1013.25  299.65  75.00  6.30e-3  2.77
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

AMP1=[ 15.0   0.00   0.00   0.00  0.00e-3  0.00
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

%*******************Initialize EGNOS & MOPS LOOK-UP TABLE
%(1)AVG:AVERAGE METEOROLOGICAL PARAMETERS
AVG2=[ 15.0  1013.25  299.65  26.31  6.30e-3  2.77
       30.0  1017.25  294.15  21.79  6.05e-3  3.15
       45.0  1015.75  283.15  11.66  5.58e-3  2.57
       60.0  1011.75  272.15  06.78  5.39e-3  1.81
       75.0  1013.00  263.65  04.11  4.53e-3  1.55];
      %lat     Po      To      E0     beta0  lambda0
%%WHERE:
%       E0: average water vapor pressure in [mbar]

AMP2=[ 15.0   0.00   0.00   0.00  0.00e-3  0.00
       30.0  -3.75   7.00   8.85  0.25e-3  0.33
       45.0  -2.25  11.00   7.24  0.32e-3  0.46
       60.0  -1.75  15.00   5.36  0.81e-3  0.74
       75.0  -0.50  14.50   3.39  0.62e-3  0.30];
     % lat     dP    dT     dE0    dbeta  dlambda
%WHERE:
%     dE0 : Seasonal variation in water vapor pressure in [mbar]
    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%*********INTERPOLATE METEOROLOGICAL PARAMETERs
if strncmpi(model,'UNB3m',5)
   %1.Compute Average & seasonal Values
   [intAVG,intAMP] =interpolate_table(lat,AVG1,AMP1);

elseif any([strncmpi(model,'EGNOS',5),strncmpi(model,'MOPS',4)])
       %1.Compute Average & seasonal Values
       [intAVG,intAMP] =interpolate_table(lat,AVG2,AMP2);       
end
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
        
%COMPUTE SURFACE TROPOSPHERIC VALUEs
P0 = intAVG(1) - intAMP(1) * COSphase;
T0 = intAVG(2) - intAMP(2) * COSphase; 
    
if strncmpi(model,'UNB3m',5)
   RH0 = intAVG(3) - intAMP(3) * COSphase;
   %*****TRANSFORM FROM RELATIVE HUMIDITY(RH) TO WVP (IERS Conventions 2003) 
   ES = 0.01 * exp(1.2378847e-5 * (T0 ^ 2) - 1.9121316e-2 * ...
    T0 + 33.93711047 - 6.3431645e3 * (T0 ^ -1));%Saturation Vapour Pressure
   FW = 1.00062 + 3.14e-6 * P0 + 5.6e-7 * ((T0 - 273.15) ^ 2);%Enhancement factor
   E0 = (RH0 / 100) * ES * FW;%Water Vapor Pressure in [mbar or hpa]
   
elseif any([strncmpi(model,'EGNOS',5),strncmpi(model,'MOPS',4)])    
   E0 = intAVG(3) - intAMP(3) * COSphase;       
end  
    
  BETA = intAVG(4) - intAMP(4) * COSphase;
LAMBDA = intAVG(5) - intAMP(5) * COSphase;
   
%***COMPUTE POWER VALUE FOR PRESSURE & WATER VAPOUR
EP = g / (R * BETA);

%*****SCALE SURFACE VALUEs TO USER'S HEIGHT
T = T0 - BETA * hgt; %Surface Temperature at user Location
P = P0 * ( T / T0 ) ^ EP; %Surface Pressure at user Location
E = E0 * ( T / T0 ) ^ ( EP * (LAMBDA+1) );%water vapor pressure at user Location
beta=BETA;
lambda=LAMBDA;

%*****COMPUTE ACCELERATION GRAVITY AT THE ATMOSPHERIC COLUMN
if strncmpi(model,'UNB3m',5)
   dgref = 1 - 2.66e-3*cos(2*latRAD) - 2.8e-07 * hgt;
    gm   = 9.784 * dgref;

elseif any([strncmpi(model,'EGNOS',5),strncmpi(model,'MOPS',4)])
       gm = 9.784;%Acceleration Gravity at the Atmospheric column in m/s^2   
end
%*********COMPUTE MEAN TEMPERATURE OF WATER VAPOR(TM)
TM  = T * (1 - (BETA * R / (gm*(LAMBDA+1))));

%*****************COMPUTE ZENITH DELAYs
if strncmpi(model,'UNB3m',5)
   %**************DEFINE CONSTANTs
   k1 = 77.604;%Refractivity Constant in [K/hPa]
   k2 = 64.79;%Refractivity Constant in [K/hPa]
   Mw = 18.0152;%Molar Mass of Water in [g/mol]
   Md = 28.9644; %Molar Mass of Dry Air  in [g/mol]
  k2p = k2 - k1*(Mw/Md);%Refractivity Constant in [K/hPa] [called k2 prime]
  k3  = 377600;%Refractivity Constant in [ K^2/hPa]       
      
   %***Compute zenith hydrostatic delay(ZHD) using Saastamoinen Model
   ZHD = ((1e-6 * R* k1) / gm)* P;

   %***Compute Zenith Wet Delay(ZWD) using Saastamoinen Model
   ZWD = ((1.0e-6 * (TM*k2p+k3)*R )/((gm*(LAMBDA + 1))- BETA*R))* E/ T;

   %Compute Zenith Total Delay(ZTD)
   ZTD =ZHD + ZWD;

elseif any([strncmpi(model,'EGNOS',5),strncmpi(model,'MOPS',4)])
       %**************DEFINE CONSTANTs
       k1 = 77.604;%Refractivity Constant in [K/hPa]
       k2 = 382000;%Refractivity Constant in [K/hPa] 
       
       %***Compute zenith hydrostatic delay(ZHD)
       ZHD = ((1e-6 * R* k1) / gm)* P;

       %***Compute Zenith Wet Delay(ZWD)
       ZWD = ((1.0e-6 *k2* R )/((gm*(LAMBDA + 1))- BETA*R))* E/ T;

       %Compute Zenith Total Delay(ZTD)
       ZTD =ZHD + ZWD;       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF getZTD_UNB3m_EGNOS_MOPS.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%D_2.SUBROUTINE TO INTERPOLATE UNB3m AVEARAGE & AMPLITUDE VALUEs
function [AVG,AMP] =interpolate_table(UserLat,avgValues,ampValues)
                                                                 
%***************************************************************************
%DESCRIPTION:                                                              * 
%            The function  "interpolate_table" Interpolates and computes   * 
%            averaged and seasonal variation Meteorological parameters...  *   
%            used for tropospheric delay prediction, given receiver...     * 
%            latitude based on the UNB3m, EGNOS & MOPS Models.Parameters   * 
%            above |lat|<=15° and |lat|>=75° are extracted directly while  *   
%            Parameters for latitudes 15°<|lat|<75° are linearly...        * 
%            interpolated between  values of two closest latitudes.        *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%1.       UserLat: Station/User geodetic latitude(s) in degrees(deg)(n x 1)*
%2.     avgValues: Average Meteorological Values                           *
%3.     ampValues: Amplitude/Seasonal Meteorological Values                *
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF interpolate_table.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%      -----------------------------------------------------------
%E.****SUBROUTINE TO COMPUTE ZENITH TROPO DELAYs USING GTrop MODEL
%      -----------------------------------------------------------
function [ZHD,ZWD,ZTD,Tm] = getZTD_GTrop(UTCtime,lat,lon,h,GTropCoeff)

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
           [ZHD(i,1),ZWD(i,1),ZTD(i,1),Tm(i,1)] = GTrop(UTCtime(i,:),lat(i),lon(i),h(i),GTropCoeff); 
      end
      
   else 
       %*****INITIALIZE OUTPUTs 
       ZHD=zeros(nrow_pos,nrow_time);
      
      [ZWD,ZTD,Tm]=deal(ZHD);%create copy of ZHD in ZWD,ZTD,Tm
    
     for i=1:nrow_time %LOOP OVER TIME
         
        for j=1:nrow_pos %LOOP OVER POSITIONS
            
            %Call the "ZenithTropDelay.m" Function
           [ZHD(j,i),ZWD(j,i),ZTD(j,i),Tm(j,i)] = GTrop(UTCtime(i,:),lat(j),lon(j),h(j),GTropCoeff); 
    
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
            [ZHD(j,i),ZWD(j,i),ZTD(j,i),Tm(j,i)] = GTrop(UTCtime(i,:),lat(j),lon(j),h(j),GTropCoeff);     
        end
    end
end
%***********************END OF getZTD_GTrop.m ***********************    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%(B.1)****SUBROUTINE TO COMPUTE ZENITH TROPOSPHERIC DELAY USING GTrop MODEL
function [ZHD,ZWD,ZTD,Tm] = GTrop(UTCtime,lat,lon,h,GTrop_Coeff)

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


%(E.1.1)****SUB-ROUTINE TO COMPUTE TROPOSPHERIC PARAMETER AT A GRID POINT
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
%(E.1.2)********SUB-ROUTINE FOR CHECKING THE EXISTENCE OF GTrop GRID FILE
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
%(E.1.2.1)**MAIN SUB-ROUTINE FOR CHECKING THE EXISTENCE OF GTrop GRID FILE
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


%  --------------------------------------------------------------------
%F.SUBROUTINE TO COMPUTE ZENITH TROPOSPHERIC DELAYS FROM VMF GRID FILES
%  --------------------------------------------------------------------
function [ZHD,ZWD,ZTD,ah,aw,VMF_model] = searchReadVMFgrids(Time,lat,lon,h,VMFgrids)
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
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%Original codes by Daniel Landskron (2017/06/28)                           *
%Modified by:                                                              +
%            OSAH SAMUEL, MSC GEOMATIC ENGINEERING (PhD STUDENT)           +
%            Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com                +
%            Phone:+233(0)246137410 / +233(0)509438484                     +         
%==========================================================================+
%**************************************************************************+
%**************************************************************************+
%CHECK IF VMFgrids IS EMPTY([]) & IMPORT GRIDS USING "SearchVMFgrids.m"
%THIS IS USER FAILS TO PROVIDE VMFgrids IN THE INPUT ARGUMENT

if isempty(VMFgrids)%Determine whether array is empty
    
   %Call the "SearchVMFgrids.m" fxn
   [VMF_grid_found,VMFgrids] = SearchVMFgrids();
   
end

%RETRIEVE VMF GRID & OROGRAPHY FILES FROM STRUCTURE ARRAY[VMFgrids]
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
        
%**********CHECK IF ALL GRID FILES ARE EMPTY([])       
if all([isempty(H0files),isempty(H6files),isempty(H12files),isempty(H18files)])
    
   %RETURN EMPTY([]) MATRICES AS OUTPUT 
   ZHD = []; ZWD = []; ZTD = []; ah= []; aw = [];  VMF_model = [];
           
   return
   
end 
                                                                   
%%%%%%%%%%%%%%%%%%%%%%END OF FILE REFORMATING
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
%*****************CONCANTENATE LAT LON & h ALL TOGETHER 
%1******RECEIVER / STATION POSTION
Rpos  = [lat lon h]; 

%CHECK THE NUMBER OF STATIONS
nrow_Rpos = size(Rpos,1);  
 
%1.GET [YEAR(Yr) MONTH(Mn) DAY]
Yr  = Time(:,1);%get Year
Mn  = Time(:,2);%get Month
Day = Time(:,3);%get Day
   
%ASSIGNMENT
Date = Time;%UTC TIME

%**************LOOK FOR VMF FILE TYPE(VMF1 OR VMF3)

%RETRIEVE STORED VMF GRID TYPE & RESOLUTION('VMF1(2°x2.5°)','VMF3(1°x1°)','VMF3(5°x5°)') from goGPS GUI
VMFgrid_type=getappdata(0,'VMFgrid_type');%VMF GRID TYPE(VMF1 or VMF3)
VMFgrid_res = getappdata(0,'VMFgrid_res');%VMF GRID RESOLUTION especially for VMF3(1°x1° OR 5°x5°)

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-        
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
if ~isempty(h0files)
           
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
         
end %//if ~isempty(h0files)
         
%2.********6 HOUR FILEs
if ~isempty(h6files)
           
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
  
end %//if ~isempty(h6files)
       
%3.********12 HOUR FILEs
if ~isempty(h12files)
           
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
 
end %//if ~isempty(h12files)
       
%4.********18 HOUR FILEs
if ~isempty(h18files)
           
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
 
end %//if ~isempty(h6files)

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
                  
        %LOOP OVER VMFG FILES
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

%(F.1)SUBROUTINE TO READ VMF1 GRID FILE & EXTRACT COEFFICIENTs(ZHD,ZWD,ZTD,ah,aw)
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

%(F.2).SUBROUTINE TO READ VMF3 GRID FILE & EXTRACT COEFFICIENTs(zhd,zwd,ztd,ah,aw)
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

%(F.3).SUBROUTINE TO IMPORT VMF3 OROGRAPHY FILE
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


%              --------------------------------------------------------
%(F.4)*********SUB-ROUTINE FOR CHECKING THE EXISTENCE OF VMF GRID FILES
%              --------------------------------------------------------

function [VMF_grid_found,VMFgrids] = SearchVMFgrids()
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *                                              
%            The "SearchVMFgrids" subroutine searches Vienna Mapping...    +
%             Function Grid (VMF)files in a given Directory/Folder         +
%            as indicated in the sub-routine. Situation where VMFG files   + 
%            are not found in the default directory/folder,it searches     +   
%            recursively through all Sub-Directories/Folders of the...     +  
%            given Directory. VMFG files are extracted by looping ...      + 
%            through all the listed files in the provided folder or ...    +
%            sub-folders.Finally,if VMFG files are still not found in the..+ 
%            default directory/folder and its sub-folders,the search for   +
%            VMFG files is extended to the current directory.              +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%%USAGE:                                                                   +
%      [VMF_grid_found,VMFgrids] = SearchVMFgrids()                        +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%NOTE:                                                                     *
%     THE FF FOLDERS/DIRECTORIES ARE SEARCHED FOR VMF GRID FILES           *
%1.   'VMF grids' folder; main folder for VMF GRIDs                        *
%2.   'VMF files' folder; main folder for VMF GRIDs                        *
%3.   'TropoGRIDS' folder; main folder for Tropospheric & Tm GRIDs         *
%4.   'data' folder; folder for goGPS data                                 *
%5.   'pwd'; Current folder/Directory. In goGPS, it is the 'goGPS' folder  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%--------------------------------------------------------------------------
%OUTPUTs                                                                   *
%1.VMF_grid_found : Flag to indicate presence of VMF grid file.1 to ...*
%                   indicate there is file in the folder & 0 to mean ...   *
%                   absence of file                                        *
%2.VMFgrids     : Structure array of VMF grid files which includes:      *
%1.    h0files  : Zero(0) Hour VMFG File Name(s) eg:'VMFG_20180101.H00'    * 
%2.    h6files  : Six(6) Hour VMFG File  Name(s) eg:'VMFG_20180101.H06'    *  
%3.    h12files : Twelve(12) Hour VMFG File Name(s) eg:'VMFG_20180101.H12' *  
%4.    h18files : Eighteen(18) Hour VMFG File Name(s)eg:'VMFG_20180101.H18'*
%5.    H0files  : Zero(0) Hour VMFG File Path(s)                           * 
%                 eg:'C:\Users\...data\VMFG_20180101.H00'                  * 
%6.    H6files  : Six(6) Hour VMFG File  Path(s)                           *  
%7.    H12files : Twelve(12) Hour VMFG File Path(s)                        *  
%8.    H18files : Eighteen(18) Hour VMFG File Path(s)                      * 
%9.    Orography: Orography file                                           *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
% WRITTEN BY: 
%            OSAH SAMUEL, MSC GEOMATIC ENGINEERING (PhD STUDENT)           +
%            Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com                +
%            Phone:+233(0)246137410 / +233(0)509438484                     +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%**************************************************************************+

%         ------------------------------------------------
%*********SEARCH FOR VMF GRID FILES FROM VARIOUS DIRECTORY
%         ------------------------------------------------

%1.SEARCH VMF GRIDs from 'VMF grids' folder (A SUB-FOLDER]
[VMF_grid_found,VMFgrids] = SearchVMFgrid('../data/TropoGRIDS/VMF files/VMF grids');

if isequal(VMF_grid_found,0)%if VMF grids are not in 'VMF grids' folder, try the 'VMF files' folder
    
   %2.SEARCH VMF GRIDs from 'VMF files' folder
   [VMF_grid_found,VMFgrids] = SearchVMFgrid('../data/TropoGRIDS/VMF files');

   if isequal(VMF_grid_found,0)%if VMF grids are not in 'VMF files' folder, try the 'TropoGRIDS' folder
       
      %3.SEARCH VMF GRIDs from 'TropoGRIDS' folder
      [VMF_grid_found,VMFgrids] = SearchVMFgrid('../data/TropoGRIDS');

      if isequal(VMF_grid_found,0)%if VMF grids are not in 'TropoGRIDS' folder, try the 'data' folder
    
         %4.SEARCH VMF GRIDs from 'data' folder 
         [VMF_grid_found,VMFgrids] = SearchVMFgrid('../data');%searching the 'data' folder
   
         if isequal(VMF_grid_found,0) %if VMF grids are not in the 'data' folder, try the 'current' directory
       
            %5.SEARCH VMF GRIDs from 'current' folder  
            [VMF_grid_found,VMFgrids] = SearchVMFgrid(pwd);%if VMF grids are not in the 'current' directory, then issue a message
      
            if isequal(VMF_grid_found,0)%if VMF grids are not in the 'current' directory, then issue a message
          
                beep %Give a beep sound
                errmsg2{1}=sprintf('No VMF grid file(s) found in goGPS directory.\n');
                errmsg2{2}=sprintf('Please Provide VMF grid file(s) & Try Again.\n');
                errmsg2{3}=sprintf('VMF grid(s) can as well be downloaded from : http://vmf.geo.tuwien.ac.at/trop_products/GRID/\n');
                warndlg(errmsg2,'VMF grid file Error','modal')
                
                 return
                    
            end   
                               
         end 
                                
      end    
      
   end   
   
end  
%===================================END OF Search FOR FOLDER WITH VMF grids
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
function [VMF_grid_found,VMFgrids] = SearchVMFgrid(directory)

%ASSIGNMENT
folder = directory;  %folder/Directory suppose to be containing VMFG FILES in goGPS

%****FIRST CHECK IF FOLDER IS EMPTY OR NOT
SizeListing=size(dir(folder)); %Size of folder content

%********CHECK IF FOLDER IS EMPTY OR NOT
if any([SizeListing(1,1)==2,SizeListing(1,1)==0])% size of 2 means folder is empty & 0 means folder doesn't exist
                  
   VMF_grid_found = 0; %Flag to indicate the absence of VMF grid file in directory
   
   %SET VARIOUS FILE NAMES & PATHS TO EMPTY([])
   h0files  = []; h6files  = []; h12files = []; h18files = [];
   H0files  = []; H6files  = []; H12files = []; H18files = [];
   orogFILE = [];
   orogPATH = [];
   VMFgrids = [];
   
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
          
       %*****GET INDEXEs FOR THE VARIOUS GRID FILES(i.e.Integer Val /[])
       H0Index = regexpi(fileList, '\w*.H00$'); %Index for 0 hour files
       H6Index = regexpi(fileList, '\w*.H06$'); %Index for 6 hour files
       H12Index = regexpi(fileList, '\w*.H12$');%Index for 12 hour files
       H18Index = regexpi(fileList, '\w*.H18$');%Index for 18 hour files
          
       %*********LOOK FOR GRID FILEs WITH .txt EXTENSION
       %NB:some browser during file download will save grid file as text
       %   file with the extension '.txt'.
       H0_1Index = regexpi(fileList, '\w*.H00.txt$'); %Index for 0 hour files
       H6_1Index = regexpi(fileList, '\w*.H06.txt$'); %Index for 6 hour files
       H12_1Index = regexpi(fileList, '\w*.H12.txt$');%Index for 12 hour files
       H18_1Index = regexpi(fileList, '\w*.H18.txt$');%Index for 18 hour files
          
       %************GET THE GRID FILES
       h0Files  = fileList(~cellfun(@isempty, H0Index));%  0 hour files(eg:'VMFG_20180101.H00')
       h6Files  = fileList(~cellfun(@isempty, H6Index)); % 6 hour files(eg:'VMFG_20180101.H06')
       h12Files = fileList(~cellfun(@isempty, H12Index));%12 hour files(eg:'VMFG_20180101.H12')
       h18Files = fileList(~cellfun(@isempty, H18Index));%18 hour files(eg:'VMFG_20180101.H18')
          
       %*********THOSE WITH '.txt' EXTENSION
       h0_1Files  = fileList(~cellfun(@isempty, H0_1Index));%  0 hour files(eg:'VMFG_20180101.H00.txt')
       h6_1Files  = fileList(~cellfun(@isempty, H6_1Index)); % 6 hour files(eg:'VMFG_20180101.H06.txt')
       h12_1Files = fileList(~cellfun(@isempty, H12_1Index));%12 hour files(eg:'VMFG_20180101.H12.txt')
       h18_1Files = fileList(~cellfun(@isempty, H18_1Index));%18 hour files(eg:'VMFG_20180101.H18.txt')
 
       %****************WORK ON VARIOUS GRID FILEs
          
       %IF Any of the Files is Found(i.e. not empty([])),save its  fullfile path 
       %============                       +                        ============
       %1.******0 HOUR GRID FILE
       if ~isempty(h0Files) & ~isempty(h0_1Files)% If h0files & h0_1files are not empty
                                          
           %CONCATENATE h0files & h0_1files(i.e. Fname0)
           h0files = [h0Files;h0_1Files];
             
           %***CONSTRUCT THE FULL PATH       
           H0files = cellfun(@(x) fullfile(folder,x),h0Files,'UniformOutput',false);%Prepend path to files  
             
       elseif ~isempty(h0Files) & isempty(h0_1Files)% If h0files is not empty & h0_1files is empty
                              
              h0files = h0Files;
              %***CONSTRUCT THE FULL PATH        
              H0files = cellfun(@(x) fullfile(folder,x),h0Files,'UniformOutput',false);%Prepend path to files
             
       elseif  isempty(h0Files) & ~isempty(h0_1Files)% If h0files is empty & h0_1files is not empty
                                 
               h0files = h0_1Files;
               %***CONSTRUCT THE FULL  PATH      
               H0files = cellfun(@(x) fullfile(folder,x),h0_1Files,'UniformOutput',false);%Prepend path to files
             
       else   %if they are both empty
              
              h0files = [];
              H0files = [];%Assign empty([]) matrix
              
       end    %if ~isempty(h0files)
           
       %2.******6 HOUR GRID FILE          
       if ~isempty(h6Files) & ~isempty(h6_1Files)% If h6files & h6_1files are not empty
                                        
          %CONCATENATE h6files & h6_1files(i.e. Fname6)
                         
          h6files = [h6Files;h6_1Files];
             
          %********CONSTRUCT THE FULL PATH 
          H6files = cellfun(@(x) fullfile(folder,x),h6Files,'UniformOutput',false);%Prepend path to files 
          
       elseif ~isempty(h6Files) & isempty(h6_1Files)% If h6files is not empty & h6_1files is empty
                             
              h6files = h6Files;
              %***CONSTRUCT THE FULL PATH        
              H6files = cellfun(@(x) fullfile(folder,x),h6Files,'UniformOutput',false);%Prepend path to files
             
       elseif  isempty(h6Files) & ~isempty(h6_1Files)% If h6files is empty & h6_1files is not empty
                             
               h6files = h6_1Files;
               %***CONSTRUCT THE FULL PATH       
               H6files = cellfun(@(x) fullfile(folder,x),h6_1Files,'UniformOutput',false);%Prepend path to files
             
       else %if they are both empty
                          
            h6files = [];
            H6files = [];%Assig empty([]) matrix
              
       end  %if ~isempty(h6files)
          
       %3.******12 HOUR GRID FILE                    
       if ~isempty(h12Files) & ~isempty(h12_1Files)% If h12files & h12_1files are not empty
                                     
          %CONCATENATE h12files & h12_1files(i.e. Fname12)
          h12files = [h12Files;h12_1Files];
             
          %********CONSTRUCT THE FULL PATH
          H12files = cellfun(@(x) fullfile(folder,x),h12Files,'UniformOutput',false);%Prepend path to files  
            
       elseif ~isempty(h12Files) & isempty(h12_1Files)% If h12files is not empty & h12_1files is empty
                               
              h12files = h12Files;
              %***CONSTRUCT THE FULL PATH        
              H12files = cellfun(@(x) fullfile(folder,x),h12Files,'UniformOutput',false);%Prepend path to files
             
       elseif  isempty(h12Files) & ~isempty(h12_1Files)% If h12files is empty & h12_1files is not empty
                             
               h12files = h12_1Files;
               %***CONSTRUCT THE FULL PATH       
               H12files = cellfun(@(x) fullfile(folder,x),h12_1Files,'UniformOutput',false);%Prepend path to files
             
       else %if they are both empty  
                          
            h12files = [];
            H12files = [];%Assig empty([]) matrix
              
       end   %if ~isempty(h12files)
          
       %4.******18 HOUR GRID FILE           
       if ~isempty(h18Files) & ~isempty(h18_1Files)% If h18files & h18_1files are not empty
                          
          %CONCATENATE h18files & h18_1files(i.e. Fname18)
          h18files = [h18Files;h18_1Files];
             
          %********CONSTRUCT THE FULL PATH
          H18files = cellfun(@(x) fullfile(folder,x),h18Files,'UniformOutput',false);%Prepend path to files  
            
       elseif  ~isempty(h18Files) & isempty(h18_1Files)% If h18files is not empty & h18_1files is empty
                               
               h18files = h18Files;
               %***CONSTRUCT THE FULL PATH        
               H18files = cellfun(@(x) fullfile(folder,x),h18Files,'UniformOutput',false);%Prepend path to files
             
       elseif  isempty(h18Files) & ~isempty(h18_1Files)% If h18files is empty & h18_1files is not empty
                            
               h18files = h18_1Files;
               %***CONSTRUCT THE FULL PATH       
               H18files = cellfun(@(x) fullfile(folder,x),h18_1Files,'UniformOutput',false);%Prepend path to files
             
       else  %if they are both empty  
                          
             h18files = [];
             H18files = [];%Assig empty([]) matrix
              
       end   %if ~isempty(h18files)
                      
                      
      %***************GET OROGRAPHY FILE AS WELL 
      if any([~isempty(H0files),~isempty(H6files),~isempty(H12files),~isempty(H18files)]) | all ([~isempty(H0files),~isempty(H6files),~isempty(H12files),~isempty(H18files)])
                     
         %*****CONSTRUCT FILEs PATH
         Filepath=fullfile(strcat(folder,'\',fileList));%Prepend path to files
                     
         u=1;%Loop index for orography file
          
         for r=1:length(fileList) 
                                       
             if strncmpi(fileList(r),'orography_ell',13)| strncmpi(fileList(r),'orography',9)| strncmpi(fileList(r),'orog',4)
                
                orogFILE=fileList(r,1);%Orography file
                                 
                %********CONSTRUCT THE FULL PATH  
                try
                    orogPATH(u,1)=Filepath(r,1);%Orography file  path
                                 
                catch     
                      orogPATH = cellfun(@(x) fullfile(folder,x),orogFILE,'UniformOutput',false);%Prepend path to files  
                end     
                 
                u=u+1;%Update index
                                                  
             end  %//if strncmpi(fileList(r),'orography_ell',13)|   
             
         end  %//for r=1:length(fileList)     
                      
      else 
           %SET OROGRAPHY FILE NAME & PATH TO EMPTY([])
           orogFILE = [];
           orogPATH = [];
              
      end   %if any([~isempty(H0files),~isempty(H6files),~isempty(H12files)
                                                       
    else
        %SET VMF GRID FILE & OROGRAPHY FILE NAMES & PATHS TO EMPTY([])  
        h0files  = []; h6files  = []; h12files = []; h18files = [];
        H0files  = []; H6files  = []; H12files = []; H18files = [];
        orogFILE = [];
        orogPATH = [];
           
    end  %if ~isempty(fileList)
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    %IF NO VMFG FILES EXIST IN THE CURRENT FOLDER, CHECK SUB-FOLDERS/              
    %SUB-DIRECTORIES FOR VMF GRID FILEs
    %                              +                        ==========
     if isempty(fileList)
              
        checkSubfolders=1; %Flag to search Sub-Directories/folders
           
     elseif all([isempty(h0Files),isempty(h6Files),isempty(h12Files),isempty(h18Files)])
                 
            checkSubfolders=1;%Flag to search Sub-Directories/folders 
                         
     else
          if any([~isempty(h0files),~isempty(h6files),~isempty(h12files),~isempty(h18files)])
                  
             checkSubfolders=0;%Flag not to search Sub-Directories/folders
             
             VMF_grid_found = 1; %Flag to indicate the presence of VMF grid file in directory
             
          end
                  
                                                                  
     end    %\\if isempty(fileList)
                
     
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
                   
              end  %try
                       
              Filelists{p,1}= fileLists;
           
              Filepaths{p,1}=Filepath;
           
              p=p+1;%Update index
                       
              %COMBINE ALL RINEX FILES FROM VARIOUS SUB-FOLDERS  
              filelists = vertcat(Filelists{:});     
              filepaths=vertcat(Filepaths{:}) ;                     
                      
              if ~isempty(fileLists)   % If fileLists is not empty
            
                 j=1;k=1;l=1;m=1; %Loop index for the VMF grid file types
                 jj=1;kk=1;ll=1;mm=1; %Loop index for the VMF grid file types
                 u=1;%Loop index for orography file
                        
                 for i = 1:length(filelists)%****LOOP OVER THE FILES & CHECK FOR VMF GRID FILES
                                            %Identify the needed files and extract them separately
                                            %into different file names(H0files,H6files,H12files & H18files )
                                            %Note that dir also lists the directories,so you have to check for them. 
                                                         
                      %***GET NAME & PATH OF FILES FROM THE LIST OF FILES FOR EACH LOOP               
                           
                      FileName=filelists{i,1};
                      filepath=filepaths{i,1};
                                          
                      %****MATCH THE END OF EACH FILE IN THE 'FileName'
                      %For 0 Hour file look for 'H00'
                      %For 6 Hour file look for 'H06' ...
                      %For 12 Hour file look for 'H12' And...
                      %For 18 Hour file look for 'H18' 
                 
                      %1.******0 HOUR GRID FILE 
                      if regexpi(FileName,'\w*.H00$') 
  
                         %When Tested & the file is a 0h file,save...
                         %its fullfile path into the H0files cell array
                         %============            +        ============
                         h0Files{j,1}=FileName;
                         H0filess{j,1}=filepath;%#ok<*AGROW> %0h file in cells                             
    
                         j=j+1;%Update index
                                    
                      %2.******6 HOUR GRID FILE 
                   
                      elseif  regexpi(FileName,'\w*.H06$' ) %Else test the file if it is a 6 hour file
      
                              %When Tested & the file is a 6h file,save...
                              %its fullfile path into the H6files cell array
                              %============            +        ============
                              h6Files{k,1}=FileName;
                              H6filess{k,1}=filepath;%6 Hour file in cells
                                   
                              k=k+1; %update index
                 
                      %3.******12 HOUR GRID FILE          
                      elseif  regexpi(FileName,'\w*.H12$')
      
                              %When Tested & the file is a 12h file,save...
                              %its fullfile path into the H12files cell array
                              %============            +        ============ 
                              h12Files{l,1}=FileName;
                              H12filess{l,1}=filepath;%12 Hour file in cells
                                 
                              l=l+1;  %update index  
                                 
                      %4.******18 HOUR GRID FILE          
                      elseif regexpi(FileName,'\w*.H18$')
      
                             %When Tested & the file is a 18h file,save...
                             %its fullfile path into the H18files cell array
                             %============            +        ============  
                             h18Files{m,1}=FileName;
                             H18filess{m,1}=filepath;%18 Hour file in cells
                                 
                             m=m+1;  %update index 
                                 
                      %*********THOSE WITH '.txt' EXTENSION       
                      %1.******0 HOUR GRID FILE 
                      elseif  regexpi(FileName,'\w*.H00.txt$') 
                              h0_1Files{jj,1}=FileName;                         
                              H0_1filess{jj,1}=filepath;%0h file in cells                             
    
                              jj=jj+1;%Update index
                                  
                      %2.******6 HOUR GRID FILE 
                      elseif  regexpi(FileName,'\w*.H06.txt$') 
                              h6_1Files{kk,1}=FileName;
                              H6_1filess{kk,1}=filepath;%6h file in cells                             
    
                              kk=kk+1;%Update index  
                                   
                      %3.******12 HOUR GRID FILE 
                      elseif  regexpi(FileName,'\w*.H12.txt$') 
                              h12_1Files{ll,1}=FileName;
                              H12_1filess{ll,1}=filepath;%12h file in cells                             
    
                              ll=ll+1;%Update index          
                                        
                      %4.******18 HOUR GRID FILE 
                      elseif  regexpi(FileName,'\w*.H18.txt$') 
                              h18_1Files{mm,1}=FileName;
                              H18_1filess{mm,1}=filepath;%18h file in cells                             
    
                              mm=mm+1;%Update index 
                              
                                                         
                      end   %//if regexpi(FileName,'\w*.H00$')
                           
                      %IMPORT OROGRAPHY FILE AS WELL
                                                                                                    
                      if strncmpi(FileName,'orography_ell',13)| strncmpi(FileName,'orography',9)| strncmpi(FileName,'orog',4)
                
                         orogFILE{u,1}=FileName;%#ok<*NASGU> %Orography file 
                                  
                         %********CONSTRUCT THE FULL PATH                                           
                         orogPATH{u,1}=filepath;%Orography file  path
                                                          
                         u=u+1;%Update index
                                                                                                
                      end       
                      
                 end %//for i = 1:length(fileLists) 
                                                                                                               
              end    %//if ~isempty(fileLists)
                      
            else                          
                 %WRITE ERROR MESSAGE
                 erMsg{q,1}=sprintf('The Folder %s and  It''s Sub-Folder, %s , are empty / Contain No VMF Grid files .\n',folder,subDirs{iDir});
                         
                 q=q+1; %update index
                        
            end   %if Sizesublist(1,1)~=2
                   
        end    %for iDir = find(validIndex)
        
      if all([~exist('H0Files','var'),~exist('H6Files','var'),~exist('H12Files','var'),~exist('H18Files','var'),...
              ~exist('h0_1Files','var'),~exist('h6_1Files','var'),~exist('h12_1Files','var'),~exist('h18_1Files','var')]) 
         
         VMF_grid_found = 0;%Flag to indicate the absence of VMF grid file in directory
         
         VMFgrids = [];
         
         return
      
      else 
          VMF_grid_found = 1; %Flag to indicate the presence of VMF grid file in directory
          
          %*********************COMBINE VMF GRID FILEs
          %********0 HOUR FILEs
          if exist('H0filess','var') & exist('H0_1filess','var')
             h0files = [h0Files;h0_1Files];
             H0files = [H0filess;H0_1filess];
               
          elseif exist('H0filess','var') & ~exist('H0_1filess','var')
                 h0files = h0Files;
                 H0files = H0filess;
                  
          elseif ~exist('H0filess','var') & exist('H0_1filess','var')
                 h0files = h0_1Files;
                 H0files = H0_1filess;
                                         
          else  
              h0files = [];
              H0files = [];
                                                                                                
          end %\\if exist('H0filess','var')
                        
          %********6 HOUR FILEs
          if exist('H6filess','var') & exist('H6_1filess','var')
             h6files = [h6Files;h6_1Files];                                         
             H6files = [H6filess;H6_1filess];
                 
          elseif exist('H6filess','var') & ~exist('H6_1filess','var')
                 h6files = h6Files;  
                 H6files = H6filess;
                 
          elseif ~exist('H6filess','var') & exist('H6_1filess','var')
                 h6files = h6_1Files;  
                 H6files = H6_1filess;
          else     
              h6files = [];  
              H6files = [];                   
                                                                                     
          end   %\\if exist('H6filess','var')
            
          %********12 HOUR FILEs
          if exist('H12filess','var') & exist('H12_1filess','var')
             h12files = [h12Files;h12_1Files];                                  
             H12files = [H12filess;H12_1filess];
                
          elseif exist('H12filess','var') & ~exist('H12_1filess','var')
                 h12files = h12Files;                              
                 H12files = H12filess;
                                      
          elseif ~exist('H12filess','var') & exist('H12_1filess','var') 
                 h12files = h12_1Files;
                 H12files = H12_1filess;  
                  
          else   
              h12files = [];
              H12files = [];
               
          end %//if exist('H12filess','var') & exist('H12_1filess','var')
                 
          %********18 HOUR FILEs
          if exist('H18filess','var') & exist('H18_1filess','var')
             h18files = [h18Files;h18_1Files];                                     
             H18files = [H18filess;H18_1filess];
                
          elseif exist('H18filess','var') & ~exist('H18_1filess','var')
                 h18files = h18Files;                            
                 H18files = H18filess;
                                      
          elseif ~exist('H18filess','var') & exist('H18_1filess','var') 
                 h18files = h18_1Files;
                 H18files = H18_1filess;   
                                      
          else  
              h18files = [];
              H18files = []; 
                                              
          end  %\\if exist('H18filess','var')
          
         
      end %//if all([~exist('H0Files','var'),~exist('H6Files','var'),~exist('H12Files','var'),~exist('H18Files','var'),...
          %          ~exist('h0_1Files','var'),~exist('h6_1Files','var'),~exist('h12_1Files','var'),~exist('h18_1Files','var')]) 
              
           
     end %//if checkSubfolders==1
     
end %//if exist('listFiles','var')

%********SORT VMF grid FILES & REMOVE DUPLICATE FILES
        
%***********0 HOUR FILES(0h)
if ~isempty(H0files)
   
   [h0files,i_h0] = unique(h0files);%GET 0 HOUR FILES with no repetitions & INDEX(i_h0)
   H0files =(H0files(i_h0));%FILES WITH PATH          
end 
                
%***********6 HOUR FILES(6h)
if ~isempty(H6files)
   [h6files,i_h6] = unique(h6files);%GET 6 HOUR FILES with no repetitions & INDEX(i_h6)
    H6files =(H6files(i_h6));%FILES WITH PATH           
end 
        
%***********12 HOUR FILES(12h)
if ~isempty(H12files)
   [h12files,i_h12] = unique(h12files);%GET 12 HOUR FILES with no repetitions & INDEX(i_h12)
   H12files =(H12files(i_h12));%FILES WITH PATH             
end 
        
%***********18 HOUR FILES(18h)
if ~isempty(H18files)
   [h18files,i_h18] = unique(h18files);%GET 18 HOUR FILES with no repetitions & INDEX(i_h18)
   H18files =(H18files(i_h18));%FILES WITH PATH
end        
                      
%***********OROGRAPHY FILES(orography_ell)
if exist('orogFILE','var') 
    
   if ~isempty(orogPATH)
       
      [orofile,i_oro] = unique(orogFILE);%GET orography file with no repetitions & INDEX(i_oro)
      oropath =(orogPATH(i_oro));%FILES WITH PATH 
       
   else
       orofile = []; %SET orography FILE TO EMPTY([]) IF IT DOES NOT EXIST IN FOLDER   
       oropath = [];%SET orography FILE PATH TO EMPTY([])
   end
   
end     

%***********CHECK IF ALL GRID FILES ARE EMPTY([]) 
if all([isempty(H0files),isempty(H6files),isempty(H12files),isempty(H18files)])
   
   VMFgrids       = [];
   VMF_grid_found = 0;%Flag to indicate the absence of VMF grid file in directory
   
else
    
     %CREATE A STRUCTURE ARRAY(STRUCT) OF ALL THE GRID FILES 
     %***00 HOUR
     h00.name = h0files;%FILE NAME[e.g.:'VMFG_20180101.H00']
     h00.path = H0files;%FILE PATH['../data/TropoGRIDS/VMF files/VMF grids/'VMFG_20180101.H00'

     %****6 HOUR
     h06.name=h6files;%FILE NAME[e.g.:'VMFG_20180101.H06']
     h06.path=H6files;%FILE PATH['../data/TropoGRIDS/VMF files/VMF grids/'VMFG_20180101.H06'

     %****12 HOUR
     h12.name = h12files;%FILE NAME[e.g.:'VMFG_20180101.H12']
     h12.path = H12files;%FILE PATH['../data/TropoGRIDS/VMF files/VMF grids/'VMFG_20180101.H12'

     %****18 HOUR
     h18.name = h18files;%FILE NAME[e.g.:'VMFG_20180101.H18']
     h18.path = H18files;%FILE PATH['../data/TropoGRIDS/VMF files/VMF grids/'VMFG_20180101.H18'

     %OROGRAPHY FILE
     Oro.name = orofile;%FILE NAME[e.g.:'orography_ell']
     Oro.path = oropath;%FILE PATH['../data/TropoGRIDS/VMF files/VMF grids/'orography_ell'

     %COMBINE ALL INTO ONE STRUCTURE(STRUCT)-VMFgrids
     VMFgrids.H00 = h00;
     VMFgrids.H06 = h06;
     VMFgrids.H12 = h12;
     VMFgrids.H18 = h18;
     VMFgrids.Orography = Oro;
end 

%=========================================END OF SearchVMFgrid.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

%  ----------------------------------------------------------------------
%G.SUBROUTINE TO COMPUTE JULIAN,MODIFIED DAY(JD,MJD) &  DAY OF YEAR(DoY)
%  ----------------------------------------------------------------------
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
 
%G_1.*********SUBROUTINE TO COMPUTE DAY OF YEAR(DoY)
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

%G_2*********SUBROUTINE TO CONVERT TWO DIGITs TO YEAR 4 DIGITS
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