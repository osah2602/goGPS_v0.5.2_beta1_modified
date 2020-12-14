%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *
%            "computePI" is a subroutine aimed at Computing PI(pie)- ...   *
%            dimensionless water conversion factor from Weighted Mean ...  *
%            Temperature(Tm).PI needed to convert Zenith Wet Delays(ZWD) to*
%            Precipitable Water Vapor(PWV) from Global Navigation          *
%            Satellite System (GNSS)signal. However,In ground-based GNSS   * 
%            meteorology, Tm is a key parameter to calculate the conversion 
%            factor(PI) that can convert the ZWD to PWV. In this regard in *
%            this subroutine, Tm is computed from various models for the   *
%            computation of PI.                                            *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%USAGE:                                                                    *
%      [PI,Tm] = computePI(Tm_model,Time,ReceiverPos,Temp,gridV_TVGG,...   *
%                          gridV_GTVR,gridV_GTrop,GPT_grid,gridRES,Timevar)*
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%1.    Tm_model : Tm Models taken from goGPS GUI FIGURE Tm popup menu      *
%2.    Time     : Time of observation                                      *
% Time format   : [year month day hour minute seconds] OR [year month day] *
%3. ReceiverPos : Receiver position of the form Latitude,Longitude & Height*
%                 [lat lon h] where latitude & longitude coordinates can   *
%                 either be in DMS or Decimal Degress(Deg)& Height in [m]. *
%4.    Temp     : Surface Temperature in degree Celcius(°C)                *
%5. gridV_TVGG  : TVGG-Tm model grid values extracted from ...             *
%                 Coeff_TVGG_ERAI.mat file                                 *
%6. gridV_GTVR  : GTVR-Tm model grid values extracted from GTVR-Tm.txt file*
%7.gridV_GTrop  : GTrop model grid values extracted from GTropCoefficient.mat
%8.    GPT_grid : grid values in cells extracted from GPT2w or GPT3 1°x1°  *  
%                 OR 5°x5° external grid file. The subroutine "readGPTgrid"*
%                 should be ran first to provide the grid inputfile(GPT_grid)
%e.g.: GPT_grid = readGPTgrid('gpt2_5w.grd','GPT2w',5); for 5°x5° grid     +
%           OR                                                             +
%      GPT_grid = readGPTgrid('gpt2_5w.mat','GPT2w',5); For 5°x5° grid     +
%           OR                                                             +
%      GPT_grid = readGPTgrid('gpt3_5.grd','GPT3',5);  For 5°x5° grid      +
%           OR                                                             +
%      GPT_grid = readGPTgrid('gpt3_5.mat','GPT3',5);  For 5°x5° grid      +

%9.    gridRES  : Grid resolution. 1 for 1°x1° and 5 for 5°x5° resolutions *
%                 respectively.                                            *
%10.    Timevar  : case 1: no time variation but static quantities         *
%                 case 0: with time variation (annual and semiannual terms)*
%NOTE:                                                                     *
%     Generally the sub-routine accepts 10 sets of inputs.However for Tm   *
%     models that do not require a particular input, that particular input *
%     should be assigned with empty([]) matrix.For example if Tm model is  *
%     Bevis i.e. Tm_model = 'Bevis et al (1992) model',Time, ReceiverPos...*
%     ,gridV_TVGG,gridV_GTVR,gridV_GTrop, GPT_grid,gridRES,Timevar can be  *                                              
%     made empty([]) & the sub-routine will still/run work.That is :  
%     Time = []; ReceiverPos = []; gridV_TVGG = []; gridV_GTVR = [];       *
%     gridV_GTrop = []; GPT_grid = []; gridRES = [];Timevar = []           *

%     AND USAGE WILL BE : [PI,Tm]=computePI('Bevis et al (1992) model',... *
%                                             [],[],Temp,[],[],[],[],[],[])*
%                                                              
%    OR: [PI,Tm] = computePI('Bevis et al (1992) model',[],[],Temp)        *

%NB : GPT_grid,gridRES,Timevar are need for GPT MODELS(SUCH AS GPT2w,GPT3) *
%     gridV_TVGG, gridV_GTVR & gridV_GTrop are needed for TVGG-Tm, GTVR-Tm *
%     GTrop models respectively.

%--------------------------------------------------------------------------
%OUTPUT:                                                                   *
%       PI = Water Conversion factor                                       *
%       Tm = Weighted Mean Temperature                                     *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================
%REFERENCE:                                                                *
%          M. Bevis, S. Businger, S. Chiswell, T.A. Herring, R.A. Anthes,  *
%          C. Rocken, et al.GPS meteorology: mapping zenith wet delays onto*
%          precipitable water J Appl Meteorol, 33 (3) (1994), pp. 379-386  *
% =========================================================================
%  Codes Written by:                                                       *
%                  OSAH SAMUEL, Msc Geomatic Engineering ,2016             *    
%                  Email: osahsamuel@yahoo.ca                              *
%                  Tel:+233 (0)246137410/+233 (0)509438484                 * 
%==========================================================================
%***************************************************************************
%***************************************************************************  
%in order to save a huge time amount.
function [PI,Tm] = computePI(Tm_model,Time,ReceiverPos,Temp,gridV_TVGG,...
                             gridV_GTVR,gridV_GTrop,GPT_grid,gridRES,Timevar)

%**********CHECK & REFORMAT INPUTs DATA
switch nargin
    
   case {10,9,8,7,6,5,4,3} %Various inputs format 
       
       %(1)*****************CHECK TIME INPUT
        %------------------------------------------------------------------
        %SOME MODELS REQUIRE TIME INPUT TO PERFORM / FUNCTION.IF TIME INPUT
        %IS EMPTY([])& ANY OF THESE MODELS IS PROVIDED,TERMINATE PROCESS
        %------------------------------------------------------------------
        if isempty(Time)
            
            UTCtime = []; %Assign empty([]) matrix to UTCtime
            
           %*********CHECK TIME DEPENDENT MODELS
           if any([strncmpi(Tm_model,'GTm-I model [Yao et al 2012]',28),strncmpi(Tm_model,'GTm-II model [Yao et al 2013]',29),strncmpi(Tm_model,'GTm-III model [Yao et al 2014]',30),...
                   strncmpi(Tm_model,'TVGG-Tm model [Jiang et al 2018]',32),strncmpi(Tm_model,'GTVR-Tm model [Yao et al 2018]',30),strncmpi(Tm_model,'UNB3m model [Leandro et al 2006]',32),...
                   strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35),strncmpi(Tm_model,'Simplified PI  model [Manandhar et al 2017]',43),...
                   strncmpi(Tm_model,'GTrop model [Sun et al 2019]',28)])
                            
              %ISSUE ERROR MESSAGE  
              beep%Give a beep sound 
              errmsg0{1}=sprintf('Observation / Reception Time is not provided i.e it''s Empty ([ ]).\n');
              errmsg0{2}='Please provide Observation / reception Time & Try again.';
              errordlg(errmsg0,'Time Input Error','modal')  
       
              PI = [];
              return
       
           end
       
        else
            %*****CHECK TIME FORMAT
            ncol=size(Time,2);% finding Number of columns 
           
            switch ncol              
                case  {6,5,4,3}
                  
                   if (any(ncol==[6,5,4,3]))
                    
                     Date = Time(:,1:3);%Assigning column 1 to 3 of time to date                     
                     %***Date
                     Yr  = Date(:,1);%Assigning 1st column to years
                     Mn  = Date(:,2);%Assigning 2ND column to months
                     Day = Date(:,3);%Assigning 3RD column to days
               
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
                         
                             Hr   = zeros(size(Year,1),1);%Assigning zeros to Hour
                             Min  = zeros(size(Year,1),1);%Assigning zeros to Minute
                             Secs = zeros(size(Year,1),1);%Assigning zeros to Seconds
                       
                     end  %if ncol==6 
                     
                     %*****CHANGE SECs, MIN, HOUR WHOSE SEC == 60
                     Min(Secs==60) = Min(Secs==60)+1;
                     Secs(Secs==60) = 0;
                     Hr(Min==60) = Hr(Min==60)+1;
                     Min(Min==60)=0;
                     
                    %CONCATENATE TIME ELEMENTS
                    UTCtime = [Yr Mn Day Hr Min Secs];%UTC time
                    
                   end   %if (any(ncol==[5,4,3,2]))
                   
            end  %switch ncol
            
        end %//if isempty(Time)
        
      %(2)*****************CHECK RECEIVER/SITE POSITION/COORDINATE  INPUT
      %-------------------------------------------------------------------
      %SOME MODELS REQUIRE RECEIVER/SITE COORDs TO PERFORM / FUNCTION.IF 
      %RECEIVER POSITION INPUT IS EMPTY([])& ANY OF THESE MODELS IS 
      %PROVIDED,TERMINATE PROCESS
      %--------------------------------------------------------------------
      if isempty(ReceiverPos)
          
         lat = [];%Assigning empty matrix ([]) to latitude coord(lat)
         lon = [];%Assigning empty matrix ([]) to longitude coord(lon)
           h = [];%Assigning empty matrix ([]) to station/receiver height(h))
          
         %*********CHECK POSITION DEPENDENT MODELS
         if any([strncmpi(Tm_model,'GTm-I model [Yao et al 2012]',28),strncmpi(Tm_model,'GTm-II model [Yao et al 2013]',29),strncmpi(Tm_model,'GTm-III model [Yao et al 2014]',30),...
                 strncmpi(Tm_model,'TVGG-Tm model [Jiang et al 2018]',32),strncmpi(Tm_model,'GTVR-Tm model [Yao et al 2018]',30),strncmpi(Tm_model,'UNB3m model [Leandro et al 2006]',32),...
                 strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35),strncmpi(Tm_model,'Simplified PI  model [Manandhar et al 2017]',43),...
                 strncmpi(Tm_model,'GTrop model [Sun et al 2019]',28)])
                            
        
            %ISSUE ERROR MESSAGE 
            beep%Give a beep sound 
            errmsg1{1}=sprintf('Reciever / Site Position is not provided i.e it''s Empty ([ ]).\n');
            errmsg1{2}='Please provide Receiver / Site position & Try again.';
            errordlg(errmsg1,'Reciever Position(s) Input Error','modal')  
           
              %RETURN EMPTY([]) MATRIX
              PI =[];
              
              return
          
         end
         
      else
          %CHECK COORD ENTRY FORMAT
          nCOL_pos=size(ReceiverPos,2); %Get number of Columns entry
                   
          switch nCOL_pos
              case 7
                    lat = dms2degrees(ReceiverPos(:,1:3));%Latitude in DMS (columns 1-3) converted to degrees 
                    lon = dms2degrees(ReceiverPos(:,4:6));%Longitude in DMS (columns 4-6) converted to degrees   
                      h = ReceiverPos(:,end);%Assigning end (7th) column to heights 
                   
              case 6   
                    lat = dms2degrees(ReceiverPos(:,1:3));%Latitude in DMS (columns 1-3) converted to degrees 
                    lon = dms2degrees(ReceiverPos(:,4:6));%Longitude in DMS (columns 4-6) converted to degrees   
                     h = zeros(size(ReceiverPos,1),1);%Assigning zeros to  heights
                             
              case 5
                    lat = dm2degrees(ReceiverPos(:,1:2));%Latitude in DM (columns 1-2) converted to degrees 
                    lon = dm2degrees(ReceiverPos(:,3:4));%Longitude in DM (columns 3-4) converted to degrees   
                      h = ReceiverPos(:,end);%Assigning end (5th) column to heights
                             
              case 4
                    lat = dm2degrees(ReceiverPos(:,1:2));%Latitude in DM (columns 1-2) converted to degrees 
                    lon = dm2degrees(ReceiverPos(:,3:4));%Longitude in DM (columns 3-4) converted to degrees   
                      h = zeros(size(ReceiverPos,1),1);%Assigning zeros to  heights  
                             
              case 3
                    lat = ReceiverPos(:,1);%Latitude in degrees(columns 1) 
                    lon = ReceiverPos(:,2);%Longitude in degrees(columns 2)
                      h = ReceiverPos(:,end);%Assigning end (3rd) column to heights       
    
              case 2
                   lat = ReceiverPos(:,1);%Latitude in degrees(columns 1) 
                   lon = ReceiverPos(:,2);%Longitude in degrees(columns 2)
                     h = zeros(size(ReceiverPos,1),1);%Assigning zeros to  heights 
            
          end   %switch nCOL_lat
          
      end %//if isempty(ReceiverPos)
             
      %(3)*****************CHECK TEMPERATURE(TEMP) INPUT
      %------------------------------------------------------------------
      %SOME MODELS REQUIRE TEMP INPUT TO PERFORM / FUNCTION.IF TEMP INPUT
      %IS EMPTY([])& ANY OF THESE MODELS IS PROVIDED,TERMINATE PROCESS
      %------------------------------------------------------------------
      if (any(nargin==[10,9,8,7,6,5,4]))
          
          %ASSIGNMENT
          T = Temp ;
          
          if isempty(Temp)
              
             %*********CHECK TEMPERATURE(TEMP) DEPENDENT MODELS
             if any([strncmpi(Tm_model,'Bevis et al (1992) model',24),strncmpi(Tm_model,'Bevis et al (1995) model',24),strncmpi(Tm_model,'Mendes et al (2000) model',25),...
                     strncmpi(Tm_model,'Schueler et al (2001) model',27),strncmpi(Tm_model,'Yao et al (2014) model',22),strncmpi(Tm_model,'TVGG-Tm model [Jiang et al 2018]',32)])
                            
                %ISSUE ERROR MESSAGE 
                beep%Give a beep sound 
                errmsg2{1}=sprintf('Surface Temperature is not provided i.e it''s Empty ([ ]).\n');
                errmsg2{2}='Please provide Surface Temperature & Try again.';
                errordlg(errmsg2,'Temperature Input Error','modal')  
           
                %RETURN EMPTY([]) MATRIX
                PI =[];
              
                return
          
             end  
          
          else   
               %CONVERT TEMPERATURE(T) IN Celcius(°C) to Kelvin(K) 
               T = T + 273.15; %Temperature in kelvin(K)
          
          end %//if isempty(Temp)
          
         %*****************WHEN 10,9 or 8 INPUTS ARE PROVIDED
         %NOTE: SOME MODELS DO NOT REQUIRE GRID FILES & RESOLUTIONS TO PERFORM
         %IF THOSE MODELS ARE NOT PROVIDED, AN ERROR MESSAGE IS ISSUED
         %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=* 
         if (any(nargin==[10,9,8])) %FOR 10,9,8 INPUTS
              
            %(4)*****************CHECK  GPT_grid INPUT
            %--------------------------------------------------------------
            %SOME MODELS REQUIRE  GPT_grid INPUT TO PERFORM / FUNCTION.IF  
            %GPT_grid INPUT IS EMPTY([])& ANY OF THESE MODELS IS PROVIDED,
            %TERMINATE PROCESS
            %-------------------------------------------------------------
             if isempty(GPT_grid) %CHECK if GPT_grid is empty([])
                 
                 %*********CHECK GPT GRID DEPENDENT MODELS 
                 if any([strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35)])
                    
                     %ISSUE ERROR MESSAGE 
                     beep%Give a beep sound
                     errmsg3{1}=sprintf('No input for GPT grid values.\n');
                     errmsg3{2}='Please provide GPT grid values & Try again.';
                     errordlg(errmsg3,'Grid values Input Error','modal') 
                     
                     %RETURN EMPTY([]) MATRIX
                     PI =[];
              
                    return
                    
                 end
                 
             end %//if isempty( GPT_grid)
             
             %*****************WHEN 10 or 9 INPUTS ARE PROVIDED
             if (any(nargin==[10,9])) %FOR 10 or 9 INPUTS
                 
               %(5)*****************CHECK  GRID RESOLUTION(gridRES) INPUT
               %-----------------------------------------------------------
               %SOME MODELS REQUIRE  gridRES INPUT TO PERFORM / FUNCTION.IF  
               %gridRES INPUT IS EMPTY([])& ANY OF THESE MODELS IS PROVIDED,
               %A DEFAULT IS A SIGNED
               %-----------------------------------------------------------
               if isempty(gridRES) %CHECK if gridRES is empty([])
                 
                  %*********CHECK gridRES DEPENDENT MODELS 
                  if any([strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35)])
              
                     if ~isempty(GPT_grid) %if GPT_grid is empty([])
                         
                        gridRES = 5; %5° grid resolution
                        
                     else
                         %ISSUE ERROR MESSAGE  
                         beep%Give a beep sound 
                         errmsg4{1}=sprintf('No input for GPT grid resolution.\n');
                         errmsg4{2}='Please provide GPT grid resolution & Try again.';
                         errordlg(errmsg4,'GPT Grid resolution Input Error','modal') 
                         
                         %RETURN EMPTY([]) MATRIX
                         PI = [];
                         
                         return
                         
                     end % if ~isempty(GPT_grid) 
                    
                  end   
                 
               end  %//if isempty(gridRES)
               
               if nargin == 10
               
                  %(6)*****************CHECK  TIME VARIATION(Timevar) INPUT
                  %--------------------------------------------------------
                  %SOME MODELS REQUIRE  Timevar INPUT TO PERFORM / FUNCTION.IF  
                  %Timevar INPUT IS EMPTY([])& ANY OF THESE MODELS IS PROVIDED,
                  %A DEFAULT IS A SIGNED
                  %--------------------------------------------------------
                  if isempty(Timevar) %CHECK if gridRES is empty([])
                 
                     %*********CHECK Timevar DEPENDENT MODELS 
                     if any([strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35)])
              
                        if ~isempty(GPT_grid) %if GPT_grid is empty([])
                         
                           Timevar = 0; %to apply time variation(annual and semiannual terms)
                        
                        end  %if ~isempty(GPT_grid) 
                    
                     end    
                 
                  end   %//if isempty(Timevar)
               
                  
              %*****************WHEN 8 INPUTS ARE PROVIDED    
               elseif nargin == 9
                      
                      %*********CHECK  Timevar DEPENDENT MODELS 
                      if any([strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35)])
              
                         if ~isempty(GPT_grid) %if GPT_grid is empty([])
                         
                             Timevar = 0; %to apply time variation(annual and semiannual terms)
                             
                         else
                             Timevar = [];%Assigning empty matrix ([]) to Timevar
                         end
                         
                      end
                      
               end  %//if nargin == 10
               
             elseif nargin == 8 %For 8 INPUTS(I.E. gridRES & Timevar ARE NOT PROVIDED)  
                 
                    %*********CHECK  gridRES & Timevar DEPENDENT MODELS 
                    if any([strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35)])
              
                      if ~isempty(GPT_grid) %if GPT_grid is empty([])
                         
                         gridRES = 5; %5° grid resolution
                         Timevar = 0; %to apply time variation(annual and semiannual terms)
                         
                      else
                          gridRES = [];%Assigning empty matrix ([]) to gridRES 
                          Timevar = [];%Assigning empty matrix ([]) to Timevar
                          
                      end
             
                    end
                    
             end  %//if (any(nargin==[10,9]))
               
                
         %*****************WHEN 7,6 or 5 INPUTS ARE PROVIDED    
         elseif  (any(nargin==[7,6,5]))
             
                 %NOTE: SOME MODELS REQUIRE GRID FILES, RESOLUTIONS &  
                 %TIME VARIATION(Timevar)TO PERFORM. IF THOSE MODELS ARE  
                 %PROVIDED, AN ERROR MESSAGE IS ISSUED
                 %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*  
                 %*********CHECK  GRID FILE & gridRES DEPENDENT MODELS 
                 if any([strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35)])
                  
                    %ISSUE ERROR MESSAGE  
                    beep%Give a beep sound 
                    errmsg5{1}=sprintf('No input for GPT grid Values, grid resolution & time variation.\n');
                    errmsg5{2}='Please provide GPT grid Values, grid resolution & time variation & Try again.';
                    errordlg(errmsg5,'GPT grid Values & resolution Input Error','modal') 
                    
                    %RETURN EMPTY([]) MATRIX
                    PI = [];
                         
                    return
                    
                 else
                      GPT_grid = [];%Assigning empty matrix ([]) to GPT_grid
                      gridRES  = [];%Assigning empty matrix ([]) to gridRES
                      Timevar  = [];%Assigning empty matrix ([]) to Timevar
                      
                 end
                 
                 if (any(nargin==[6,5,4]))
                  
                    gridV_GTrop = [];%Assigning empty matrix ([]) to gridV_GTrop
                 
                    if (any(nargin==[5,4]))
                     
                       gridV_GTVR = [];%Assigning empty matrix ([]) to gridV_GTVR
                       
                       if nargin == 4
                     
                          gridV_TVGG = [];%Assigning empty matrix ([]) to gridV_TVGG
                          
                       end %if nargin == 4 
                 
                    end  %//if (any(nargin==[5,4]))
                 
                 end %if (any(nargin==[6,5,4]))
                 
         end   %//if (any(nargin==[10,9,8]))
       
      %*****************WHEN 3 INPUTS ARE PROVIDED   
      elseif nargin == 3
        %(Tm_model,Time,ReceiverPos)  
             T           = [];%Assigning empty matrix ([]) to Temperature(T)
             gridV_TVGG  = [];%Assigning empty matrix ([]) to gridV_TVGG
             gridV_GTVR  = [];%Assigning empty matrix ([]) to gridV_GTVR
             gridV_GTrop = [];%Assigning empty matrix ([]) to gridV_GTrop
             GPT_grid    = [];%Assigning empty matrix ([]) to GPT_grid
             gridRES     = [];%Assigning empty matrix ([]) to gridRES
             Timevar     = [];%Assigning empty matrix ([]) to Timevar
             
             %*********CHECK TEMPERATURE(TEMP) DEPENDENT MODELS
             if any([strncmpi(Tm_model,'Bevis et al (1992) model',24),strncmpi(Tm_model,'Bevis et al (1995) model',24),strncmpi(Tm_model,'Mendes et al (2000) model',25),...
                     strncmpi(Tm_model,'Schueler et al (2001) model',27),strncmpi(Tm_model,'Yao et al (2014) model',22),strncmpi(Tm_model,'TVGG-Tm model [Jiang et al 2018]',32)])
                
                %COMPUTE SURFACE TEMPERATURE USING GPT MODEL IF TIME & STATION 
                %POSITION INPUTS ARE NOT EMPTY=============================
                if all([~isempty(Time),~isempty(ReceiverPos)])
                  
                   if exist('getMETpara_GPT.m','file')
          
                      try
                         %Call the "getMETpara_GPT.m" external fxn 
                         [~,T] = getMETpara_GPT(UTCtime,[lat,lon, h]);
                         
                      catch  
                           %Call the "getTEMP_GPT.m" internal fxn 
                           T = getTEMP_GPT(UTCtime,lat,lon, h);
                      end 
         
                   else 
                       %Call the "getTEMP_GPT.m" internal fxn 
                       T = getTEMP_GPT(UTCtime,lat,lon, h);
          
                   end  %if exist('getMETpara_GPT.m','file') 
                    
                   %CONVERT TEMPERATURE(T) IN Celcius(°C) to Kelvin(K) 
                   T = T + 273.15; %Temperature in kelvin(K)
                    
                else
                  
                    %ISSUE ERROR MESSAGE 
                    beep%Give a beep sound 
                    errmsg6{1}=sprintf('Surface Temperature is not provided.\n');
                    errmsg6{2}='Please provide Surface Temperature & Try again.';
                    errordlg(errmsg6,'Temperature Input Error','modal')
                    
                    %RETURN EMPTY([]) MATRIX
                    PI =[];
              
                    return
                    
                end %//if all([~isempty(Time),~isempty(ReceiverPos)])
              
             %*********CHECK  GRID FILE & gridRES DEPENDENT MODELS   
             elseif any([strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35)])
                  
                    %ISSUE ERROR MESSAGE  
                    beep%Give a beep sound 
                    errmsg7{1}=sprintf('No input for GPT grid Values & grid resolution.\n');
                    errmsg7{2}='Please provide GPT grid Values & grid resolution & Try again.';
                    errordlg(errmsg7,'GPT grid Values & resolution Input Error','modal') 
                    
                    %RETURN EMPTY([]) MATRIX
                    PI = [];
                         
                    return
          
             end %//if any([strncmpi(Tm_model,'Bevis et al (1992) model',24),strncmpi(Tm_model,'Bevis et al (1995) model',24),strncmpi(Tm_model,'Mendes et al (2000) model',25),...
                 %          strncmpi(Tm_model,'Schueler et al (2001) model',27),strncmpi(Tm_model,'Yao et al (2014) model',22),strncmpi(Tm_model,'TVGG-Tm model [Jiang et al 2018]',32)])
                     
             
      end %//if (any(nargin==[10,9,8,7,6,5,4]))
              
    otherwise   
             beep
             fprintf('\n\nInsuficient Inputs for computation of PI. Check inputs & Try again.\n\n')     
             
             %RETURN EMPTY([]) MATRIX
             PI = [];
             
             return
        
end   %switch nargin

%************FINAL CHECK ON TEMPERATURE(T)
if any([~exist('T','var'),isempty(T)])
    
   %*********CHECK TEMPERATURE(TEMP) DEPENDENT MODELS
   if any([strncmpi(Tm_model,'Bevis et al (1992) model',24),strncmpi(Tm_model,'Bevis et al (1995) model',24),strncmpi(Tm_model,'Mendes et al (2000) model',25),...
           strncmpi(Tm_model,'Schueler et al (2001) model',27),strncmpi(Tm_model,'Yao et al (2014) model',22),strncmpi(Tm_model,'TVGG-Tm model [Jiang et al 2018]',32)])
                
      %COMPUTE SURFACE TEMPERATURE USING GPT MODEL IF TIME & STATION 
      %POSITION INPUTS ARE NOT EMPTY=============================
      if all([~isempty(Time),~isempty(ReceiverPos)])
                  
         if exist('getMETpara_GPT.m','file')
          
            try
               %Call the "getMETpara_GPT.m" external fxn 
               [~,T] = getMETpara_GPT(UTCtime,[lat,lon, h]);
                         
            catch   
                 %Call the "getTEMP_GPT.m" internal fxn 
                 T = getTEMP_GPT(UTCtime,lat,lon, h);
            end  
         
         else  
              %Call the "getTEMP_GPT.m" internal fxn 
              T = getTEMP_GPT(UTCtime,lat,lon, h);
          
         end  %if exist('getMETpara_GPT.m','file') 
                    
         %CONVERT TEMPERATURE(T) IN Celcius(°C) to Kelvin(K) 
         T = T + 273.15; %Temperature in kelvin(K)  

      else
          %ISSUE ERROR MESSAGE 
           beep%Give a beep sound 
           errmsg6{1}=sprintf('Surface Temperature is not provided.\n');
           errmsg6{2}='Please provide Surface Temperature & Try again.';
           errordlg(errmsg6,'Temperature Input Error','modal')
                    
           %RETURN EMPTY([]) MATRIX
           PI =[];
              
           return
                    
      end  %//if all([~isempty(Time),~isempty(ReceiverPos)])

   end
   
end %//if any([~exist('T','var'),isempty(T)])

%===============================END OF INPUT REFORMATTING
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

%                    -------------------------------------
%********************COMPUTE WEIGHTED MEAN TEMPERATURE(Tm) 
%                    --------------------------------------

%*************CHECK FOR Tm MODELS
if any([strncmpi(Tm_model,'Bevis et al (1992) model',24),strncmpi(Tm_model,'Bevis et al (1995) model',24),...
        strncmpi(Tm_model,'Mendes et al (2000) model',25),strncmpi(Tm_model,'Schueler et al (2001) model',27),...
        strncmpi(Tm_model,'Yao et al (2014) model',22),strncmpi(Tm_model,'GTm-I model [Yao et al 2012]',28),...
        strncmpi(Tm_model,'GTm-II model [Yao et al 2013]',29),strncmpi(Tm_model,'GTm-III model [Yao et al 2014]',30),...
        strncmpi(Tm_model,'TVGG-Tm model [Jiang et al 2018]',32),strncmpi(Tm_model,'GTVR-Tm model [Yao et al 2018]',30),...
        strncmpi(Tm_model,'UNB3m model [Leandro et al 2006]',32),strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),...
        strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35),strncmpi(Tm_model,'GTrop model [Sun et al 2019]',28)])
    
  if strncmpi(Tm_model,'Bevis et al (1992) model',24)
      
     %COMPUTE Tm using BEVIS et al(1992,1994) MODEL           
     Tm = 70.2 + 0.72.*T;
     
  elseif strncmpi(Tm_model,'Bevis et al (1995) model',24)

         %COMPUTE Tm using BEVIS et al(1995) MODEL      
         Tm = 86.63 + 0.668.*T;

  elseif strncmpi(Tm_model,'Mendes et al (2000) model',25)
      
         %COMPUTE Tm using Mendes et al (2000) MODEL 
         Tm = 50.4 + 0.789.*T ;
         
  elseif strncmpi(Tm_model,'Schueler et al (2001) model',27)
      
         %COMPUTE Tm using Schueler et al (2001) MODEL 
         Tm = 86.9 + 0.647.*T ;       

  elseif strncmpi(Tm_model,'Yao et al (2014) model',22)
      
         %COMPUTE Tm using Yao et al (2014) MODEL 
         Tm = 43.69 + 0.8116.*T;  
         
  elseif strncmpi(Tm_model,'GTm-I model [Yao et al 2012]',28)
      
         %COMPUTE Tm using GTm-I model [Yao et al 2012]  
         %Call the "Tm_GTm1.m" Function
         Tm = Tm_GTm1(UTCtime,lat,lon,h);
         
  elseif strncmpi(Tm_model,'GTm-II model [Yao et al 2013]',29)
      
         %COMPUTE Tm using GTm-II model [Yao et al 2013] 
         %Call the "Tm_GTm2.m" Function
         Tm = Tm_GTm2(UTCtime,lat,lon,h); 
                  
elseif strncmpi(Tm_model,'GTm-III model [Yao et al 2014]',30)
      
         %COMPUTE Tm using GTm-III model [Yao et al 2012]  
         %Call the "Tm_GTm3.m" Function
         Tm = Tm_GTm3(UTCtime,lat,lon,h);
                 
 elseif strncmpi(Tm_model,'TVGG-Tm model [Jiang et al 2018]',32)
       
        %CONVERT TEMPERATURE(T) IN Kelvin(K) to Celcius(°C) 
        T = T - 273.15;
        
         %COMPUTE Tm using TVGG-Tm model [Jiang et al 2018] 
         %Call the "Tm_TVGG_ERAI.m" Function
         Tm = Tm_TVGG_ERAI(UTCtime,lat,lon,T,gridV_TVGG);       
         
elseif strncmpi(Tm_model,'GTVR-Tm model [Yao et al 2018]',30)
        
         %COMPUTE Tm using GTVR-Tm model [Yao et al 2018] 
         %Call the "Tm_GTVR.m" Function
         Tm = Tm_GTVR(UTCtime,lat,lon,h,gridV_GTVR);
         
  elseif strncmpi(Tm_model,'UNB3m model [Leandro et al 2006]',32)
    
         %COMPUTE Tm using UNB3m model [Leandro et al 2006] 
         %Call the "getMETpara_UNB3m.m" Function       
         [~, ~, ~,Tm] = getMETpara_UNB3m(UTCtime,[lat,lon,h]);
         
  elseif strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28)
      
         %COMPUTE Tm using GPT2w model [Böhm et al 2014]  
         %Call the "getMETpara_GPT2w.m" Function
         [~,~,~,Tm]=getMETpara_GPT2w(UTCtime,[lat,lon,h],GPT_grid,gridRES,Timevar);
         
  elseif strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35)
      
         %COMPUTE Tm using GPT3 model [Landskron & Böhm, 2018] 
         %Call the "getMETpara_GPT3.m" Function 
         [~,~,~,Tm] = getMETpara_GPT3(UTCtime,[lat,lon,h],GPT_grid,gridRES,Timevar);
         
         
  elseif strncmpi(Tm_model,'GTrop model [Sun et al 2019]',28)
      
         %COMPUTE Tm using GTrop model [Sun et al 2019]
         %Call the "TropModel_GTrop.m" Function 
         [~,~,~,~,~,~,Tm] = TropModel_GTrop(UTCtime,[lat,lon,h],gridV_GTrop);
           
  end %//if strncmpi(Tm_model,'Bevis et al (1992) model',24)
  
  %***********NOW COMPUTE PI (from Bevis et al., 1994)
  
  %DEFINE CONSTANT TERMS(from Bevis et al., 1994)
  k1  = 77.6; %refractivity constant in [ K/hPa]
  k2  = 70.4; %refractivity constant in [K/hPa]
  k3  = 3.739E5;%refractivity constant in [K^2/hPa]
  Mw  = 18.0152;%Molar Mass of Water in [g/mol]
  Md  = 28.9644; %Molar Mass of Dry Air  in [g/mol]
  k2p = k2 - k1*Mw/Md;%refractivity constant in [ K/hPa]  {called k2 prime}          
  Rv  = 461.51; %the specific gas constant for water vapor in [J/kg/K]
  
  %k2p = 22.1;% k2-k1*Mw/Md(k1=77.6,k2=70.4,Mw=18.0152,Md=28.9644)

  PI = (10^5)./( ( (k3./Tm) + k2p).*Rv ); %Dimensionless
  
  %****************************END OF USING Tm MODELS TO COMPUTE PI 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  %-------------------------------------------------------------
  %USING Simplified PI  model INSTEAD OF Tm MODELS TO COMPUTE PI
  %--------------------------------------------------------------
  
elseif strncmpi(Tm_model,'Simplified PI  model [Manandhar et al 2017]',43)
    
       %COMPUTE PI Simplified PI  model [Manandhar et al 2017] MODEL 
       %Call the "simplified_PI_model.m" Function 
       [PI] = simplified_PI_model(UTCtime,lat,h);
      
         Tm = zeros(size(PI));%Create zero matrix for Tm 
           
end
%===============================END OF computePI.m
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
%1.      T  : Surface temperature in degree Celcius(°C)                    *
%2.      P  : Surface pressure in millibar(mbar)                           *
%3.    undu : geoid undulation in meters                                   *
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