%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            "MF" is a subroutine that Computes the wet and                *         
%             hydrostatic Tropospheric Mapping Functions from  the various *  
%             Mapping Function  Models                                     *
%USAGE:[MFh, MFw]=MF(UTCtime,lat,lon,hgt,satEL,MFh_type,MFw_type,T,P,e)    *

%INPUT:                                                                    *
%1.   UTCtime  : Receiver reception time in[year,Month,Day,Hour,Min,Sec]   *
%2.      latD  : Station/site geodetic latitude coordinates in [degrees]   *
%3.      lonD  : Station/site geodetic longitude coordinates in [degrees]  *
%4.       hgt  : Station/site Ellipsoidal Hight in [meters]                *
%5.     satEL  : Satellite Elevation in [degrees]                          *
%6.  MFh_type  : Hydrostatic mapping function model type e.g.'Herring'     *
%7.  MFw_type  : Wet mapping function model type e.g.'Chao'                *
%8.         T  : Temperature in [Kelvin]                                   *
%9.         P  : Pressure in [mbar / hpa]                                  *
%10.        e  : Water vapor pressure in [mbar / hpa]                      *
%11.  grid_MFh : Hydrostatic MF grid values in cells extracted from the    *
%                grid file
%12.gridRES_MFh: Hydrostatic MF grid resolution (°) (possible: 1 or 5)     *
%13.  grid_MFw : Wet MF grid values in cells extracted from the grid file  *             
%14.gridRES_MFw: Wet MF grid resolution (°) (possible: 1 or 5)             *              
%15.    Timevar: case 1: no time variation but static quantities           *
%                case 0: with time variation (annual and semiannual terms) *
%--------------------------------------------------------------------------*
%NOTE:
%     INPUTES 10 - 14 CAN BE IGNORED IF MF MODELS ARE NON OF THE VMF/GPT   *
%     MODELS OR REPLACE EACH OF THEM WITH EMPTY([]) MATRIX. I.E.:          *
%     grid_MFh = []; gridRES_MFh = [] ;grid_MFw = []; gridRES_MFw = []     *
%     Timevar = []   AND                                                   *
%USAGE BECOMES:
%[MFh, MFw] = MF(UTCtime,latD,lonD,hgt,satEL,MFh_model,MFw_model,T,P,e) OR *
%[MFh, MFw] = MF(UTCtime,latD,lonD,hgt,satEL,MFh_model,MFw_model,T,P,e,[],[],[],[],[])
%--------------------------------------------------------------------------*
%OUTPUT:                                                                   *
%1.     MFh : Hydrostatic/dry mapping function                             *
%2.     MFw : Non Hydrostatic/wet mapping function                         *
% =========================================================================
% WRITTEN BY:                                                              +
%            OSAH SAMUEL, MSC GEOMATIC ENGINEERING (PhD STUDENT)           +
%            Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com                +
%            Phone:+233(0)246137410 / +233(0)509438484                     +  
%==========================================================================
%***************************************************************************
%***************************************************************************

function [MFh, MFw] = MF(UTCtime,latD,lonD,hgt,satEL,MFh_model,MFw_model,T,P,es,grid_MFh,gridRES_MFh,grid_MFw,gridRES_MFw,Timevar,VMFgrids)

%***CHECK INPUTs
if nargin == 14
    
   VMFgrids = [];%Assigning empty matrix ([]) to VMFgrids
   
   %If flag for Time variation(0 or 1) if not provided
   
   %GET TIME VARIATION STATUS FROM goGPS GUI
   Timevar  = getappdata(0,'Timevar_dry');%time variation  
 
   if isempty(Timevar)
     
      %GET TIME VARIATION STATUS FROM goGPS GUI
      Timevar  = getappdata(0,'Timevar_wet');%time variation 
      
   else
       Timevar = 0;%assign default of 0. meaning time variation (annual and semiannual terms) 
   end
   
elseif nargin == 15
    
       VMFgrids = [];%Assigning empty matrix ([]) to VMFgrids
 
       if isempty(Timevar) 
           
          %GET TIME VARIATION STATUS FROM goGPS GUI
          Timevar  = getappdata(0,'Timevar_dry');%time variation  
 
          if isempty(Timevar)
     
             %GET TIME VARIATION STATUS FROM goGPS GUI
             Timevar  = getappdata(0,'Timevar_wet');%time variation  
          end
          
       else %IF Timevar IS NOT EMPTY([])
           
           %CHECK # OF COLUMNS IN Timevar
           if size(Timevar,2) == 2 %IF #COLUMNS = 2, THEN ... 
              Timevar_h = Timevar(1);  
              Timevar_w = Timevar(2); 
              
           else
               Timevar_h = Timevar; 
               Timevar_w = Timevar; 
              
           end
         
       end %//if isempty(Timevar)  
       
elseif nargin == 16 %IF ALL INPUTS ARE PROVIDED
       
       %CHECK # OF COLUMNS IN Timevar
       if size(Timevar,2) == 2 %IF # COLUMNS = 2, THEN ... 
          Timevar_h = Timevar(1);  
          Timevar_w = Timevar(2); 
          
       else 
            Timevar_h = Timevar; 
            Timevar_w = Timevar; 
              
       end 
          
end %if nargin == 14

    if isempty(Timevar) 
       Timevar = 0;%assign default of 0. meaning time variation (annual and semiannual terms)   
    end 
    
    %ASSIGNMENT
    if any([~exist('Timevar_h','var'),~exist('Timevar_w','var')]) 
       Timevar_h = Timevar;
       Timevar_w = Timevar;
    end
    
  
%***CONVERT lat and lon TO radian
%****DEFINE pi
pi = 3.14159265359d0;
latR = latD.*pi/180; %[radian]
lonR = lonD.*pi/180; %[radian] 
%******************COMPUTE HYDROSTATIC MAPPING FUNCTION
%                  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

%CHECK FOR GLOBAL MAPPING FUNCTION(GMF)
if strncmpi(MFh_model,'GMF',3)
    
    %Call the "GMF.m" sub-routine
   [MFh,MFw_gmf] = GMF(UTCtime,latR,lonR,hgt,satEL);
 
      %CHECK FOR VIENNA MAPPING FUNCTION(VMF)
elseif any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),...
            strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12),...
            strncmpi(MFh_model,'GPT2 (5° x 5°)',14),strncmpi(MFh_model,'GPT2w (1° x 1°)',14),...
            strncmpi(MFh_model,'GPT2w (5° x 5°)',14),strncmpi(MFh_model,'GPT3 (1° x 1°)',14),...
            strncmpi(MFh_model,'GPT3 (5° x 5°)',14),strncmpi(MFh_model,'VMF gridded ZHD',15)])
    
       %Call the "VMF.m" sub-routine
       [MFh,MFw_vmf] = VMF(UTCtime,latD,lonD,hgt,satEL,VMFgrids,grid_MFh,grid_MFw,gridRES_MFh,gridRES_MFw,Timevar_h,Timevar_w,MFh_model,MFw_model);
 
      %CHECK FOR NEILL MAPPING FUNCTION(NMF)      
elseif any([strncmpi(MFh_model,'NMF',3),strncmpi(MFh_model,'Neill',5),...
            strncmpi(MFh_model,'UNB3m',5),strncmpi(MFh_model,'UNB3',4),...
            strncmpi(MFh_model,'UNB',3),any([strncmpi(MFh_model,'GTrop [Sun et al 2019]',22),strfind(MFh_model,'GTrop')])])
       
       %Call the "NMF.m" sub-routine
       [MFh,MFw_nmf] = NMF(UTCtime,latD,hgt,satEL) ; 
       
     %CHECK FOR BLACK & EISNER MAPPING FUNCTION
elseif any([strncmpi(MFh_model,'Black & Eisner',14),strncmpi(MFh_model,'Black and Eisner',16),...
            strncmpi(MFh_model,'Askne & Nordius',15),strncmpi(MFh_model,'Askne and Nordius',17),...
            strncmpi(MFh_model,'EGNOS',5),strncmpi(MFh_model,'MOPS',4)])
       
       %Call the "Black_Eisner_MF.m" sub-routine
       [MFh,MFw_be]=Black_Eisner_MF(satEL);
      
      %CHECK FOR BLACK MAPPING FUNCTION 
elseif strncmpi(MFh_model,'Black',5)
    
       %Call the "Black_MF.m" sub-routine
       [MFh,MFw_b]=Black_MF(satEL,T);
       
       %CHECK FOR IFADIS MAPPING FUNCTION  
elseif  strncmpi(MFh_model,'Ifadis',6)  
    
        %Call the "Ifadis_MF.m" sub-routine
        [MFh,MFw_ifa]=Ifadis_MF(satEL,T,P,es);      
       
       %CHECK FOR HERRING MAPPING FUNCTION 
elseif strncmpi(MFh_model,'Herring',7)
    
       %Call the "Herring_MF.m" sub-routine
       [MFh,MFw_her]=Herring_MF(latR,hgt,satEL,T);
      
       %CHECK FOR CHAO MAPPING FUNCTION(GMF) 
elseif strncmpi(MFh_model,'Chao',4)
    
       %Call the "Chao_MF.m" sub-routine
       [MFh,MFw_chao]=Chao_MF(satEL);
       
       %CHECK FOR MOFFETT MAPPING FUNCTION
elseif strncmpi(MFh_model,'Moffett',7) 
    
       %Call the "Moffett_MF.m" sub-routine
       [MFh,MFw_mof]=Moffett_MF(satEL);
       
       %CHECK FORHOPFIELD MAPPING FUNCTION 
elseif  strncmpi(MFh_model,'Hopfield',8)
    
        %Call the "Hopfield_MF.m" sub-routine
        [MFh,MFw_hop]=Hopfield_MF(satEL);
        
        %CHECK FOR DAVIS MAPPING FUNCTION
elseif any([strncmpi(MFh_model,'Davis',5),strncmpi(MFh_model,'Davis et al',11)])    
    
       %Call the "Davis_MF.m" sub-routine
       MFh=Davis_MF(satEL,T,P,es);
       
       %CHECK FOR SAASTAMOINEN MAPPING FUNCTION 
elseif any([strncmpi(MFh_model,'Saastamoinen',12),strncmpi(MFh_model,'Simple',6)]) 
    
       %Call the "Saastamoinen_MF.m" sub-routine
       [MFh,MFw_saas]=Saastamoinen_MF(satEL); 
       
elseif any([strncmpi(MFh_model,'Modified Hopfield(Goads & Goodman)',34),...
            strncmpi(MFh_model,'Goad & Goodman',14),strncmpi(MFh_model,'GG',2),...
            strncmpi(MFh_model,'Modified Hop',12),strncmpi(MFh_model,'Modified Hopfield',17)])       
       
       %Call the "Goad_Goodman_MF.m" sub-routine 
       [MFh,MFw_gg]=Goad_Goodman_MF(T,P,es,satELEV) ;
       
end %//if strncmpi(MFh_model,'GMF',3)
  
%******************COMPUTE WET MAPPING FUNCTION(MFw)
%                  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%CHECK FOR GLOBAL MAPPING FUNCTION(GMF)
if strncmpi(MFw_model,'GMF',3)
    
    if exist('MFw_gmf','var')
       
       if strncmpi(MFw_model,MFh_model,3) 
            
          MFw = MFw_gmf;
       else 
          %Call the "GMF.m" sub-routine
          [~,MFw] = GMF(UTCtime,latR,lonR,hgt,satEL);
       end 
       
    else 
        %Call the "GMF.m" sub-routine
        [~,MFw] = GMF(UTCtime,latR,lonR,hgt,satEL);
       
    end 
 
      %CHECK FOR VIENNA MAPPING FUNCTION(VMF)
elseif any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),...
            strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12),...
            strncmpi(MFw_model,'GPT2 (5° x 5°)',14),strncmpi(MFw_model,'GPT2w (1° x 1°)',14),...
            strncmpi(MFw_model,'GPT2w (5° x 5°)',14),strncmpi(MFw_model,'GPT3 (1° x 1°)',14),...
            strncmpi(MFw_model,'GPT3 (5° x 5°)',14),strncmpi(MFw_model,'VMF gridded ZWD',15)])
    
       if exist('MFw_vmf','var')
           
          MFw = MFw_vmf;
          
       else    
            %Call the "VMF.m" sub-routine
            [~,MFw] = VMF(UTCtime,latD,lonD,hgt,satEL,VMFgrids,grid_MFh,grid_MFw,gridRES_MFh,gridRES_MFw,Timevar_h,Timevar_w,MFh_model,MFw_model);
       end  
          
      %CHECK FOR NEILL MAPPING FUNCTION(NMF)      
elseif any([strncmpi(MFw_model,'NMF',3),strncmpi(MFw_model,'Neill',5),...
            strncmpi(MFw_model,'UNB3m',5),strncmpi(MFw_model,'UNB3',4),...
            strncmpi(MFw_model,'UNB',3),any([strncmpi(MFw_model,'GTrop [Sun et al 2019]',22),strfind(MFw_model,'GTrop')])])
    
       if exist('MFw_nmf','var')
           
          if strncmpi(MFw_model,MFh_model,5)
       
             MFw = MFw_nmf; 
          else  
              %Call the "NMF.m" sub-routine
              [~,MFw] = NMF(UTCtime,latD,hgt,satEL) ; 
          end 
          
       else  
           %Call the "NMF.m" sub-routine
           [~,MFw] = NMF(UTCtime,latD,hgt,satEL) ; 
       end 
       
     %CHECK FOR BLACK & EISNER MAPPING FUNCTION   
elseif any([strncmpi(MFw_model,'Black & Eisner',14),strncmpi(MFw_model,'Black and Eisner',16),...
            strncmpi(MFw_model,'Askne & Nordius',15),strncmpi(MFw_model,'Askne and Nordius',17),...
            strncmpi(MFw_model,'Callahan',8), strncmpi(MFw_model,'Berman',6), strncmpi(MFw_model,'Mendes',6),...
            strncmpi(MFw_model,'EGNOS',5),strncmpi(MFw_model,'MOPS',4)])
    
      if exist('MFw_be','var')
       
        if strncmpi(MFw_model,MFh_model,17)
            
           MFw = MFw_be;
      else 
           %Call the "Black_Eisner_MF.m" sub-routine
           [~,MFw] = Black_Eisner_MF(satEL);
        end 
        
      else
           %Call the "Black_Eisner_MF.m" sub-routine
           [~,MFw] = Black_Eisner_MF(satEL);
      end
      
      %CHECK FOR BLACK MAPPING FUNCTION 
elseif strncmpi(MFw_model,'Black',5)
    
       if exist('MFw_b','var')
       
          if strncmpi(MFw_model,MFh_model,5)
              
             MFw = MFw_b;
        
          else  
              %Call the "Black_MF.m" sub-routine
              [~,MFw] = Black_MF(satEL,T);
          end
       
       else  
           %Call the "Black_MF.m" sub-routine
           [~,MFw] = Black_MF(satEL,T);
           
       end
       
       %CHECK FOR IFADIS MAPPING FUNCTION  
elseif  strncmpi(MFw_model,'Ifadis',6)  
    
        if exist('MFw_ifa','var')
           
           if strncmpi(MFw_model,MFh_model,6) 
              
              MFw = MFw_ifa;
           else  
               %Call the "Ifadis_MF.m" sub-routine
               [~,MFw] = Ifadis_MF(satEL,T,P,es); 
           end
           
        else  
            %Call the "Ifadis_MF.m" sub-routine
            [~,MFw] = Ifadis_MF(satEL,T,P,es); 
        end  
       
       %CHECK FOR HERRING MAPPING FUNCTION 
elseif strncmpi(MFw_model,'Herring',7)
    
       if exist('MFw_her','var')
       
          if strncmpi(MFw_model,MFh_model,7) 
              
             MFw = MFw_her;
          else  
              %Call the "Herring_MF.m" sub-routine
              [~,MFw] = Herring_MF(latR,hgt,satEL,T);
          end 
          
       else  
           %Call the "Herring_MF.m" sub-routine
           [~,MFw] = Herring_MF(latR,hgt,satEL,T);
       end 
      
       %CHECK FOR CHAO MAPPING FUNCTION(GMF) 
elseif strncmpi(MFw_model,'Chao',4)
    
       if exist('MFw_chao','var')
          
          if strncmpi(MFw_model,MFh_model,4) 
              
             MFw = MFw_chao;
          else  
              %Call the "Chao_MF.m" sub-routine
              [~,MFw] = Chao_MF(satEL);
          end 
          
       else  
           %Call the "Chao_MF.m" sub-routine
           [~,MFw] = Chao_MF(satEL);
       end
       
       %CHECK FOR MOFFETT MAPPING FUNCTION
elseif strncmpi(MFw_model,'Moffett',7) 
    
       if exist('MFw_mof','var')
           
          if strncmpi(MFw_model,MFh_model,7)
       
             MFw = MFw_mof;
          else  
              %Call the "Moffett_MF.m" sub-routine
              [~,MFw]=Moffett_MF(satEL);
          end 
          
       else  
           %Call the "Moffett_MF.m" sub-routine
           [~,MFw]=Moffett_MF(satEL);
       end 
       
       %CHECK FORHOPFIELD MAPPING FUNCTION 
elseif  strncmpi(MFw_model,'Hopfield',8)
        
        if exist('MFw_hop','var')
            
            if strncmpi(MFw_model,MFh_model,8)
       
               MFw = MFw_hop;
            else  
                %Call the "Hopfield_MF.m" sub-routine
                [~,MFw]=Hopfield_MF(satEL);
            end  
            
        else  
            %Call the "Hopfield_MF.m" sub-routine
            [~,MFw]=Hopfield_MF(satEL);
        end   
        
        
       %CHECK FOR SAASTAMOINEN MAPPING FUNCTION 
elseif any([strncmpi(MFw_model,'Saastamoinen',12),strncmpi(MFw_model,'Simple',6)]) 
    
       if exist('MFw_saas','var')
           
          if strncmpi(MFw_model,MFh_model,12)
       
             MFw = MFw_saas;
          else  
              %Call the "Saastamoinen_MF.m" sub-routine
              [~,MFw]=Saastamoinen_MF(satEL);  
          end 
          
       else  
           %Call the "Saastamoinen_MF.m" sub-routine
           [~,MFw]=Saastamoinen_MF(satEL);  
       end 
       
elseif any([strncmpi(MFh_model,'Modified Hopfield(Goads & Goodman)',34),...
            strncmpi(MFh_model,'Goad & Goodman',14),strncmpi(MFh_model,'GG',2),...
            strncmpi(MFh_model,'Modified Hop',12),strncmpi(MFh_model,'Modified Hopfield',17)])       
       
       if exist('MFw_gg','var')
           
          if strncmpi(MFw_model,MFh_model,34) 
       
             MFw = MFw_gg;
          else  
               %Call the "Goad_Goodman_MF.m" sub-routine
               [~,MFw] = Goad_Goodman_MF(T,P,es,satEL) ;    
          end 
          
       else  
            %Call the "Goad_Goodman_MF.m" sub-routine
            [~,MFw] = Goad_Goodman_MF(T,P,es,satEL) ;    
       end
       
end %//if strncmpi(MFw_model,'GMF',3)
%******************************************END OF MF.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%A.SUBROUTINE TO COMPUTE GLOBAL MAPPING FUNCTION(GMF)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
function [GMFh,GMFw] = GMF(UTCtime,lat,lon,hgt,satELEV)

%**************************************************************************
%***DESCRIPTION:
%              This subroutine determines the Global Mapping Functions GMF *
%              given receiver/station position in [lat long hgt],...       * 
%              Satellites Elevation Angle(el),and Modified Julian Day(mjd).*

%             Original source codesis by Böhm et al (2006), modified by ...* 
%             Osah Samuel.                                                 *                                                              
%USAGE:                                                                    *
%      General:                                                            *
%              [GMFh,GMFw] = UTCtime(UTCtime,lat,lon,hgt,satEL)            *       
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%****INPUTs:                                                               *
%1.       UTCtime:.........UTC time in [Year,Month,Day,Hour,Minute,Seconds]*                    
%2.           lat:.........Station geodetic latitude in [radians]          *
%3.           lon:.........Station geodetic in [radians]                   *
%4.           hgt:.........Station ellipsoidal height in [m]               *
%5.         satEL:.........Satellites Elevation Angle in [degrees]         *
% 
%***OUTPUTs:                                                               *
%1.         GMFh:..........Global Hydrostatic mapping function             *
%2.         GMFw:..........Global Wet mapping function                     *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
       
      %1.*****INITIALIZING OUTPUT VARIABLEs
      satel=zeros(size(satELEV,1),size(satELEV,2));%Assign zeros of nxm to satel
      [GMFh,GMFw]=deal(satel);%copy the contents of satel to GMFh,GMFw
      
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
      
                %CONVERT ELEVATION ANGLES IN DEGREES TO RADIANS
                satEL =  abs(satELEV(j,i)) * pi / 180; 
                
                %Call the "gmf.m" function
                [GMFh(j,i),GMFw(j,i)] = gmf(UTCtime(i,:),lat(i),lon(i),hgt(i),satEL);
                
             end
             
          end
          
   else %for different times
       
      [nrow,ncol]=size(satELEV);%Get size of satELEV
       GMFh = cell(nrow,ncol,nrow_time);%Create a cell array of GMFh Output
       GMFw = deal(GMFh);%Create copy of GMFh 
       
      for k =1:nrow_time %Loop over UTCtime for different sets of time
          
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
      
                %CONVERT ELEVATION ANGLES IN DEGREES TO RADIANS
                satEL =  abs(satELEV(j,i)) * pi / 180; 
                
                 %Call the "gmf.m" function
                [GMFh{j,i,k},GMFw{j,i,k}] = gmf(UTCtime(k,:),lat(i),lon(i),hgt(i),satEL);
                                                             
             end
         end
      end
      
   end %//if isequal(Identical,1)
   
else  
     %FIND NUMBER OF ROWS & COLUMNs IN satELEV
     [nrow,ncol]=size(satELEV);%Get size of satELEV  
     
     nrow_time=size(UTCtime,1); %Get the # of Rows in UTCtime 
      
     if nrow_time > 1 %(Indication for different sets of time)
          
       GMFh = cell(nrow,ncol,nrow_time);%Create a cell array of GMFh Output
       GMFw = deal(GMFh);%Create copy of GMFh 
       
       for k =1:nrow_time %Loop over sets if time
          
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
      
              %CONVERT ELEVATION ANGLES IN DEGREES TO RADIANS
              satEL =  abs(satELEV(j,i)) * pi / 180; 
              
              %Call the "gmf.m" function 
              [GMFh{j,i,k},GMFw{j,i,k}] = gmf(UTCtime(k,:),lat(i),lon(i),hgt(i),satEL); 
                
          end                             
          
          end 
          
       end 
       
     else %for a single time info
         %INITIALIZE OUTPUT
         satel=zeros(size(satELEV,1),size(satELEV,2));%Assign zeros of nxm to satel
         
         [GMFh,GMFw]=deal(satel);%copy the contents of satel to GMFh,GMFw

            
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
      
             %CONVERT ELEVATION ANGLES IN DEGREES TO RADIANS 
             satEL =  abs(satELEV(j,i)) * pi / 180;  
             
             %Call the "vmf1.m" function
             [GMFh(j,i),GMFw(j,i)] = gmf(UTCtime,lat(i),lon(i),hgt(i),satEL);
                                                                            
           end
           
        end
        
     end %//if nrow_time > 1
                
end %//if isequal(nrow_time,nrow_pos)

%******************************************END OF GMF.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%A.1.SUBROUTINE TO COMPUTE GLOBAL MAPPING FUNCTION(GMF)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
function [gmfh,gmfw] = gmf(UTCtime,lat,lon,hgt,satEL)

%**************************************************************************
%***DESCRIPTION:
%              This subroutine determines the Global Mapping Functions GMF *
%              given receiver/station position in [lat long hgt],...       * 
%              Satellites Elevation Angle(el),and Modified Julian Day(mjd).*

%             Original source codesis by Böhm et al (2006), modified by ...* 
%             Osah Samuel.                                                 *                                                              
%USAGE:                                                                    *
%      General:                                                            *
%              [gmfh,gmfw] = gmf(UTCtime,lat,lon,hgt,satEL)                *       
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%****INPUTs:                                                               *
%1.       UTCtime:.........UTC time in [Year,Month,Day,Hour,Minute,Seconds]*                    
%2.           lat:.........Station geodetic latitude in [radians]          *
%3.           lon:.........Station geodetic in [radians]                   *
%4.           hgt:.........Station ellipsoidal height in [m]               *
%5.         satEL:.........Satellites Elevation Angle in [radians]         *
% 
%***OUTPUTs:                                                               *
%1.         gmfh:..........Global Hydrostatic mapping function             *
%2.         gmfw:..........Global Wet mapping function                     *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%GET UTC TIME COMPONENTs
%1.******UTCdate
  Yr = UTCtime(:,1);%get Hour
  Mn = UTCtime(:,2);%get Month
 Day = UTCtime(:,3);%get Day
%2.UTCtime
   H = UTCtime(:,4);%get Hour
 MIN = UTCtime(:,5);%get Minute
SECs = UTCtime(:,6);%get Seconds

%******************DEFINE MAPPING COEFFICIENTs
%1.HYDROSTATIC COEFFICIENTs
ah_mean = ...
    [+1.2517d+02, +8.503d-01, +6.936d-02, -6.760d+00, +1.771d-01, ...
      +1.130d-02, +5.963d-01, +1.808d-02, +2.801d-03, -1.414d-03, ...
      -1.212d+00, +9.300d-02, +3.683d-03, +1.095d-03, +4.671d-05, ...
      +3.959d-01, -3.867d-02, +5.413d-03, -5.289d-04, +3.229d-04, ...
      +2.067d-05, +3.000d-01, +2.031d-02, +5.900d-03, +4.573d-04, ...
      -7.619d-05, +2.327d-06, +3.845d-06, +1.182d-01, +1.158d-02, ...
      +5.445d-03, +6.219d-05, +4.204d-06, -2.093d-06, +1.540d-07, ...
      -4.280d-08, -4.751d-01, -3.490d-02, +1.758d-03, +4.019d-04, ...
      -2.799d-06, -1.287d-06, +5.468d-07, +7.580d-08, -6.300d-09, ...
      -1.160d-01, +8.301d-03, +8.771d-04, +9.955d-05, -1.718d-06, ...
      -2.012d-06, +1.170d-08, +1.790d-08, -1.300d-09, +1.000d-10];
  
bh_mean = ...
     [+0.000d+00, +0.000d+00, +3.249d-02, +0.000d+00, +3.324d-02, ...
      +1.850d-02, +0.000d+00, -1.115d-01, +2.519d-02, +4.923d-03, ...
      +0.000d+00, +2.737d-02, +1.595d-02, -7.332d-04, +1.933d-04, ...
      +0.000d+00, -4.796d-02, +6.381d-03, -1.599d-04, -3.685d-04, ...
      +1.815d-05, +0.000d+00, +7.033d-02, +2.426d-03, -1.111d-03, ...
      -1.357d-04, -7.828d-06, +2.547d-06, +0.000d+00, +5.779d-03, ...
      +3.133d-03, -5.312d-04, -2.028d-05, +2.323d-07, -9.100d-08, ...
      -1.650d-08, +0.000d+00, +3.688d-02, -8.638d-04, -8.514d-05, ...
      -2.828d-05, +5.403d-07, +4.390d-07, +1.350d-08, +1.800d-09, ...
      +0.000d+00, -2.736d-02, -2.977d-04, +8.113d-05, +2.329d-07, ...
      +8.451d-07, +4.490d-08, -8.100d-09, -1.500d-09, +2.000d-10];
       
ah_amp = ...
     [-2.738d-01, -2.837d+00, +1.298d-02, -3.588d-01, +2.413d-02, ...
      +3.427d-02, -7.624d-01, +7.272d-02, +2.160d-02, -3.385d-03, ...
      +4.424d-01, +3.722d-02, +2.195d-02, -1.503d-03, +2.426d-04, ...
      +3.013d-01, +5.762d-02, +1.019d-02, -4.476d-04, +6.790d-05, ...
      +3.227d-05, +3.123d-01, -3.535d-02, +4.840d-03, +3.025d-06, ...
      -4.363d-05, +2.854d-07, -1.286d-06, -6.725d-01, -3.730d-02, ...
      +8.964d-04, +1.399d-04, -3.990d-06, +7.431d-06, -2.796d-07, ...
      -1.601d-07, +4.068d-02, -1.352d-02, +7.282d-04, +9.594d-05, ...
      +2.070d-06, -9.620d-08, -2.742d-07, -6.370d-08, -6.300d-09, ...
      +8.625d-02, -5.971d-03, +4.705d-04, +2.335d-05, +4.226d-06, ...
      +2.475d-07, -8.850d-08, -3.600d-08, -2.900d-09, +0.000d+00];
       
bh_amp = ...
     [+0.000d+00, +0.000d+00, -1.136d-01, +0.000d+00, -1.868d-01, ...
      -1.399d-02, +0.000d+00, -1.043d-01, +1.175d-02, -2.240d-03, ...
      +0.000d+00, -3.222d-02, +1.333d-02, -2.647d-03, -2.316d-05, ...
      +0.000d+00, +5.339d-02, +1.107d-02, -3.116d-03, -1.079d-04, ...
      -1.299d-05, +0.000d+00, +4.861d-03, +8.891d-03, -6.448d-04, ...
      -1.279d-05, +6.358d-06, -1.417d-07, +0.000d+00, +3.041d-02, ...
      +1.150d-03, -8.743d-04, -2.781d-05, +6.367d-07, -1.140d-08, ...
      -4.200d-08, +0.000d+00, -2.982d-02, -3.000d-03, +1.394d-05, ...
      -3.290d-05, -1.705d-07, +7.440d-08, +2.720d-08, -6.600d-09, ...
      +0.000d+00, +1.236d-02, -9.981d-04, -3.792d-05, -1.355d-05, ...
      +1.162d-06, -1.789d-07, +1.470d-08, -2.400d-09, -4.000d-10];

%2.WET COEFFICIENTs  
aw_mean = ...
     [+5.640d+01, +1.555d+00, -1.011d+00, -3.975d+00, +3.171d-02, ...
      +1.065d-01, +6.175d-01, +1.376d-01, +4.229d-02, +3.028d-03, ...
      +1.688d+00, -1.692d-01, +5.478d-02, +2.473d-02, +6.059d-04, ...
      +2.278d+00, +6.614d-03, -3.505d-04, -6.697d-03, +8.402d-04, ...
      +7.033d-04, -3.236d+00, +2.184d-01, -4.611d-02, -1.613d-02, ...
      -1.604d-03, +5.420d-05, +7.922d-05, -2.711d-01, -4.406d-01, ...
      -3.376d-02, -2.801d-03, -4.090d-04, -2.056d-05, +6.894d-06, ...
      +2.317d-06, +1.941d+00, -2.562d-01, +1.598d-02, +5.449d-03, ...
      +3.544d-04, +1.148d-05, +7.503d-06, -5.667d-07, -3.660d-08, ...
      +8.683d-01, -5.931d-02, -1.864d-03, -1.277d-04, +2.029d-04, ...
      +1.269d-05, +1.629d-06, +9.660d-08, -1.015d-07, -5.000d-10];
       
bw_mean = ...
     [+0.000d+00, +0.000d+00, +2.592d-01, +0.000d+00, +2.974d-02, ...
      -5.471d-01, +0.000d+00, -5.926d-01, -1.030d-01, -1.567d-02, ...
      +0.000d+00, +1.710d-01, +9.025d-02, +2.689d-02, +2.243d-03, ...
      +0.000d+00, +3.439d-01, +2.402d-02, +5.410d-03, +1.601d-03, ...
      +9.669d-05, +0.000d+00, +9.502d-02, -3.063d-02, -1.055d-03, ...
      -1.067d-04, -1.130d-04, +2.124d-05, +0.000d+00, -3.129d-01, ...
      +8.463d-03, +2.253d-04, +7.413d-05, -9.376d-05, -1.606d-06, ...
      +2.060d-06, +0.000d+00, +2.739d-01, +1.167d-03, -2.246d-05, ...
      -1.287d-04, -2.438d-05, -7.561d-07, +1.158d-06, +4.950d-08, ...
      +0.000d+00, -1.344d-01, +5.342d-03, +3.775d-04, -6.756d-05, ...
      -1.686d-06, -1.184d-06, +2.768d-07, +2.730d-08, +5.700d-09];
       
aw_amp = ...
     [+1.023d-01, -2.695d+00, +3.417d-01, -1.405d-01, +3.175d-01, ...
      +2.116d-01, +3.536d+00, -1.505d-01, -1.660d-02, +2.967d-02, ...
      +3.819d-01, -1.695d-01, -7.444d-02, +7.409d-03, -6.262d-03, ...
      -1.836d+00, -1.759d-02, -6.256d-02, -2.371d-03, +7.947d-04, ...
      +1.501d-04, -8.603d-01, -1.360d-01, -3.629d-02, -3.706d-03, ...
      -2.976d-04, +1.857d-05, +3.021d-05, +2.248d+00, -1.178d-01, ...
      +1.255d-02, +1.134d-03, -2.161d-04, -5.817d-06, +8.836d-07, ...
      -1.769d-07, +7.313d-01, -1.188d-01, +1.145d-02, +1.011d-03, ...
      +1.083d-04, +2.570d-06, -2.140d-06, -5.710d-08, +2.000d-08, ...
      -1.632d+00, -6.948d-03, -3.893d-03, +8.592d-04, +7.577d-05, ...
      +4.539d-06, -3.852d-07, -2.213d-07, -1.370d-08, +5.800d-09];
       
bw_amp = ...
     [+0.000d+00, +0.000d+00, -8.865d-02, +0.000d+00, -4.309d-01, ...
      +6.340d-02, +0.000d+00, +1.162d-01, +6.176d-02, -4.234d-03, ...
      +0.000d+00, +2.530d-01, +4.017d-02, -6.204d-03, +4.977d-03, ...
      +0.000d+00, -1.737d-01, -5.638d-03, +1.488d-04, +4.857d-04, ...
      -1.809d-04, +0.000d+00, -1.514d-01, -1.685d-02, +5.333d-03, ...
      -7.611d-05, +2.394d-05, +8.195d-06, +0.000d+00, +9.326d-02, ...
      -1.275d-02, -3.071d-04, +5.374d-05, -3.391d-05, -7.436d-06, ...
      +6.747d-07, +0.000d+00, -8.637d-02, -3.807d-03, -6.833d-04, ...
      -3.861d-05, -2.268d-05, +1.454d-06, +3.860d-07, -1.068d-07, ...
      +0.000d+00, -2.658d-02, -1.947d-03, +7.131d-04, -3.506d-05, ...
      +1.885d-07, +5.792d-07, +3.990d-08, +2.000d-08, -5.700d-09];

%****DEFINE pi
pi = 3.14159265359d0;

%********COMPUTE MODIFIED JULIAN DATE(MJD) & DAY OF YEAR(doy)

%MODIFIED JULIAN DATE(MJD)
%Call the "utc2JulianDay_DoY.m" function
[~, dmjd,~]=utc2JulianDay_DoY(Yr,Mn,Day,H,MIN,SECs);
    
%****DAY OF YEAR(doy)
%NOTE:Reference day is 28 January. This is taken from Niell (1996) 
%to be consistent     
doy = dmjd  - 44239.d0 + 1 - 28;    
  
%1.COMPUTE ZENITH ANGLE(zd) FROM ELEVATION ANGLEs       
zd =(pi/2)-satEL; % zd in [radian]
     
%****DEFINE degree n and order m
nmax = 9;
mmax = 9;

%COMPUTE unit vectors(x,y,z)
x = cos(lat)*cos(lon);
y = cos(lat)*sin(lon);
z = sin(lat);
  
%COMPUTE  Legendre polynomials
V(1,1) = 1.0D0;
W(1,1) = 0.0D0;
V(2,1) = z * V(1,1);
W(2,1) = 0.0;

for n = 2:nmax
    V(n+1,1) = ((2*n-1) * z * V(n,1) - (n-1) * V(n-1,1)) / n;
    W(n+1,1) = 0.0D0;
end  
for m = 1:nmax
    V(m+1,m+1) = (2*m-1) * (x*V(m,m) - y*W(m,m));
    W(m+1,m+1) = (2*m-1) * (x*W(m,m) + y*V(m,m));
    if (m < nmax) 
       V(m+2,m+1) = (2*m+1) * z * V(m+1,m+1);
       W(m+2,m+1) = (2*m+1) * z * W(m+1,m+1);
    end  
    for n = m+2:nmax
        V(n+1,m+1) = ((2*n-1)*z*V(n,m+1) - (n+m-1)*V(n-1,m+1)) / (n-m);
        W(n+1,m+1) = ((2*n-1)*z*W(n,m+1) - (n+m-1)*W(n-1,m+1)) / (n-m);
    end  
end  
%*****HYDROSTATIC MAPPING FUNCTION(MFh)
bh = 0.0029;
c0h = 0.062;

if (lat < 0) %Southern Hemisphere
  phh  = pi;
  c11h = 0.007;
  c10h = 0.002;
else              %Northern Hemisphere
     phh  = 0;
     c11h = 0.005;
     c10h = 0.001;
end  

ch = c0h + ((cos(doy/365.25d0*2*pi + phh)+1)*c11h/2 + c10h)*(1-cos(lat));

ahm = 0.d0;
aha = 0.d0;
  i = 0;
for n = 0:nmax
    for m = 0:n
        i = i+1;
        ahm = ahm + (ah_mean(i)*V(n+1,m+1) + bh_mean(i)*W(n+1,m+1));
        aha = aha + (ah_amp(i) *V(n+1,m+1) + bh_amp(i) *W(n+1,m+1));
    end  
end  
ah  = (ahm + aha*cos(doy/365.25d0*2.d0*pi))*1d-5;

sine   = sin(pi/2 - zd);
cose   = cos(pi/2 - zd);
beta   = bh./( sine + ch  );
gamma  = ah./( sine + beta);
topcon = (1.d0 + ah./(1.d0 + bh./(1.d0 + ch)));
gmfh   = topcon./(sine+gamma);

%HEIGHT CORRECTION FOR HYDROSTATIC MF from Niell (1996)
a_ht   = 2.53d-5;
b_ht   = 5.49d-3;
c_ht   = 1.14d-3;
hs_km  = hgt/1000.d0;

beta   = b_ht./( sine + c_ht );
gamma  = a_ht./( sine + beta);
topcon = (1.d0 + a_ht./(1.d0 + b_ht./(1.d0 + c_ht)));
ht_corr_coef = 1./sine - topcon./(sine + gamma);
ht_corr      = ht_corr_coef .* hs_km;
gmfh         = gmfh + ht_corr;

%WET MAPPING FUNCTION(MFw)
bw  = 0.00146;
cw  = 0.04391;
awm = 0.d0;
awa = 0.d0;
i = 0;
for n = 0:nmax
    for m = 0:n
        i = i+1;
        awm = awm + (aw_mean(i)*V(n+1,m+1) + bw_mean(i)*W(n+1,m+1));
        awa = awa + (aw_amp(i) *V(n+1,m+1) + bw_amp(i) *W(n+1,m+1));
    end  
end  
aw =  (awm + awa*cos(doy/365.25d0*2*pi))*1d-5;

beta   = bw./( sine + cw );
gamma  = aw./( sine + beta);
topcon = (1.d0 + aw./(1.d0 + bw./(1.d0 + cw)));
gmfw   = topcon./(sine+gamma);

%******************************************END OF gmf.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%                 *END OF GLOBAL MAPPING FUNCTION*
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%B.SUBROUTINE TO COMPUTE NEILL MAPPING FUNCTION(NMF)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
function [NMFh,NMFw] = NMF(UTCtime,lat,hgt,satELEV)

%**************************************************************************
%***DESCRIPTION:
%              This subroutine determines the Global Mapping Functions GMF *
%              given receiver/station position in [lat long hgt],...       * 
%              Satellites Elevation Angle(el),and Modified Julian Day(mjd).*

%             Original source codesis by Böhm et al (2006), modified by ...* 
%             Osah Samuel.                                                 *                                                              
%USAGE:                                                                    *
%      General:                                                            *
%              [NMFh,NMFw] = UTCtime(UTCtime,lat,hgt,satEL)                *       
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%****INPUTs:                                                               *
%1.       UTCtime:.........UTC time in [Year,Month,Day,Hour,Minute,Seconds]*                    
%2.           lat:.........Station geodetic latitude in [Degrees]          *
%3.           hgt:.........Station ellipsoidal height in [m]               *
%4.         satEL:.........Satellites Elevation Angle in [degrees]         *
% 
%***OUTPUTs:                                                               *
%1.         NMFh:..........Neill Hydrostatic mapping function              *
%2.         NMFw:..........Neill Wet mapping function                      *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
       
      %1.*****INITIALIZING OUTPUT VARIABLEs
      satel=zeros(size(satELEV,1),size(satELEV,2));%Assign zeros of nxm to satel
      [NMFh,NMFw]=deal(satel);%copy the contents of satel to NMFh,NMFw
      
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
      
                %CONVERT ELEVATION ANGLES IN DEGREES TO RADIANS
                satEL =  abs(satELEV(j,i)) * pi / 180; 
                
                %Call the "nmf.m" function
                [NMFh(j,i),NMFw(j,i)] = nmf(UTCtime(i,:),lat(i),hgt(i),satEL);
                
             end
             
          end
          
   else %for different times
       
      [nrow,ncol]=size(satELEV);%Get size of satELEV
       NMFh = cell(nrow,ncol,nrow_time);%Create a cell array of NMFh Output
       NMFw = deal(NMFh);%Create copy of NMFh 
       
      for k =1:nrow_time %Loop over UTCtime for different sets of time
          
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
      
                %CONVERT ELEVATION ANGLES IN DEGREES TO RADIANS
                satEL =  abs(satELEV(j,i)) * pi / 180; 
                
                 %Call the "nmf.m" function
                [NMFh{j,i,k},NMFw{j,i,k}] = nmf(UTCtime(k,:),lat(i),hgt(i),satEL);
                                                             
             end
         end
      end
      
   end %//if isequal(Identical,1)
   
else  
     %FIND NUMBER OF ROWS & COLUMNs IN satELEV
     [nrow,ncol]=size(satELEV);%Get size of satELEV  
     
     nrow_time=size(UTCtime,1); %Get the # of Rows in UTCtime 
      
     if nrow_time > 1 %(Indication for different sets of time)
          
       NMFh = cell(nrow,ncol,nrow_time);%Create a cell array of NMFh Output
       NMFw = deal(NMFh);%Create copy of NMFh 
       
       for k =1:nrow_time %Loop over sets if time
          
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
      
              %CONVERT ELEVATION ANGLES IN DEGREES TO RADIANS
              satEL =  abs(satELEV(j,i)) * pi / 180; 
              
              %Call the "nmf.m" function 
              [NMFh{j,i,k},NMFw{j,i,k}] = nmf(UTCtime(k,:),lat(i),hgt(i),satEL); 
                
          end                             
          
          end 
          
       end 
       
     else %for a single time info
         %INITIALIZE OUTPUT
         satel=zeros(size(satELEV,1),size(satELEV,2));%Assign zeros of nxm to satel
         
         [NMFh,NMFw]=deal(satel);%copy the contents of satel to NMFh,NMFw

            
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
      
             %CONVERT ELEVATION ANGLES IN DEGREES TO RADIANS 
             satEL =  abs(satELEV(j,i)) * pi / 180;  
             
             %Call the "nmf.m" function
             [NMFh(j,i),NMFw(j,i)] = nmf(UTCtime,lat(i),hgt(i),satEL);
                                                                            
           end
           
        end
        
     end %//if nrow_time > 1
                
end %//if isequal(nrow_time,nrow_pos)

%******************************************END OF NMF.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%B.1.SUBROUTINE TO COMPUTE NEILL MAPPING FUNCTION(NMF)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
function [nmfh,nmfw] = nmf(UTCtime,lat,hgt,satEL) 

%***************************************************************************
%DESCRIPTION:                                                              * 
%            "nmf" is a subroutine that Computes the wet and dry ...       *
%             Niell mapping function[NMF] 1996                             *
%             The coefficients of the hydrostatic mapping function depend  * 
%             on the latitude and height above sea level of the receiver   *
%             station, and on the day of the year(DOY). The wet mapping    * 
%             function on the other hand, depends only on latitude.        *
%USAGE:                                                                    *
%      General:                                                            *
%              [nmfh,nmfw] = nmf(UTCtime,lat,hgt,satEL)                    *       
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%****INPUTs:                                                               *
%1.       UTCtime:.........UTC time in [Year,Month,Day,Hour,Minute,Seconds]*                         
%2.           lat:.........Station geodetic latitude in [degrees]          *
%4.           hgt:.........Station ellipsoidal height in [meters]          *
%5.         satEL:.........Satellites Elevation Angle in [radians]         *
% 
%***OUTPUTs:                                                               *
%1.         nmfh:..........Neill Hydrostatic mapping function              *
%2.         nmfw:..........Neill Wet mapping function                      *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-           
%==========================================================================+
%REFERENCE:                                                                *
%          Niell, A.E.(1996).Global mapping functions for the atmosphere   * 
%          delay at radio wavelengths. J. Geophys. Res.,101, pp.3227-3246  *
%==========================================================================+
%Written by: Samuel Osah,Msc Geomatic Engineering ,2016                    *    
%             Email: osahsamuel@yahoo.ca                                   *
%             Tel:+233 (0)246137410/+233 (0)509438484                      *
%==========================================================================
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

%Coefficient for changing DOY to radian 
doy2rad = 2*pi/365.25d0; 
 
%LATITUDE ARRAY FOR THE HYDROSTATIC & WET MAPPING FUNCTION COEFFICIENTs
TableLat = [15,30,45,60,75]';

%*********************Initialize NEILL look-up table
%1.*******HYDROSTATIC MAPPING COEFFICIENTs
%AVERAGE OF COEFFICIENTs a,b and c CORRESPONDING 
%             TO THE GIVEN LATITUDE
abc_avg = [1.2769934  2.9153695  62.610505  %15
           1.2683230  2.9152299  62.837393  %30
           1.2465397  2.9288445  63.721774  %45
           1.2196049  2.9022565  63.824265  %60
           1.2045996  2.9024912  64.258455]*1.0e-3 ;%75
           %    a          b         c  

%SEASONAL VARIATION/AMPLITUDE OF COEFFICIENTs a,b and c CORRESPONDING 
%                  TO THE GIVEN LATITUDE
abc_amp =  [   0          0          0              %15
           1.2709626  2.1414979  9.0128400          %30 
           2.6523662  3.0160779  4.3497037          %45
           3.4000452  7.2562722  84.795348          %60
           4.1202191  11.723375  170.37206]*1.0e-5; %75
           %    a          b         c           
 
%HEIGHT CORRECTION COEFFICIENTs   
a_ht = 2.53d-5;
b_ht = 5.49d-3;
c_ht = 1.14d-3;
                                        
%2.***************WET MAPPING COEFFICIENTs
%******AVERAGE OF COEFFICIENTs a,b and c CORRESPONDING 
%                    TO THE GIVEN LATITUDE
abc_w2po = [5.8021897*1.0e-4  1.4275268*1.0e-3  4.3472961*1.0e-2    %15
            5.6794847*1.0e-4  1.5138625*1.0e-3  4.6729510*1.0e-2    %30 
            5.8118019*1.0e-4  1.4572752*1.0e-3  4.3908931*1.0e-2    %45
            5.9727542*1.0e-4  1.5007428*1.0e-3  4.4626982*1.0e-2    %60 
            6.1641693*1.0e-4  1.7599082*1.0e-3  5.4736038*1.0e-2];  %75 
            %    a                   b               c              
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%********INTERPOLATE AVERAGE & AMPLITUDE COEFFICIENT VALUEs 
%Call the "int_Neill_Mapping_COE.m" function
[avg,amp,abc_wet] = int_Neill_Mapping_COE(lat,TableLat,abc_avg,abc_amp,abc_w2po);

%***************COMPUTING MAPPING FUNCTIONs (NMFh,NMFw)
%GET Wet Mapping Coefficients
aw=abc_wet(:,1);%Assigning all rows,1st column to aw
bw=abc_wet(:,2);%Assigning all rows,2nd column to bw
cw=abc_wet(:,3);%Assigning all rows,3rd column to cw

%1.*****INITIALIZING OUTPUT VARIABLEs
%Mapping Functions
satel=zeros(size(satEL,1),size(satEL,2));%Assign zeros of nxm to satel
[nmfh,nmfw]=deal(satel);%copy the contents of satel to MFh,MFw

%FIND NUMBER OF ROWS & COLUMNs IN satEL
[nrow,ncol]=size(satEL);

%**************COMPUTE DAY OF YEAR(DOY)   
%Call the "utc2JulianDay_DoY.m" function
[~, ~,DOY]=utc2JulianDay_DoY(Yr,Mn,Day,H,MIN,SECs);
           
for i=1:ncol %Loop over Stations
    
   for j=1:nrow %Loop over Satellites 
      
   %*************COMPUTE HYDROSTATIC MAPPING COEFFICIENTs (a,b,c)     
    %************DEAL WITH SOUTHERN HEMISPHERE & YEARLY VARIATION
    %NOTE:  
         %Since no data from the southern hemisphere were used in developing
         %these mapping functions, the inversion of the seasons has been
         %accounted for simply by adding half a year to the phase for
         %southern latitudes. 
    if (lat(i) < 0) %if Southern Hemisphere
      DOY = DOY + 365.25/2; %Add half a year to the phase
    end 
    %***COMPUTE COS OF PHASE ANGLE
    COSphase = cos((DOY - 28) * doy2rad );

    %********NOW, COMPUTE HYDROSTATIC COFFICIENT VALUEs(ah,bh,ch)
    ah = avg(i,1) - amp(i,1) * COSphase;
    bh = avg(i,2) - amp(i,2) * COSphase;
    ch = avg(i,3) - amp(i,3) * COSphase;

    %*****COMPUTE sine ELEVATION ANGLE
    sine = sin(satEL(j,i));

    %********COMPUTE NEILL HYDROSTATIC MAPPING FUNCTION value
    beta   = bh/( sine + ch  );
    gamma  = ah/( sine + beta);
    topcon = (1+ ah/(1 + bh/(1 + ch)));
    nmfh   = topcon/(sine+gamma); 
    
    %*********APPLY HEIGHT CORRECTION     
    hs_km  = hgt(i)/1000.d0;%convert height to km
    beta         = b_ht/( sine + c_ht);
    gamma        = a_ht/( sine + beta);
    topcon       = (1 + a_ht/(1 + b_ht/(1 + c_ht)));
    ht_corr_coef = 1/sine - topcon/(sine + gamma);
    ht_corr      = ht_corr_coef * hs_km;
    nmfh(j,i)        = nmfh + ht_corr; 
    
    %********COMPUTE NEILL WET MAPPING FUNCTION value
    beta   = bw(i)/( sine + cw(i) );
    gamma  = aw(i)/( sine + beta);
    topcon = (1 + aw(i)/(1.d0 + bw(i)/(1 + cw(i))));
    nmfw(j,i)   = topcon/(sine+gamma);
      
   end  %for j=1:nrow
   
end %for i=1:nrow
%******************************************END OF nmf.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%B.2.SUBROUTINE TO INTERPOLATE NEILL MAPPING COEFFICIENTs
function [avg,amp,abc_wet] = int_Neill_Mapping_COE(UserLat,TableLat,HydavgValues...
                                                   ,ampValues,WetavgValues)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            The function "int_Neill_Mapping_CO" Interpolates and Computes * 
%            Hydrostatic & Wet mapping function coefficients(a,b,c) from...*  
%            averaged  and seasonal variation/Amplitude values used ...    * 
%            for tropospheric delay prediction given receiver latitude,... * 
%            based on Neill Mapping Function(NMF)Model.Parameters above... *    
%            |lat|<=15° and |lat|>=75° are extracted directly while ...    * 
%            Parameters for latitudes 15°<|lat|<75° are linearly...        *
%            interpolated between values of two  closest latitudes.        *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%INPUT:                                                                    *
%       UserLat: Station/User geodetic latitude(s) in degrees(deg) (n x 1) *
%      TableLat: Latitude array for the Hydrostatic & Wet mapping function *
%                coefficients in degrees(deg).                             *
%  HydavgValues: Hydrostatic Average Values                                *
%     ampValues: Amplitude Values for Hydrostatic Mapping function         *
%  WetavgValues: Wet Average Values                                        *

%AVERAGE OF COEFFICIENTs a,b and c CORRESPONDING TO THE GIVEN LATITUDE     *
%HydavgValues =[1.2769934  2.9153695  62.610505         %15                *
%               1.2683230  2.9152299  62.837393         %30                *
%               1.2465397  2.9288445  63.721774         %45                *
%               1.2196049  2.9022565  63.824265         %60                *
%               1.2045996  2.9024912  64.258455]*1.0e-3;%75                *
%                   a          b         c                                 *
%AMPLITUDE OF COEFFICIENTs a,b and c CORRESPONDING TO THE GIVEN LATITUDE   *
%ampValues =  [   0          0          0               %15                *
%              1.2709626  2.1414979  9.0128400          %30                *
%              2.6523662  3.0160779  4.3497037          %45                *
%              3.4000452  7.2562722  84.795348          %60                *
%              4.1202191  11.723375  170.37206]*1.0e-5; %75                *
%                  a          b         c                                  *         
%2.********WET MAPPING FUNCTION                                            *
%WetavgValues =[5.8021897*1.0e-4  1.4275268*1.0e-3  4.3472961*1.0e-2  %15  *
%               5.6794847*1.0e-4  1.5138625*1.0e-3  4.6729510*1.0e-2  %30  *
%               5.8118019*1.0e-4  1.4572752*1.0e-3  4.3908931*1.0e-2  %45  *
%               5.9727542*1.0e-4  1.5007428*1.0e-3  4.4626982*1.0e-2  %60  *
%               6.1641693*1.0e-4  1.7599082*1.0e-3  5.4736038*1.0e-2];%75  *
%                    a                  b                c               
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%OUTPUT:                                                                   *
%       avg - Average(a,b,c) hydrostatic MF coefficients (nx3)             *
%       amp - seasonal  variation/Amplitude coefficients of a, b and c(nx3)*
%       abc_wet-coefficients of a, b and c(nx3)                            *
%WHERE:
%      nx3-->Means n Rows,three(3) Columns (e.g. 4x3)
%      n-->Indicates the Number of User/Receiver Input Coordinates(Latitude)
%          That is, each row in (avg,amp,abc_wet) represents Receiver Coord*
%     The Columns arrangements are as follows:
%          1       2       3       
%     avg:[a       b       c  ]
%     amp:[a       b       c  ]
% abc_wet:[a       b       c  ]
% =========================================================================
%REFERENCE:                                                                *
%          A.E.Niell, Global mapping functions for the atmosphere delay at * 
%          radio wavelengths, 1996                                         *
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************
%***************************************************************************
% Latitude=[15;16;25;30;35;45;60;75;80];
%*********ASSIGNMENT
abc_avg=HydavgValues;
abc_amp=ampValues;
abc_w2po=WetavgValues;

%****************INTERPOLATIONS  
%******INITIALIZE OUTPUT
 
 %***HYDROSTATIC
 avg=zeros(length(UserLat),size(abc_avg,2));%Assigning zeros to avg 
 amp=deal(avg);%copy the contents avg to amp
 %***WET
 abc_wet=deal(avg);%copy the contents avg to abc_wet
 
 for i=1:length(UserLat) %Loop through user given Latitude coordinates
     
     for j=1:length(TableLat) %Loop through Look-up Table
         
        if (j==1 & UserLat(i)<= TableLat(j))
           
           %***HYDROSTATIC
           avg(i,:) = abc_avg(j,:); %Average(a,b,c) hydrostatic MF coefficients
           amp(i,:) = abc_amp(j,:); %Amplitude of coefficients a, b and c
           
           %***WET
           abc_wet(i,:) = abc_w2po(j,:); %Wet coefficients of a, b and c 
           
        elseif (j==length(TableLat) & UserLat(i)>= TableLat(j))
               
               %***HYDROSTATIC
               avg(i,:) = abc_avg(j,:); %Average(a,b,c) hydrostatic MF coefficients
               amp(i,:) = abc_amp(j,:); %Amplitude of coefficients a, b and c
               
               %***WET
               abc_wet(i,:) = abc_w2po(j,:); %Wet coefficients of a, b and c 
        else     
            if (UserLat(i)== TableLat(j))
              
              %***HYDROSTATIC  
              avg(i,:) = abc_avg(j,:); %Average a hydrostatic MF coefficients
              amp(i,:) = abc_amp(j,:); %Amplitude of coefficients a, b and c
              
              %***WET
              abc_wet(i,:) = abc_w2po(j,:); %Wet coefficients of a, b and c 
              
            elseif (UserLat(i)>TableLat(j,1) && TableLat(j,1)<TableLat(j+1,1))
                  %********HYDROSTATIC
                  diff_avg= abc_avg(j+1,:)-abc_avg(j,:);%Difference in average Met parameters between two closest latitudes
                  diff_amp= abc_amp(j+1,:)-abc_amp(j,:);%Difference in seasonal variation of Met parameters between two closest latitudes
                  diff_latRMF=UserLat(i)-TableLat(j,1);%Difference in latitude values between user Lat and closest latitude of Met parameters
                  diff_latMF=TableLat(j+1,1)-TableLat(j,1);%Difference in latitude values between user Lat and closest latitude of Met parameters
                 
                  %********WET
                  diff_wet= abc_w2po(j+1,:)-abc_w2po(j,:);
                  
                  %***********PERFORM INTERPOLATION
                  %***HYDROSTATIC
                  avg(i,:) = ((diff_latRMF/diff_latMF).*diff_avg)+abc_avg(j,:); %Average Meteorological parameters
                  amp(i,:) = ((diff_latRMF/diff_latMF).*diff_amp)+ abc_amp(j,:);%Seasonal variation Meteorological parameters          
                 
                  %***WET
                  abc_wet(i,:) = ((diff_latRMF/diff_latMF).*diff_wet)+ abc_w2po(j,:);
                  
            end %if (latR(i)== lat_MF(j))
        end %if (j==1 & latR(i)<= lat_MF(j))
        
       if (j+1)>length(TableLat)
             
          break
       end  %if j+1>length(lat_MF) 
        
     end %for j=1:length(lat_MF)
 end %for i=1:length(latR)
%******************************END OF int_Neill_Mapping_COE.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%                 *END OF NEILL MAPPING FUNCTION*
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

%C.SUBROUTINE TO COMPUTE VIENNA MAPPING FUNCTION(VMF)
function [VMFh,VMFw] = VMF(UTCtime,latD,lonD,hell,satELEV,VMFgrids,grid_h,grid_w,grid_resh,grid_resw,Timevar_h,Timevar_w,MFh_model,MFw_model)

%***********COMPUTE VMF COEFFICIENTS

%IF HYDROSTATIC & WET MAPPING FUNCTION(MF) MODELS ARE THE SAME
if any([strncmpi(MFh_model,MFw_model,14),strfind(MFh_model,MFw_model),any([strncmpi(MFh_model,'VMF gridded ZHD',15),...
        strncmpi(MFw_model,'VMF gridded ZWD',15)])])
    
   %CHECK FOR VIENNA MAPPING FUNCTION TYPE (VMF1 / VMF3)
   %********IF VMF1
   if all([any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),strncmpi(MFh_model,'VMF g',5),...
                strfind(MFh_model,'VMF g'),strncmpi(MFh_model,'GPT2 (5° x 5°)',14),strncmpi(MFh_model,'GPT2w (1° x 1°)',14),...
                strncmpi(MFh_model,'GPT2w (5° x 5°)',14)]),...
           any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),strncmpi(MFw_model,'VMF Z',5),...
                strfind(MFw_model,'VMF Z'),strncmpi(MFw_model,'GPT2 (5° x 5°)',14),strncmpi(MFw_model,'GPT2w (1° x 1°)',14),...
                strncmpi(MFw_model,'GPT2w (5° x 5°)',14)])])
       
       %Call the "VMF1coe.m" sub-routine
       [ah1,aw1] = VMF1coe(UTCtime,latD,lonD,hell,VMFgrids,grid_h,grid_w,grid_resh,grid_resw,Timevar_h,Timevar_w,MFh_model,MFw_model);
       
       ah1_flag = 1;
       aw1_flag = 1;
       ah1_aw1  = 1; %flag to indicate 
       
elseif all([any([strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12),strncmpi(MFh_model,'VMF g',5),...
                 strfind(MFh_model,'VMF g'),strncmpi(MFh_model,'GPT3 (1° x 1°)',14),strncmpi(MFh_model,'GPT3 (5° x 5°)',14)]),...
            any([strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12),strncmpi(MFw_model,'VMF g',5),...
                 strfind(MFw_model,'VMF g'),strncmpi(MFw_model,'GPT3 (1° x 1°)',14),strncmpi(MFw_model,'GPT2w (5° x 5°)',14),])])
     
       %Call the "VMF3coe.m" sub-routine  
       [ah3,aw3] = VMF3coe(UTCtime,latD,lonD,hell,VMFgrids,grid_h,grid_w,grid_resh,grid_resw,Timevar_h,Timevar_w,MFh_model,MFw_model); 
       
       ah3_flag = 1;
       aw3_flag = 1;
       ah3_aw3  = 1;
   end

else
    if all([any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),strncmpi(MFh_model,'VMF g',5),...
                 strfind(MFh_model,'VMF g'),strncmpi(MFh_model,'GPT2 (5° x 5°)',14),strncmpi(MFh_model,'GPT2w (1° x 1°)',14),...
                 strncmpi(MFh_model,'GPT2w (5° x 5°)',14)]),...
           ~any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),strncmpi(MFw_model,'VMF g',5),...
                 strfind(MFw_model,'VMF g'),strncmpi(MFw_model,'GPT2 (5° x 5°)',14),strncmpi(MFw_model,'GPT2w (1° x 1°)',14),...
                 strncmpi(MFw_model,'GPT2w (5° x 5°)',14)])])
            
            
       %Call the "VMF1coe.m" sub-routine
       [ah1,aw1] = VMF1coe(UTCtime,latD,lonD,hell,VMFgrids,grid_h,grid_w,grid_resh,grid_resw,Timevar_h,Timevar_w,MFh_model,MFw_model);
       
       ah1_flag = 1;
       aw1_flag = 0;
       ah1_aw1  = 0;
       
    elseif all([~any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),strncmpi(MFh_model,'VMF g',5),...
                      strfind(MFh_model,'VMF g'),strncmpi(MFh_model,'GPT2 (5° x 5°)',14),strncmpi(MFh_model,'GPT2w (1° x 1°)',14),...
                      strncmpi(MFh_model,'GPT2w (5° x 5°)',14)]),...
                 any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),strncmpi(MFw_model,'VMF g',5),...
                      strfind(MFw_model,'VMF g'),strncmpi(MFw_model,'GPT2 (5° x 5°)',14),strncmpi(MFw_model,'GPT2w (1° x 1°)',14),...
                      strncmpi(MFw_model,'GPT2w (5° x 5°)',14)])])
                        
           %Call the "VMF1coe.m" sub-routine
           [ah1,aw1] = VMF1coe(UTCtime,latD,lonD,hell,VMFgrids,grid_h,grid_w,grid_resh,grid_resw,Timevar_h,Timevar_w,MFh_model,MFw_model);
       
           aw1_flag = 1;
           ah1_flag = 0;
           ah1_aw1  = 0;
           
           
    elseif all([any([strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12),strncmpi(MFh_model,'VMF g',5),...
                     strfind(MFh_model,'VMF g'),strncmpi(MFh_model,'GPT3 (1° x 1°)',14),strncmpi(MFh_model,'GPT3 (5° x 5°)',14)]),...
               ~any([strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12),strncmpi(MFw_model,'VMF g',5),...
                     strfind(MFw_model,'VMF g'),strncmpi(MFw_model,'GPT3 (1° x 1°)',14),strncmpi(MFw_model,'GPT2w (5° x 5°)',14),])])
     
           %Call the "VMF3coe.m" sub-routine  
           [ah3,aw3] = VMF3coe(UTCtime,latD,lonD,hell,VMFgrids,grid_h,grid_w,grid_resh,grid_resw,Timevar_h,Timevar_w,MFh_model,MFw_model); 
       
           ah3_flag = 1;
           aw3_flag = 0;
           ah3_aw3  = 0;   
       
    elseif all([~any([strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12),strncmpi(MFh_model,'VMF g',5),...
                      strfind(MFh_model,'VMF g'),strncmpi(MFh_model,'GPT3 (1° x 1°)',14),strncmpi(MFh_model,'GPT3 (5° x 5°)',14)]),...
                 any([strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12),strncmpi(MFw_model,'VMF g',5),...
                      strfind(MFw_model,'VMF g'),strncmpi(MFw_model,'GPT3 (1° x 1°)',14),strncmpi(MFw_model,'GPT2w (5° x 5°)',14),])])
     
           %Call the "VMF3coe.m" sub-routine  
           [ah3,aw3] = VMF3coe(UTCtime,latD,lonD,hell,VMFgrids,grid_h,grid_w,grid_resh,grid_resw,Timevar_h,Timevar_w,MFh_model,MFw_model); 
       
           ah3_flag = 0;
           aw3_flag = 1;
           ah3_aw3  = 0;
           
           
    end %//if all([any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),strncmpi(MFh_model,'VMF Z',5),...
        %               strfind(MFh_model,'VMF Z'),strncmpi(MFh_model,'GPT2 (5° x 5°)',14),strncmpi(MFh_model,'GPT2w (1° x 1°)',14),...
        %               strncmpi(MFh_model,'GPT2w (5° x 5°)',14)]),...
        %         ~any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),strncmpi(MFw_model,'VMF Z',5),...
        %               strfind(MFw_model,'VMF Z'),strncmpi(MFw_model,'GPT2 (5° x 5°)',14),strncmpi(MFw_model,'GPT2w (1° x 1°)',14),...
        %               strncmpi(MFw_model,'GPT2w (5° x 5°)',14)])])
    
    
end %//any([strncmpi(MFh_model,MFw_model,14),strfind(MFh_model,MFw_model),any([strncmpi(MFh_model,'VMF ZHD',7),...
    %       strncmpi(MFw_model,'VMF ZWD',7)])])
    
if any([exist('ah1_flag','var'),exist('aw1_flag','var'),exist('ah1_aw1','var')])
                    
   if any([all([ah1_flag == 1,aw1_flag == 1]),ah1_aw1  == 1])
       
       ah = ah1;
       
   elseif any([ah1_flag == 1,all([aw1_flag == 0,ah1_aw1 == 0])])
       
          ah = ah1;
          
   elseif any([aw1_flag == 1,all([ah1_flag == 0,ah1_aw1 == 0])])  
       
          ah = aw1;
          
   end
    
elseif any([exist('ah3_flag','var'),exist('aw3_flag','var'),exist('ah3_aw3','var')])    
                   
       if any([all([ah3_flag == 1,aw3_flag == 1]),ah3_aw3  == 1])
           
           ah = ah3;
           
       elseif any([ah3_flag == 1,all([aw3_flag == 0,ah3_aw3 == 0])])
       
              ah = ah3;
          
       elseif  any([aw3_flag == 1,all([ah3_flag == 0,ah3_aw3 == 0])])  
       
               ah = aw3;
          
       end 
       
end
                      
%***CONVERT latD and lonD IN DEGREES TO RADIAN 
latR = (latD./180).*pi; %[radian]
lonR = (lonD./180).*pi ;%[radian]        

%GET SIZE OF USER INPUT TIME & POSITION
nrow_time=size(UTCtime,1);
nrow_pos=size(latD,1);

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
                
                %CONVERT ELEVATION ANGLE IN DEGREES TO RADIANS
                satEL =  abs(satELEV(j,i)) * pi / 180; 
                
                %CHECK FOR VIENNA MAPPING FUNCTION TYPE (VMF1 / VMF3)
                %********IF VMF1
                if any([exist('ah1_flag','var'),exist('aw1_flag','var'),exist('ah1_aw1','var')])
                    
                   if any([all([ah1_flag == 1,aw1_flag == 1]),ah1_aw1  == 1])
       
                      %COMPUTE BOTH HYDROSTATIC & WET COE USING vmf1 SUB-ROUTINE  
                      %Call the "vmf1.m" function
                      [VMFh(j,i),VMFw(j,i)] = vmf1(ah1(i),aw1(i),UTCtime(i,:),latR(i),hell(i),satEL);
                
                   elseif any([ah1_flag == 1,all([aw1_flag == 0,ah1_aw1 == 0])])
                          
                          %COMPUTE HYDROSTATIC COE USING vmf1 SUB-ROUTINE
                          %Call the "vmf1.m" function
                          [VMFh(j,i)] = vmf1(ah1(i),aw1(i),UTCtime(i,:),latR(i),hell(i),satEL);
                          
                   elseif any([aw1_flag == 1,all([ah1_flag == 0,ah1_aw1 == 0])])
                       
                          %COMPUTE WET COE USING vmf1 SUB-ROUTINE
                          %Call the "vmf1.m" function
                          [~,VMFw(j,i)] = vmf1(ah1(i),aw1(i),UTCtime(i,:),latR(i),hell(i),satEL);
                          
                   end %//if any([all([ah1_flag == 1,aw1_flag == 1]),ah1_aw1  == 1])
                   
                  
                %********IF VMF3  
                elseif any([exist('ah3_flag','var'),exist('aw3_flag','var'),exist('ah3_aw3','var')])    
                   
                       if any([all([ah3_flag == 1,aw3_flag == 1]),ah3_aw3  == 1])
                    
                          %COMPUTE BOTH HYDROSTATIC & WET COE USING vmf3 SUB-ROUTINE   
                          %Call the "vmf3.m" function
                          [VMFh(j,i),VMFw(j,i)] = vmf3(ah3(i),aw3(i),UTCtime(i,:),latR(i),lonR(i),hell(i),satEL);
                          
                       elseif any([ah3_flag == 1,all([aw3_flag == 0,ah3_aw3 == 0])])
                          
                              %COMPUTE HYDROSTATIC COE USING vmf3 SUB-ROUTINE
                              %Call the "vmf3.m" function
                              [VMFh(j,i)] = vmf3(ah3(i),aw3(i),UTCtime(i,:),latR(i),lonR(i),hell(i),satEL);
                              
                       elseif any([aw3_flag == 1,all([ah3_flag == 0,ah3_aw3 == 0])])
                           
                              %COMPUTE WET COE USING vmf3 SUB-ROUTINE
                              %Call the "vmf3.m" function
                              [~,VMFw(j,i)] = vmf3(ah3(i),aw3(i),UTCtime(i,:),latR(i),lonR(i),hell(i),satEL);
                             
                       end  
                
                end %//if any([exist('ah1_flag','var'),exist('aw1_flag','var'),exist('ah1_aw1','var')])
             
             end %//for j=1:nrow
             
          end %//for i=1:ncol 
          
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
      
                %CONVERT ELEVATION ANGLE IN DEGREES TO RADIANS
                satEL =  abs(satELEV(j,i)) * pi / 180; 
                
                %CHECK FOR VIENNA MAPPING FUNCTION TYPE (VMF1 / VMF3)
                %********IF VMF1
                if any([exist('ah1_flag','var'),exist('aw1_flag','var'),exist('ah1_aw1','var')])
                    
                   if any([all([ah1_flag == 1,aw1_flag == 1]),ah1_aw1  == 1])
       
                      %COMPUTE BOTH HYDROSTATIC & WET COE USING vmf1 SUB-ROUTINE  
                      %Call the "vmf1.m" function
                      [VMFh{j,i,k},VMFw{j,i,k}] = vmf1(ah1(i,k),aw1(i,k),UTCtime(k,:),latR(i),hell(i),satEL);
                
                   elseif any([ah1_flag == 1,all([aw1_flag == 0,ah1_aw1 == 0])])
                          
                          %COMPUTE HYDROSTATIC COE USING vmf1 SUB-ROUTINE
                          %Call the "vmf1.m" function
                          [VMFh{j,i,k}] = vmf1(ah1(i,k),aw1(i,k),UTCtime(k,:),latR(i),hell(i),satEL);
                          
                   elseif any([aw1_flag == 1,all([ah1_flag == 0,ah1_aw1 == 0])])
                       
                          %COMPUTE WET COE USING vmf1 SUB-ROUTINE
                          %Call the "vmf1.m" function
                          [~,VMFw{j,i,k}] = vmf1(ah1(i,k),aw1(i,k),UTCtime(k,:),latR(i),hell(i),satEL);
                          
                   end %//if any([all([ah1_flag == 1,aw1_flag == 1]),ah1_aw1  == 1])
                      
                        
                 %********IF VMF3
                elseif any([exist('ah3_flag','var'),exist('aw3_flag','var'),exist('ah3_aw3','var')])    
                   
                       if any([all([ah3_flag == 1,aw3_flag == 1]),ah3_aw3  == 1])
                    
                          %COMPUTE BOTH HYDROSTATIC & WET COE USING vmf3 SUB-ROUTINE   
                          %Call the "vmf3.m" function
                          [VMFh{j,i,k},VMFw{j,i,k}] = vmf3(ah3(i,k),aw3(i,k),UTCtime(k,:),latR(i),lonR(i),hell(i),satEL);
                          
                       elseif any([ah3_flag == 1,all([aw3_flag == 0,ah3_aw3 == 0])])
                          
                              %COMPUTE HYDROSTATIC COE USING vmf3 SUB-ROUTINE
                              %Call the "vmf3.m" function
                              [VMFh{j,i,k}] = vmf3(ah3(i,k),aw3(i,k),UTCtime(k,:),latR(i),lonR(i),hell(i),satEL);
                              
                       elseif any([aw3_flag == 1,all([ah3_flag == 0,ah3_aw3 == 0])])
                           
                              %COMPUTE WET COE USING vmf3 SUB-ROUTINE
                              %Call the "vmf3.m" function
                              [~,VMFw{j,i,k}] = vmf3(ah3(i,k),aw3(i,k),UTCtime(k,:),latR(i),lonR(i),hell(i),satEL);
                       
                       end %//if any([all([ah3_flag == 1,aw3_flag == 1]),ah3_aw3  == 1])   
                 
                end %//if any([exist('ah1_flag','var'),exist('aw1_flag','var'),exist('ah1_aw1','var')])
                
             end %//for j=1:nrow
             
         end %//for i=1:ncol
         
      end %//for k =1:ncol_ah
      
   end %//if isequal(Identical,1)
   
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
              
              %CONVERT ELEVATION ANGLE IN DEGREES TO RADIANS
              satEL =  abs(satELEV(j,i)) * pi / 180; %Elevation Angles converted to radian
              
              %CHECK FOR VIENNA MAPPING FUNCTION TYPE (VMF1 / VMF3)
              %********IF VMF1
              if any([exist('ah1_flag','var'),exist('aw1_flag','var'),exist('ah1_aw1','var')])
                    
                   if any([all([ah1_flag == 1,aw1_flag == 1]),ah1_aw1  == 1])
       
                      %COMPUTE BOTH HYDROSTATIC & WET COE USING vmf1 SUB-ROUTINE  
                      %Call the "vmf1.m" function
                      [VMFh{j,i,k},VMFw{j,i,k}] = vmf1(ah1(i,k),aw1(i,k),UTCtime(k,:),latR(i),hell(i),satEL);
                
                   elseif any([ah1_flag == 1,all([aw1_flag == 0,ah1_aw1 == 0])])
                          
                          %COMPUTE HYDROSTATIC COE USING vmf1 SUB-ROUTINE
                          %Call the "vmf1.m" function
                          [VMFh{j,i,k}] = vmf1(ah1(i,k),aw1(i,k),UTCtime(k,:),latR(i),hell(i),satEL);
                          
                   elseif any([aw1_flag == 1,all([ah1_flag == 0,ah1_aw1 == 0])])
                       
                          %COMPUTE WET COE USING vmf1 SUB-ROUTINE
                          %Call the "vmf1.m" function
                          [~,VMFw{j,i,k}] = vmf1(ah1(i,k),aw1(i,k),UTCtime(k,:),latR(i),hell(i),satEL);
                          
                   end %//if any([all([ah1_flag == 1,aw1_flag == 1]),ah1_aw1  == 1])
                      
                        
                 %********IF VMF3
              elseif  any([exist('ah3_flag','var'),exist('aw3_flag','var'),exist('ah3_aw3','var')])    
                   
                       if any([all([ah3_flag == 1,aw3_flag == 1]),ah3_aw3  == 1])
                    
                          %COMPUTE BOTH HYDROSTATIC & WET COE USING vmf3 SUB-ROUTINE   
                          %Call the "vmf3.m" function
                          [VMFh{j,i,k},VMFw{j,i,k}] = vmf3(ah3(i,k),aw3(i,k),UTCtime(k,:),latR(i),lonR(i),hell(i),satEL);
                          
                       elseif any([ah3_flag == 1,all([aw3_flag == 0,ah3_aw3 == 0])])
                          
                              %COMPUTE HYDROSTATIC COE USING vmf3 SUB-ROUTINE
                              %Call the "vmf3.m" function
                              [VMFh{j,i,k}] = vmf3(ah3(i,k),aw3(i,k),UTCtime(k,:),latR(i),lonR(i),hell(i),satEL);
                              
                       elseif any([aw3_flag == 1,all([ah3_flag == 0,ah3_aw3 == 0])])
                           
                              %COMPUTE WET COE USING vmf3 SUB-ROUTINE
                              %Call the "vmf3.m" function
                              [~,VMFw{j,i,k}] = vmf3(ah3(i,k),aw3(i,k),UTCtime(k,:),latR(i),lonR(i),hell(i),satEL);
                       
                       end %//if any([all([ah3_flag == 1,aw3_flag == 1]),ah3_aw3  == 1])   
                 
              end  %//if any([exist('ah1_flag','var'),exist('aw1_flag','var'),exist('ah1_aw1','var')])
                
          end %//for j=1:nrow
             
          end %//for i=1:ncol
         
       end %//for k =1:ncol_ah
      
       
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
      
              %CONVERT ELEVATION ANGLE IN DEGREES TO RADIANS
              satEL =  abs(satELEV(j,i)) * pi / 180; %Elevation Angles converted to radian 
             
             %CHECK FOR VIENNA MAPPING FUNCTION TYPE (VMF1 / VMF3)
             %********IF VMF1
             if any([exist('ah1_flag','var'),exist('aw1_flag','var'),exist('ah1_aw1','var')])
                    
                   if any([all([ah1_flag == 1,aw1_flag == 1]),ah1_aw1  == 1])
       
                      %COMPUTE BOTH HYDROSTATIC & WET COE USING vmf1 SUB-ROUTINE  
                      %Call the "vmf1.m" function
                      [VMFh(j,i),VMFw(j,i)] = vmf1(ah1(i),aw1(i),UTCtime(i,:),latR(i),hell(i),satEL);
                
                   elseif any([ah1_flag == 1,all([aw1_flag == 0,ah1_aw1 == 0])])
                          
                          %COMPUTE HYDROSTATIC COE USING vmf1 SUB-ROUTINE
                          %Call the "vmf1.m" function
                          [VMFh(j,i)] = vmf1(ah1(i),aw1(i),UTCtime(i,:),latR(i),hell(i),satEL);
                          
                   elseif any([aw1_flag == 1,all([ah1_flag == 0,ah1_aw1 == 0])])
                       
                          %COMPUTE WET COE USING vmf1 SUB-ROUTINE
                          %Call the "vmf1.m" function
                          [~,VMFw(j,i)] = vmf1(ah1(i),aw1(i),UTCtime(i,:),latR(i),hell(i),satEL);
                          
                   end %//if any([all([ah1_flag == 1,aw1_flag == 1]),ah1_aw1  == 1])
                   
                  
                %********IF VMF3  
                elseif any([exist('ah3_flag','var'),exist('aw3_flag','var'),exist('ah3_aw3','var')])    
                   
                       if any([all([ah3_flag == 1,aw3_flag == 1]),ah3_aw3  == 1])
                    
                          %COMPUTE BOTH HYDROSTATIC & WET COE USING vmf3 SUB-ROUTINE   
                          %Call the "vmf3.m" function
                          [VMFh(j,i),VMFw(j,i)] = vmf3(ah3(i),aw3(i),UTCtime(i,:),latR(i),lonR(i),hell(i),satEL);
                          
                       elseif any([ah3_flag == 1,all([aw3_flag == 0,ah3_aw3 == 0])])
                          
                              %COMPUTE HYDROSTATIC COE USING vmf3 SUB-ROUTINE
                              %Call the "vmf3.m" function
                              [VMFh(j,i)] = vmf3(ah3(i),aw3(i),UTCtime(i,:),latR(i),lonR(i),hell(i),satEL);
                              
                       elseif any([aw3_flag == 1,all([ah3_flag == 0,ah3_aw3 == 0])])
                           
                              %COMPUTE WET COE USING vmf3 SUB-ROUTINE
                              %Call the "vmf3.m" function
                              [~,VMFw(j,i)] = vmf3(ah3(i),aw3(i),UTCtime(i,:),latR(i),lonR(i),hell(i),satEL);
                             
                       end  
                
             end  %//if any([exist('ah1_flag','var'),exist('aw1_flag','var'),exist('ah1_aw1','var')])
             
           end  %//for j=1:nrow
             
        end  %//for i=1:ncol 
        
     end %//if ncol_ah > 1
                
end %//if isequal(nrow_time,nrow_pos)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF VMF.m  %%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%C.1.SUBROUTINE TO COMPUTE VIENNA MAPPING FUNCTION1(VMF1)
function [vmf1h,vmf1w] = vmf1(ah,aw,UTCtime,lat,hgt,satELEV)

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
vmf1h         = vmf1h + ht_corr; 
 
%(2)*********COMPUTE WET MAPPING FUNCTION(VMFw)
%DEFINE COEFFICIENTs
bw = 0.00146;
cw = 0.04391;
beta   = bw/( sine + cw );
gamma  = aw/( sine + beta);
topcon = (1.d0 + aw/(1.d0 + bw/(1.d0 + cw)));
vmf1w   = topcon/(sine+gamma);  
%******************************************END OF vmf1.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%C.2.SUBROUTINE TO COMPUTE VIENNA MAPPING FUNCTION 3(VMF3)
function [vmf3h,vmf3w] = vmf3(ah,aw,UTCtime,lat,lon,h_ell,satELEV)

%**************************************************************************
%***DESCRIPTION:
%              This subroutine determines the Tropospheric Mapping         *
%              Function using the Vienna Mapping Functions 3(VMF3) for     * 
%              Specific site given receiver/station position in[lat lon hgt*
%              Satellites Elevation Angle(el),and receiver reception time  *
%              in the utc time format[Yr Mon Day Hr Min Sec]                                      * 
%              The a coefficients have to be inserted from discrete data,  *
%              while the b and c coefficients are of empirical nature      *
%              containing a geographical and temporal dependence,          *
%              represented  in spherical harmonics.The spherical harmonics *
%              coefficients are developed to degree and order 12 and are   *
%              based on a 5°x5° grid containing ray-tracing data from      *
%              2001-2010. All input quantities have to be scalars!         *
%USAGE:                                                                    *
%      General:                                                            *
%              [VMFh3,VMFw3] = vmf3(ah,aw,UTCtime,lat,lon,hgt,satELEV)     *       
%****INPUTs:                                                               *
%1.         ah: hydrostatic MF coefficient a (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)
%2.         aw: wet MF coefficient a         (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)                           *
%3.    UTCtime: UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *
%4.        lat: station ellipsoidal latitude in [radians]                  *
%          lon: station ellipsoidal longitude (radians)                    *
%5.      h_ell: station ellipsoidal height in [meters]                     *
%6.    satELEV: Satellites Elevation Angle in [radian]                    *
% 
%***OUTPUTs:                                                               *
%         VMFh:Vienna Hydrostatic mapping function 3                       *
%         VMFw:Vienna Wet mapping function  3                              *
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
[~, mjd,~]=utc2JulianDay_DoY(Yr,Mn,Day,H,MIN,SECs);
    
%*********COMPUTE DAY OF YEAR(doy)
% first check if the specified year is leap year or not (logical output)
leapYear = ((mod(Yr,4) == 0 & mod(Yr,100) ~= 0) | mod(Yr,400) == 0);

days = [31 28 31 30 31 30 31 31 30 31 30 31];
doy = sum(days(1:Mn-1)) + Day;
if leapYear == 1 && Mn > 2
    doy = doy + 1;
end
doy = doy + mjd-floor(mjd);   % add decimal places

    
%LEGENDRE FUNCTIONS FOR for bh, bw, ch and cw 
anm_bh = [0.00271285863109945 -1.39197786008938e-06 1.34955672002719e-06 2.71686279717968e-07 1.56659301773925e-06;9.80476624811974e-06 -5.83922611260673e-05 -2.07307023860417e-05 1.14628726961148e-06 4.93610283608719e-06;-1.03443106534268e-05 -2.05536138785961e-06 2.09692641914244e-06 -1.55491034130965e-08 -1.89706404675801e-07;-3.00353961749658e-05 2.37284447073503e-05 2.02236885378918e-05 1.69276006349609e-06 8.72156681243892e-07;-7.99121077044035e-07 -5.39048313389504e-06 -4.21234502039861e-06 -2.70944149806894e-06 -6.80894455531746e-07;7.51439609883296e-07 3.85509708865520e-07 4.41508016098164e-08 -2.07507808307757e-08 4.95354985050743e-08;2.21790962160087e-05 -5.56986238775212e-05 -1.81287885563308e-05 -4.41076013532589e-06 4.93573223917278e-06;-4.47639989737328e-06 -2.60452893072120e-06 2.56376320011189e-06 4.41600992220479e-07 2.93437730332869e-07;8.14992682244945e-07 2.03945571424434e-07 1.11832498659806e-08 3.25756664234497e-08 3.01029040414968e-08;-7.96927680907488e-08 -3.66953150925865e-08 -6.74742632186619e-09 -1.30315731273651e-08 -2.00748924306947e-09;-2.16138375166934e-05 1.67350317962556e-05 1.93768260076821e-05 1.99595120161850e-06 -2.42463528222014e-06;5.34360283708044e-07 -3.64189022040600e-06 -2.99935375194279e-06 -2.06880962903922e-06 -9.40815692626002e-07;6.80235884441822e-07 1.33023436079845e-07 -1.80349593705226e-08 2.51276252565192e-08 -1.43240592002794e-09;-7.13790897253802e-08 7.81998506267559e-09 1.13826909570178e-09 -5.89629600214654e-09 -4.20760865522804e-09;-5.80109372399116e-09 1.13702284491976e-09 7.29046067602764e-10 -9.10468988754012e-10 -2.58814364808642e-10;1.75558618192965e-05 -2.85579168876063e-05 -1.47442190284602e-05 -6.29300414335248e-06 -5.12204538913460e-07;-1.90788558291310e-06 -1.62144845155361e-06 7.57239241641566e-07 6.93365788711348e-07 6.88855644570695e-07;2.27050351488552e-07 1.03925791277660e-07 -3.31105076632079e-09 2.88065761026675e-08 -8.00256848229136e-09;-2.77028851807614e-08 -5.96251132206930e-09 2.95987495527251e-10 -5.87644249625625e-09 -3.28803981542337e-09;-1.89918479865558e-08 3.54083436578857e-09 8.10617835854935e-10 4.99207055948336e-10 -1.52691648387663e-10;1.04022499586096e-09 -2.36437143845013e-10 -2.25110813484842e-10 -7.39850069252329e-11 7.95929405440911e-11;-3.11579421267630e-05 -3.43576336877494e-06 5.81663608263384e-06 8.31534700351802e-07 4.02619520312154e-06;6.00037066879001e-07 -1.12538760056168e-07 -3.86745332115590e-07 -3.88218746020826e-07 -6.83764967176388e-07;-9.79583981249316e-08 9.14964449851003e-08 4.77779838549237e-09 2.44283811750703e-09 -6.26361079345158e-09;-2.37742207548109e-08 -5.53336301671633e-09 -3.73625445257115e-09 -1.92304189572886e-09 -7.18681390197449e-09;-6.58203463929583e-09 9.28456148541896e-10 2.47218904311077e-10 1.10664919110218e-10 -4.20390976974043e-11;9.45857603373426e-10 -3.29683402990254e-11 -8.15440375865127e-11 -1.21615589356628e-12 -9.70713008848085e-12;1.61377382316176e-10 6.84326027598147e-12 -4.66898885683671e-12 2.31211355085535e-12 2.39195112937346e-12;2.99634365075821e-07 8.14391615472128e-06 6.70458490942443e-06 -9.92542646762000e-07 -3.04078064992750e-06;-6.52697933801393e-07 2.87255329776428e-07 -1.78227609772085e-08 2.65525429849935e-07 8.60650570551813e-08;-1.62727164011710e-07 1.09102479325892e-07 4.97827431850001e-09 7.86649963082937e-11 -6.67193813407656e-09;-2.96370000987760e-09 1.20008401576557e-09 1.75885448022883e-09 -1.74756709684384e-09 3.21963061454248e-09;-9.91101697778560e-10 7.54541713140752e-10 -2.95880967800875e-10 1.81009160501278e-10 8.31547411640954e-11;1.21268051949609e-10 -5.93572774509587e-11 -5.03295034994351e-11 3.05383430975252e-11 3.56280438509939e-11;6.92012970333794e-11 -9.02885345797597e-12 -3.44151832744880e-12 2.03164894681921e-12 -5.44852265137606e-12;5.56731263672800e-12 3.57272150106101e-12 2.25885622368678e-12 -2.44508240047675e-13 -6.83314378535235e-13;3.96883487797254e-06 -4.57100506169608e-06 -3.30208117813256e-06 3.32599719134845e-06 4.26539325549339e-06;1.10123151770973e-06 4.58046760144882e-07 1.86831972581926e-07 -1.60092770735081e-07 -5.58956114867062e-07;-3.40344900506653e-08 2.87649741373047e-08 -1.83929753066251e-08 -9.74179203885847e-09 -2.42064137485043e-09;-6.49731596932566e-09 -3.07048108404447e-09 -2.84380614669848e-09 1.55123146524283e-09 4.53694984588346e-10;5.45175793803325e-10 -3.73287624700125e-10 -1.16293122618336e-10 7.25845618602690e-11 -4.34112440021627e-11;1.89481447552805e-10 3.67431482211078e-12 -1.72180065021194e-11 1.47046319023226e-11 1.31920481414062e-11;2.10125915737167e-12 -3.08420783495975e-12 -4.87748712363020e-12 1.16363599902490e-14 1.26698255558605e-13;-8.07894928696254e-12 9.19344620512607e-13 3.26929173307443e-13 2.00438149416495e-13 -9.57035765212079e-15;1.38737151773284e-12 1.09340178371420e-13 5.15714202449053e-14 -5.92156438588931e-14 -3.29586752336143e-14;6.38137197198254e-06 4.62426300749908e-06 4.42334454191034e-06 1.15374736092349e-06 -2.61859702227253e-06;-2.25320619636149e-07 3.21907705479353e-07 -3.34834530764823e-07 -4.82132753601810e-07 -3.22410936343355e-07;3.48894515496995e-09 3.49951261408458e-08 -6.01128959281142e-09 4.78213900943443e-09 1.46012816168576e-08;-9.66682871952083e-11 3.75806627535317e-09 2.38984004956705e-09 2.07545049877203e-09 1.58573595632766e-09;1.06834370693917e-09 -4.07975055112153e-10 -2.37598937943957e-10 5.89327007480137e-11 1.18891820437634e-10;5.22433722695807e-11 6.02011995016293e-12 -7.80605402956048e-12 1.50873145627341e-11 -1.40550093106311e-12;2.13396242187279e-13 -1.71939313965536e-12 -3.57625378660975e-14 -5.01675184988446e-14 -1.07805487368797e-12;-1.24352330043311e-12 8.26105883301606e-13 4.63606970128517e-13 6.39517888984486e-14 -7.35135439920086e-14;-5.39023859065631e-13 2.54188315588243e-14 1.30933833278664e-14 6.06153473304781e-15 -4.24722717533726e-14;3.12767756884813e-14 -2.29517847871632e-15 2.53117304424948e-16 7.07504914138118e-16 -1.20089065310688e-15;2.08311178819214e-06 -1.22179185044174e-06 -2.98842190131044e-06 3.07310218974299e-06 2.27100346036619e-06;-3.94601643855452e-07 -5.44014825116083e-07 -6.16955333162507e-08 -2.31954821580670e-07 1.14010813005310e-07;6.11067575043044e-08 -3.93240193194272e-08 -1.62979132528933e-08 1.01339204652581e-08 1.97319601566071e-08;2.57770508710055e-09 1.87799543582899e-09 1.95407654714372e-09 1.15276419281270e-09 2.25397005402120e-09;7.16926338026236e-10 -3.65857693313858e-10 -1.54864067050915e-11 6.50770211276549e-11 -7.85160007413546e-12;4.90007693914221e-12 3.31649396536340e-12 4.81664871165640e-13 7.26080745617085e-12 2.30960953372164e-12;9.75489202240545e-13 -1.68967954531421e-13 7.38383391334110e-13 -3.58435515913239e-13 -3.01564710027450e-13;-3.79533601922805e-13 2.76681830946617e-13 1.21480375553803e-13 -1.57729077644850e-14 -8.87664977818700e-14;-3.96462845480288e-14 2.94155690934610e-14 6.78413205760717e-15 -4.12135802787361e-15 -1.46373307795619e-14;-8.64941937408121e-15 -1.91822620970386e-15 -8.01725413560744e-16 5.02941051180784e-16 -1.07572628474344e-15;-4.13816294742758e-15 -7.43602019785880e-17 -5.54248556346072e-17 -4.83999456005158e-17 -1.19622559730466e-16;-8.34852132750364e-07 -7.45794677612056e-06 -6.58132648865533e-06 -1.38608110346732e-06 5.32326534882584e-07;-2.75513802414150e-07 3.64713745106279e-08 -7.12385417940442e-08 -7.86206067228882e-08 2.28048393207161e-08;-4.26696415431918e-08 -4.65599668635087e-09 7.35037936327566e-09 1.17098354115804e-08 1.44594777658035e-08;1.12407689274199e-09 7.62142529563709e-10 -6.72563708415472e-10 -1.18094592485992e-10 -1.17043815733292e-09;1.76612225246125e-10 -1.01188552503192e-10 7.32546072616968e-11 1.79542821801610e-11 -2.23264859965402e-11;-9.35960722512375e-12 1.90894283812231e-12 -6.34792824525760e-13 3.98597963877826e-12 -4.47591409078971e-12;-3.34623858556099e-12 4.56384903915853e-14 2.72561108521416e-13 -3.57942733300468e-15 1.99794810657713e-13;-6.16775522568954e-14 8.25316968328823e-14 7.19845814260518e-14 -2.92415710855106e-14 -5.49570017444031e-15;-8.50728802453217e-15 8.38161600916267e-15 3.43651657459983e-15 -8.19429434115910e-16 -4.08905746461100e-15;4.39042894275548e-15 -3.69440485320477e-16 1.22249256876779e-16 -2.09359444520984e-16 -3.34211740264257e-16;-5.36054548134225e-16 3.29794204041989e-17 2.13564354374585e-17 -1.37838993720865e-18 -1.29188342867753e-17;-3.26421841529845e-17 7.38235405234126e-18 2.49291659676210e-18 8.18252735459593e-19 1.73824952279230e-20;4.67237509268208e-06 1.93611283787239e-06 9.39035455627622e-07 -5.84565118072823e-07 -1.76198705802101e-07;-3.33739157421993e-07 4.12139555299163e-07 1.58754695700856e-07 1.37448753329669e-07 1.04722936936873e-07;6.64200603076386e-09 1.45412222625734e-08 1.82498796118030e-08 2.86633517581614e-09 1.06066984548100e-09;5.25549696746655e-09 -1.33677183394083e-09 7.60804375937931e-11 -1.07918624219037e-10 8.09178898247941e-10;1.89318454110039e-10 9.23092164791765e-11 5.51434573131180e-11 3.86696392289240e-11 -1.15208165047149e-11;-1.02252706006226e-12 -7.25921015411136e-13 -1.98110126887620e-12 -2.18964868282672e-13 -7.18834476685625e-13;-2.69770025318548e-12 -2.17850340796321e-14 4.73040820865871e-13 1.57947421572149e-13 1.86925164972766e-13;1.07831718354771e-13 2.26681841611017e-14 2.56046087047783e-14 -1.14995851659554e-14 -2.27056907624485e-14;6.29825154734712e-15 8.04458225889001e-16 9.53173540411138e-16 1.16892301877735e-15 -1.04324684545047e-15;-5.57345639727027e-16 -2.93949227634932e-16 7.47621406284534e-18 -5.36416885470756e-17 -2.87213280230513e-16;1.73219775047208e-16 2.05017387523061e-17 9.08873886345587e-18 -2.86881547225742e-18 -1.25303645304992e-17;-7.30829109684568e-18 2.03711261415353e-18 7.62162636124024e-19 -7.54847922012517e-19 -8.85105098195030e-19;5.62039968280587e-18 -1.38144206573507e-19 1.68028711767211e-20 1.81223858251981e-19 -8.50245194985878e-20];
anm_bw = [0.00136127467401223 -6.83476317823061e-07 -1.37211986707674e-06 7.02561866200582e-07 -2.16342338010651e-07;-9.53197486400299e-06 6.58703762338336e-06 2.42000663952044e-06 -6.04283463108935e-07 2.02144424676990e-07;-6.76728911259359e-06 6.03830755085583e-07 -8.72568628835897e-08 2.21750344140938e-06 1.05146032931020e-06;-3.21102832397338e-05 -7.88685357568093e-06 -2.55495673641049e-06 -1.99601934456719e-06 -4.62005252198027e-07;-7.84639263523250e-07 3.11624739733849e-06 9.02170019697389e-07 6.37066632506008e-07 -9.44485038780872e-09;2.19476873575507e-06 -2.20580510638233e-07 6.94761415598378e-07 4.80770865279717e-07 -1.34357837196401e-07;2.18469215148328e-05 -1.80674174262038e-06 -1.52754285605060e-06 -3.51212288219241e-07 2.73741237656351e-06;2.85579058479116e-06 1.57201369332361e-07 -2.80599072875081e-07 -4.91267304946072e-07 -2.11648188821805e-07;2.81729255594770e-06 3.02487362536122e-07 -1.64836481475431e-07 -2.11607615408593e-07 -6.47817762225366e-08;1.31809947620223e-07 -1.58289524114549e-07 -7.05580919885505e-08 5.56781440550867e-08 1.23403290710365e-08;-1.29252282695869e-05 -1.07247072037590e-05 -3.31109519638196e-06 2.13776673779736e-06 -1.49519398373391e-07;1.81685152305722e-06 -1.17362204417861e-06 -3.19205277136370e-08 4.09166457255416e-07 1.53286667406152e-07;1.63477723125362e-06 -2.68584775517243e-08 4.94662064805191e-09 -7.09027987928288e-08 4.44353430574937e-08;-2.13090618917978e-07 4.05836983493219e-08 2.94495876336549e-08 -1.75005469063176e-08 -3.03015988647002e-09;-2.16074435298006e-09 9.37631708987675e-09 -2.05996036369828e-08 6.97068002894092e-09 -8.90988987979604e-09;1.38047798906967e-05 2.05528261553901e-05 1.59072148872708e-05 7.34088731264443e-07 1.28226710383580e-06;7.08175753966264e-07 -9.27988276636505e-07 1.60535820026081e-07 -3.27296675122065e-07 -2.20518321170684e-07;1.90932483086199e-07 -7.44215272759193e-08 1.81330673333187e-08 4.37149649043616e-08 4.18884335594172e-08;-5.37009063880924e-08 2.22870057779431e-08 1.73740123037651e-08 -4.45137302235032e-09 9.44721910524571e-09;-6.83406949047909e-08 -1.95046676795923e-10 2.57535903049686e-09 4.82643164083020e-09 3.37657333705158e-09;3.96128688448981e-09 -6.63809403270686e-10 2.44781464212534e-10 5.92280853590699e-11 -4.78502591970721e-10;1.75859399041414e-05 -2.81238050668481e-06 -2.43670534594848e-06 3.58244562699714e-06 -1.76547446732691e-06;-1.06451311473304e-07 1.54336689617184e-06 -2.00690000442673e-07 1.38790047911880e-09 -1.62490619890017e-07;-2.72757421686155e-07 1.71139266205398e-07 -2.55080309401917e-08 -8.40793079489831e-09 -1.01129447760167e-08;2.92966025844079e-08 -2.07556718857313e-08 5.45985315647905e-09 8.76857690274150e-09 1.06785510440474e-08;-1.22059608941331e-08 6.52491630264276e-09 -1.79332492326928e-10 3.75921793745396e-10 -7.06416506254786e-10;1.63224355776652e-09 4.95586028736232e-10 -3.07879011759040e-10 -7.78354087544277e-11 1.43959047067250e-10;3.86319414653663e-10 -2.06467134617933e-10 4.37330971382694e-11 -5.00421056263711e-11 -9.40237773015723e-12;-1.23856142706451e-05 7.61047394008415e-06 -1.99104114578138e-07 6.86177748886858e-07 -1.09466747592827e-07;2.99866062403128e-07 1.87525561397390e-07 4.99374806994715e-08 4.86229763781404e-07 4.46570575517658e-07;-5.05748332368430e-07 1.95523624722285e-08 -9.17535435911345e-08 -2.56671607433547e-08 -7.11896201616653e-08;-2.66062200406494e-08 -5.40470019739274e-09 -2.29718660244954e-09 -3.73328592264404e-09 3.38748313712376e-09;5.30855327954894e-10 5.28851845648032e-10 -2.22278913745418e-10 -5.52628653064771e-11 -9.24825145219684e-10;6.03737227573716e-10 -3.52190673510919e-12 -1.30371720641414e-10 -9.12787239944822e-12 6.42187285537238e-12;1.78081862458539e-10 2.93772078656037e-12 -1.04698379945322e-11 -2.82260024833024e-11 -5.61810459067525e-12;9.35003092299580e-12 -8.23133834521577e-13 5.54878414224198e-13 -3.62943215777181e-13 2.38858933771653e-12;-1.31216096107331e-05 -5.70451670731759e-06 -5.11598683573971e-06 -4.99990779887599e-06 1.27389320221511e-07;-1.23108260369048e-06 5.53093245213587e-07 8.60093183929302e-07 2.65569700925696e-07 1.95485134805575e-07;-2.29647072638049e-07 -5.45266515081825e-08 2.85298129762263e-08 1.98167939680185e-08 5.52227340898335e-09;-2.73844745019857e-08 -4.48345173291362e-10 -1.93967347049382e-09 -1.41508853776629e-09 -1.75456962391145e-09;-2.68863184376108e-11 -2.20546981683293e-09 6.56116990576877e-10 1.27129855674922e-10 -2.32334506413213e-10;1.98303136881156e-10 6.04782006047075e-11 2.91291115431570e-11 6.18098615782757e-11 -3.82682292530379e-11;9.48294455071158e-12 -3.05873596453015e-13 5.31539408055057e-13 -7.31016438665600e-12 -1.19921002209198e-11;-2.25188050845725e-11 -3.91627574966393e-13 -6.80217235976769e-13 5.91033607278405e-13 5.02991534452191e-13;1.29532063896247e-12 1.66337285851564e-13 3.25543028344555e-13 1.89143357962363e-13 3.32288378169726e-13;-2.45864358781728e-06 4.49460524898260e-06 1.03890496648813e-06 -2.73783420376785e-06 7.12695730642593e-07;-9.27805078535168e-07 -4.97733876686731e-07 9.18680298906510e-08 -2.47200617423980e-07 6.16163630140379e-08;-1.39623661883136e-08 -1.12580495666505e-07 2.61821435950379e-08 -2.31875562002885e-08 5.72679835033659e-08;-9.52538983318497e-09 -5.40909215302433e-09 1.88698793952475e-09 -4.08127746406372e-09 1.09534895853812e-10;3.79767457525741e-09 1.11549801373366e-10 -6.45504957274111e-10 3.05477141010356e-10 1.26261210565856e-10;5.08813577945300e-11 1.43250547678637e-11 8.81616572082448e-12 2.58968878880804e-11 3.83421818249954e-11;8.95094368142044e-12 -3.26220304555971e-12 -1.28047847191896e-12 2.67562170258942e-12 2.72195031576670e-12;-6.47181697409757e-12 1.13776457455685e-12 2.84856274334969e-13 -7.63667272085395e-14 -1.34451657758826e-13;-1.25291265888343e-12 8.63500441050317e-14 -1.21307856635548e-13 5.12570529540511e-14 3.32389276976573e-14;3.73573418085813e-14 -5.37808783042784e-16 -4.23430408270850e-16 -4.75110565740493e-15 6.02553212780166e-15;8.95483987262751e-06 -3.90778212666235e-06 -1.12115019808259e-06 1.78678942093383e-06 1.46806344157962e-06;-4.59185232678613e-07 1.09497995905419e-07 1.31663977640045e-07 4.20525791073626e-08 -9.71470741607431e-08;1.63399802579572e-07 1.50909360648645e-08 -1.11480472593347e-08 -1.84000857674573e-08 7.82124614794256e-09;1.22887452385094e-08 -4.06647399822746e-10 -6.49120327585597e-10 8.63651225791194e-10 -2.73440085913102e-09;2.51748630889583e-09 4.79895880425564e-10 -2.44908073860844e-10 2.56735882664876e-10 -1.64815306286912e-10;4.85671381736718e-11 -2.51742732115131e-11 -2.60819437993179e-11 6.12728324086123e-12 2.16833310896138e-11;4.11389702320298e-12 -8.09433180989935e-13 -1.19812498226024e-12 1.46885737888520e-12 3.15807685137836e-12;-1.47614580597013e-12 4.66726413909320e-13 1.72089709006255e-13 1.13854935381418e-13 2.77741161317003e-13;-1.02257724967727e-13 1.10394382923502e-13 -3.14153505370805e-15 2.41103099110106e-14 2.13853053149771e-14;-3.19080885842786e-14 -9.53904307973447e-15 2.74542788156379e-15 2.33797859107844e-15 -2.53192474907304e-15;-5.87702222126367e-15 -1.80133850930249e-15 -3.09793125614454e-16 -1.04197538975295e-16 3.72781664701327e-16;1.86187054729085e-06 8.33098045333428e-06 3.18277735484232e-06 -7.68273797022231e-07 -1.52337222261696e-06;-5.07076646593648e-07 -8.61959553442156e-07 -3.51690005432816e-07 -4.20797082902431e-07 -3.07652993252673e-07;-7.38992472164147e-08 -8.39473083080280e-08 -2.51587083298935e-08 7.30691259725451e-09 -3.19457155958983e-08;-1.99777182012924e-09 -3.21265085916022e-09 -4.84477421865675e-10 -1.82924814205799e-09 -3.46664344655997e-10;-7.05788559634927e-11 1.21840735569025e-10 7.97347726425926e-11 1.08275679614409e-10 -1.17891254809785e-10;1.10299718947774e-11 -3.22958261390263e-11 -1.43535798209229e-11 6.87096504209595e-12 -6.64963212272352e-12;-6.47393639740084e-12 1.03156978325120e-12 -9.20099775082358e-14 -2.40150316641949e-13 1.14008812047857e-12;-1.23957846397250e-13 2.85996703969692e-13 1.91579874982553e-13 5.20597174693064e-14 -4.06741434883370e-14;-2.35479068911236e-14 1.97847338186993e-14 1.58935977518516e-15 -2.32217195254742e-15 -8.48611789490575e-15;1.03992320391626e-14 1.54017082092642e-15 1.05950035082788e-16 -1.17870898461353e-15 -1.10937420707372e-15;-1.09011948374520e-15 -6.04168007633584e-16 -9.10901998157436e-17 1.98379116989461e-16 -1.03715496658498e-16;-1.38171942108278e-16 -6.33037999097522e-17 -1.38777695011470e-17 1.94191397045401e-17 5.70055906754485e-18;1.92989406002085e-06 -3.82662130483128e-06 -4.60189561036048e-07 2.24290587856309e-06 1.40544379451550e-06;6.49033717633394e-08 2.41396114435326e-07 2.73948898223321e-07 1.10633664439332e-07 -3.19555270171075e-08;-2.91988966963297e-08 -6.03828192816571e-09 1.18462386444840e-08 1.32095545004128e-08 -5.06572721528914e-09;7.31079058474148e-09 -8.42775299751834e-10 1.10190810090667e-09 1.96592273424306e-09 -2.13135932785688e-09;7.06656405314388e-11 1.43441125783756e-10 1.46962246686924e-10 7.44592776425197e-11 -3.64331892799173e-11;-2.52393942119372e-11 1.07520964869263e-11 5.84669886072094e-12 6.52029744217103e-12 1.82947123132059e-12;-4.15669940115121e-12 -1.95963254053648e-13 2.16977822834301e-13 -2.84701408462031e-13 4.27194601040231e-13;3.07891105454129e-13 1.91523190672955e-13 1.05367297580989e-13 -5.28136363920236e-14 -3.53364110005917e-14;7.02156663274738e-15 9.52230536780849e-15 -3.41019408682733e-15 -3.59825303352899e-15 -2.62576411636150e-15;-1.75110277413804e-15 5.29265220719483e-16 4.45015980897919e-16 -3.80179856341347e-16 -4.32917763829695e-16;1.16038609651443e-16 -6.69643574373352e-17 2.65667154817303e-17 -9.76010333683956e-17 4.07312981076655e-17;5.72659246346386e-18 1.30357528108671e-18 2.49193258417535e-18 1.76247014075584e-18 7.59614374197688e-19;1.03352170833303e-17 -2.30633516638829e-18 2.84777940620193e-18 -7.72161347944693e-19 6.07028034506380e-19];
anm_ch = [0.0571481238161787 3.35402081801137e-05 3.15988141788728e-05 -1.34477341887086e-05 -2.61831023577773e-07;5.77367395845715e-05 -0.000669057185209558 -6.51057691648904e-05 -1.61830149147091e-06 8.96771209464758e-05;-8.50773002452907e-05 -4.87106614880272e-05 4.03431160775277e-05 2.54090162741464e-06 -5.59109319864264e-06;0.00150536423187709 0.000611682258892697 0.000369730024614855 -1.95658439780282e-05 -3.46246726553700e-05;-2.32168718433966e-05 -0.000127478686553809 -9.00292451740728e-05 -6.07834315901830e-05 -1.04628419422714e-05;-1.38607250922551e-06 -3.97271603842309e-06 -8.16155320152118e-07 5.73266706046665e-07 2.00366060212696e-07;6.52491559188663e-05 -0.00112224323460183 -0.000344967958304075 -7.67282640947300e-05 0.000107907110551939;-0.000138870461448036 -7.29995695401936e-05 5.35986591445824e-05 9.03804869703890e-06 8.61370129482732e-06;-9.98524443968768e-07 -6.84966792665998e-08 1.47478021860771e-07 1.94857794008064e-06 7.17176852732910e-07;1.27066367911720e-06 1.12113289164288e-06 2.71525688515375e-07 -2.76125723009239e-07 -1.05429690305013e-07;-0.000377264999981652 0.000262691217024294 0.000183639785837590 3.93177048515576e-06 -6.66187081899168e-06;-4.93720951871921e-05 -0.000102820030405771 -5.69904376301748e-05 -3.79603438055116e-05 -3.96726017834930e-06;-2.21881958961135e-06 -1.40207117987894e-06 1.60956630798516e-07 2.06121145135022e-06 6.50944708093149e-07;2.21876332411271e-07 1.92272880430386e-07 -6.44016558013941e-09 -1.40954921332410e-07 -4.26742169137667e-07;-3.51738525149881e-08 2.89616194332516e-08 -3.40343352397886e-08 -2.89763392721812e-08 -6.40980581663785e-10;3.51240856823468e-05 -0.000725895015345786 -0.000322514037108045 -0.000106143759981636 4.08153152459337e-05;-2.36269716929413e-05 -4.20691836557932e-05 1.43926743222922e-05 2.61811210631784e-05 2.09610762194903e-05;-7.91765756673890e-07 1.64556789159745e-06 -9.43930166276555e-07 6.46641738736139e-07 -5.91509547299176e-07;3.92768838766879e-07 -1.98027731703690e-07 -5.41303590057253e-08 -4.21705797874207e-07 -6.06042329660681e-08;-1.56650141024305e-08 7.61808165752027e-08 -1.81900460250934e-08 1.30196216971675e-08 1.08616031342379e-08;-2.80964779829242e-08 -7.25951488826103e-09 -2.59789823306225e-09 -2.79271942407154e-09 4.10558774868586e-09;-0.000638227857648286 -0.000154814045363391 7.78518327501759e-05 -2.95961469342381e-05 1.15965225055757e-06;4.47833146915112e-06 1.33712284237555e-05 3.61048816552123e-06 -2.50717844073547e-06 -1.28100822021734e-05;-2.26958070007455e-06 2.57779960912242e-06 1.08395653197976e-06 1.29403393862805e-07 -1.04854652812567e-06;-3.98954043463392e-07 -2.26931182815454e-07 -1.09169545045028e-07 -1.49509536031939e-07 -3.98376793949903e-07;2.30418911071110e-08 1.23098508481555e-08 -1.71161401463708e-08 2.35829696577657e-09 1.31136164162040e-08;3.69423793101582e-09 3.49231027561927e-10 -1.18581468768647e-09 5.43180735828820e-10 5.43192337651588e-10;-1.38608847117992e-09 -1.86719145546559e-10 -8.13477384765498e-10 2.01919878240491e-10 1.00067892622287e-10;-4.35499078415956e-05 0.000450727967957804 0.000328978494268850 -3.05249478582848e-05 -3.21914834544310e-05;1.24887940973241e-05 1.34275239548403e-05 1.11275518344713e-06 7.46733554562851e-06 -2.12458664760353e-06;9.50250784948476e-07 2.34367372695203e-06 -5.43099244798980e-07 -4.35196904508734e-07 -8.31852234345897e-07;5.91775478636535e-09 -1.48970922508592e-07 2.99840061173840e-08 -1.30595933407792e-07 1.27136765045597e-07;-1.78491083554475e-08 1.76864919393085e-08 -1.96740493482011e-08 1.21096708004261e-08 2.95518703155064e-10;1.75053510088658e-09 -1.31414287871615e-09 -1.44689439791928e-09 1.14682483668460e-09 1.74488616540169e-09;1.08152964586251e-09 -3.85678162063266e-10 -2.77851016629979e-10 3.89890578625590e-11 -2.54627365853495e-10;-1.88340955578221e-10 5.19645384002867e-11 2.14131326027631e-11 1.24027770392728e-11 -9.42818962431967e-12;0.000359777729843898 -0.000111692619996219 -6.87103418744904e-05 0.000115128973879551 7.59796247722486e-05;5.23717968000879e-05 1.32279078116467e-05 -5.72277317139479e-07 -7.56326558610214e-06 -1.95749622214651e-05;1.00109213210139e-06 -2.75515216592735e-07 -1.13393194050846e-06 -4.75049734870663e-07 -3.21499480530932e-07;-2.07013716598890e-07 -7.31392258077707e-08 -3.96445714084160e-08 3.21390452929387e-08 -1.43738764991525e-08;2.03081434931767e-09 -1.35423687136122e-08 -4.47637454261816e-09 2.18409121726643e-09 -3.74845286805217e-09;3.17469255318367e-09 2.44221027314129e-10 -2.46820614760019e-10 7.55851003884434e-10 6.98980592550891e-10;9.89541493531067e-11 -2.78762878057315e-11 -2.10947962916771e-10 3.77882267360636e-11 -1.20009542671532e-12;5.01720575730940e-11 1.66470417102135e-11 -7.50624817938091e-12 9.97880221482238e-12 4.87141864438892e-12;2.53137945301589e-11 1.93030083090772e-12 -1.44708804231290e-12 -1.77837100743423e-12 -8.10068935490951e-13;0.000115735341520738 0.000116910591048350 8.36315620479475e-05 1.61095702669207e-05 -7.53084853489862e-05;-9.76879433427199e-06 9.16968438003335e-06 -8.72755127288830e-06 -1.30077933880053e-05 -9.78841937993320e-06;1.04902782517565e-07 2.14036988364936e-07 -7.19358686652888e-07 1.12529592946332e-07 7.07316352860448e-07;7.63177265285080e-08 1.22781974434290e-07 8.99971272969286e-08 5.63482239352990e-08 4.31054352285547e-08;3.29855763107355e-09 -6.95004336734441e-09 -6.52491370576354e-09 1.97749180391742e-09 3.51941791940498e-09;3.85373745846559e-10 1.65754130924183e-10 -3.31326088103057e-10 5.93256024580436e-10 1.27725220636915e-10;-1.08840956376565e-10 -4.56042860268189e-11 -4.77254322645633e-12 -2.94405398621875e-12 -3.07199979999475e-11;2.07389879095010e-11 1.51186798732451e-11 9.28139802941848e-12 5.92738269687687e-12 9.70337402306505e-13;-2.85879708060306e-12 1.92164314717053e-13 4.02664678967890e-14 5.18246319204277e-13 -7.91438726419423e-13;6.91890667590734e-13 -8.49442290988352e-14 -5.54404947212402e-15 9.71093377538790e-15 -5.33714333415971e-14;-5.06132972789792e-05 -4.28348772058883e-05 -6.90746551020305e-05 8.48380415176836e-05 7.04135614675053e-05;-1.27945598849788e-05 -1.92362865537803e-05 -2.30971771867138e-06 -8.98515975724166e-06 5.25675205004752e-06;-8.71907027470177e-07 -1.02091512861164e-06 -1.69548051683864e-07 4.87239045855761e-07 9.13163249899837e-07;-6.23651943425918e-08 6.98993315829649e-08 5.91597766733390e-08 4.36227124230661e-08 6.45321798431575e-08;-1.46315079552637e-10 -7.85142670184337e-09 1.48788168857903e-09 2.16870499912160e-09 -1.16723047065545e-09;3.31888494450352e-10 1.90931898336457e-10 -3.13671901557599e-11 2.60711798190524e-10 8.45240112207997e-11;1.36645682588537e-11 -5.68830303783976e-12 1.57518923848140e-11 -1.61935794656758e-11 -4.16568077748351e-12;9.44684950971905e-13 7.30313977131995e-12 3.14451447892684e-12 6.49029875639842e-13 -9.66911019905919e-13;-8.13097374090024e-13 5.23351897822186e-13 8.94349188113951e-14 -1.33327759673270e-13 -4.04549450989029e-13;-3.76176467005839e-14 -6.19953702289713e-14 -3.74537190139726e-14 1.71275486301958e-14 -3.81946773167132e-14;-4.81393385544160e-14 3.66084990006325e-15 3.10432030972253e-15 -4.10964475657416e-15 -6.58644244242900e-15;-7.81077363746945e-05 -0.000254773632197303 -0.000214538508009518 -3.80780934346726e-05 1.83495359193990e-05;5.89140224113144e-06 -3.17312632433258e-06 -3.81872516710791e-06 -2.27592226861647e-06 1.57044619888023e-06;-1.44272505088690e-06 -1.10236588903758e-07 2.64336813084693e-07 4.76074163332460e-07 4.28623587694570e-07;3.98889120733904e-08 -1.29638005554027e-08 -4.13668481273828e-08 1.27686793719542e-09 -3.54202962042383e-08;1.60726837551750e-09 -2.70750776726156e-09 2.79387092681070e-09 -3.01419734793998e-10 -1.29101669438296e-10;-2.55708290234943e-10 2.27878015173471e-11 -6.43063443462716e-12 1.26531554846856e-10 -1.65822147437220e-10;-3.35886470557484e-11 -3.51895009091595e-12 5.80698399963198e-12 -2.84881487149207e-12 8.91708061745902e-12;-3.12788523950588e-12 3.35366912964637e-12 2.52236848033838e-12 -8.12801050709184e-13 -2.63510394773892e-13;6.83791881183142e-14 2.41583263270381e-13 8.58807794189356e-14 -5.12528492761045e-14 -1.40961725631276e-13;-1.28585349115321e-14 -2.11049721804969e-14 5.26409596614749e-15 -4.31736582588616e-15 -1.60991602619068e-14;-9.35623261461309e-15 -3.94384886372442e-16 5.04633016896942e-16 -5.40268998456055e-16 -1.07857944298104e-15;8.79756791888023e-16 4.52529935675330e-16 1.36886341163227e-16 -1.12984402980452e-16 6.30354561057224e-18;0.000117829256884757 2.67013591698442e-05 2.57913446775250e-05 -4.40766244878807e-05 -1.60651761172523e-06;-1.87058092029105e-05 1.34371169060024e-05 5.59131416451555e-06 4.50960364635647e-06 2.87612873904633e-06;2.79835536517287e-07 8.93092708148293e-07 8.37294601021795e-07 -1.99029785860896e-08 -8.87240405168977e-08;4.95854313394905e-08 -1.44694570735912e-08 2.51662229339375e-08 -3.87086600452258e-09 2.29741919071270e-08;4.71497840986162e-09 2.47509999454076e-09 1.67323845102824e-09 8.14196768283530e-10 -3.71467396944165e-10;-1.07340743907054e-10 -8.07691657949326e-11 -5.99381660248133e-11 2.33173929639378e-12 -2.26994195544563e-11;-3.83130441984224e-11 -5.82499946138714e-12 1.43286311435124e-11 3.15150503353387e-12 5.97891025146774e-12;-5.64389191072230e-13 9.57258316335954e-13 1.12055192185939e-12 -4.42417706775420e-13 -9.93190361616481e-13;1.78188860269677e-13 7.82582024904950e-14 5.18061650118009e-14 2.13456507353387e-14 -5.26202113779510e-14;-8.18481324740893e-15 -3.71256746886786e-15 4.23508855164371e-16 -2.91292502923102e-15 -1.15454205389350e-14;6.16578691696810e-15 6.74087154080877e-16 5.71628946437034e-16 -2.05251213979975e-16 -7.25999138903781e-16;9.35481959699383e-17 6.23535830498083e-17 3.18076728802060e-18 -2.92353209354587e-17 7.65216088665263e-19;2.34173078531701e-17 -8.30342420281772e-18 -4.33602329912952e-18 1.90226281379981e-18 -7.85507922718903e-19];
anm_cw = [0.0395329695826997 -0.000131114380761895 -0.000116331009006233 6.23548420410646e-05 5.72641113425116e-05;-0.000441837640880650 0.000701288648654908 0.000338489802858270 3.76700309908602e-05 -8.70889013574699e-06;1.30418530496887e-05 -0.000185046547597376 4.31032103066723e-05 0.000105583334124319 3.23045436993589e-05;3.68918433448519e-05 -0.000219433014681503 3.46768613485000e-06 -9.17185187163528e-05 -3.69243242456081e-05;-6.50227201116778e-06 2.07614874282187e-05 -5.09131314798362e-05 -3.08053225174359e-05 -4.18483655873918e-05;2.67879176459056e-05 -6.89303730743691e-05 2.11046783217168e-06 1.93163912538178e-05 -1.97877143887704e-06;0.000393937595007422 -0.000452948381236406 -0.000136517846073846 0.000138239247989489 0.000133175232977863;5.00214539435002e-05 3.57229726719727e-05 -9.38010547535432e-07 -3.52586798317563e-05 -7.01218677681254e-06;3.91965314099929e-05 1.02236686806489e-05 -1.95710695226022e-05 -5.93904795230695e-06 3.24339769876093e-06;6.68158778290653e-06 -8.10468752307024e-06 -9.91192994096109e-06 -1.89755520007723e-07 -3.26799467595579e-06;0.000314196817753895 -0.000296548447162009 -0.000218410153263575 -1.57318389871000e-05 4.69789570185785e-05;0.000104597721123977 -3.31000119089319e-05 5.60326793626348e-05 4.71895007710715e-05 3.57432326236664e-05;8.95483021572039e-06 1.44019305383365e-05 4.87912790492931e-06 -3.45826387853503e-06 3.23960320438157e-06;-1.35249651009930e-05 -2.49349762695977e-06 -2.51509483521132e-06 -9.14254874104858e-07 -8.57897406100890e-07;-1.68143325235195e-06 1.72073417594235e-06 1.38765993969565e-06 4.09770982137530e-07 -6.60908742097123e-07;-0.000639889366487161 0.00120194042474696 0.000753258598887703 3.87356377414663e-05 1.31231811175345e-05;2.77062763606783e-05 -9.51425270178477e-06 -6.61068056107547e-06 -1.38713669012109e-05 9.84662092961671e-06;-2.69398078539471e-06 6.50860676783123e-06 3.80855926988090e-06 -1.98076068364785e-06 1.17187335666772e-06;-2.63719028151905e-06 5.03149473656743e-07 7.38964893399716e-07 -8.38892485369078e-07 1.30943917775613e-06;-1.56634992245479e-06 -2.97026487417045e-08 5.06602801102463e-08 -4.60436007958792e-08 -1.62536449440997e-07;-2.37493912770935e-07 1.69781593069938e-08 8.35178275224265e-08 -4.83564044549811e-08 -4.96448864199318e-08;0.00134012259587597 -0.000250989369253194 -2.97647945512547e-05 -6.47889968094926e-05 8.41302130716859e-05;-0.000113287184900929 4.78918993866293e-05 -3.14572113583139e-05 -2.10518256626847e-05 -2.03933633847417e-05;-4.97413321312139e-07 3.72599822034753e-06 -3.53221588399266e-06 -1.05232048036416e-06 -2.74821498198519e-06;4.81988542428155e-06 4.21400219782474e-07 1.02814808667637e-06 4.40299068486188e-09 3.37103399036634e-09;1.10140301678818e-08 1.90257670180182e-07 -1.00831353341885e-08 1.44860642389714e-08 -5.29882089987747e-08;6.12420414245775e-08 -4.48953461152996e-09 -1.38837603709003e-08 -2.05533675904779e-08 1.49517908802329e-09;9.17090243673643e-10 -9.24878857867367e-09 -2.30856560363943e-09 -4.36348789716735e-09 -4.45808881183025e-10;-0.000424912699609112 -0.000114365438471564 -0.000403200981827193 4.19949560550194e-05 -3.02068483713739e-05;3.85435472851225e-05 -5.70726887668306e-05 4.96313706308613e-07 1.02395703617082e-05 5.85550000567006e-06;-7.38204470183331e-06 -4.56638770109511e-06 -3.94007992121367e-06 -2.16666812189101e-06 -4.55694264113194e-06;5.89841165408527e-07 1.40862905173449e-08 1.08149086563211e-07 -2.18592601537944e-07 -3.78927431428119e-07;4.85164687450468e-08 8.34273921293655e-08 1.47489605513673e-08 6.01494125001291e-08 6.43812884159484e-09;1.13055580655363e-08 3.50568765400469e-09 -5.09396162501750e-09 -1.83362063152411e-09 -4.11227251553035e-09;3.16454132867156e-09 -1.39634794131087e-09 -7.34085003895929e-10 -7.55541371271796e-10 -1.57568747643705e-10;1.27572900992112e-09 -3.51625955080441e-10 -4.84132020565098e-10 1.52427274930711e-10 1.27466120431317e-10;-0.000481655666236529 -0.000245423313903835 -0.000239499902816719 -0.000157132947351028 5.54583099258017e-05;-1.52987254785589e-05 2.78383892116245e-05 4.32299123991860e-05 1.70981319744327e-05 -1.35090841769225e-06;-8.65400907717798e-06 -6.51882656990376e-06 -2.43810171017369e-07 8.54348785752623e-07 2.98371863248143e-07;-1.68155571776752e-06 -3.53602587563318e-07 -1.00404435881759e-07 -2.14162249012859e-08 -2.42131535531526e-07;-1.08048603277187e-08 -9.78850785763030e-08 -2.32906554437417e-08 2.22003630858805e-08 -2.27230368089683e-09;-5.98864391551041e-09 7.38970926486848e-09 3.61322835311957e-09 3.70037329172919e-09 -3.41121137081362e-09;-7.33113754909726e-10 -9.08374249335220e-11 -1.78204392133739e-10 8.28618491929026e-11 -1.32966817912373e-10;-5.23340481314676e-10 1.36403528233346e-10 -7.04478837151279e-11 -6.83175201536443e-12 -2.86040864071134e-12;3.75347503578356e-11 -1.08518134138781e-11 -2.53583751744508e-12 1.00168232812303e-11 1.74929602713312e-11;-0.000686805336370570 0.000591849814585706 0.000475117378328026 -2.59339398048415e-05 3.74825110514968e-05;3.35231363034093e-05 2.38331521146909e-05 7.43545963794093e-06 -3.41430817541849e-06 7.20180957675353e-06;3.60564374432978e-07 -3.13300039589662e-06 -6.38974746108020e-07 -8.63985524672024e-07 2.43367665208655e-06;-4.09605238516094e-07 -2.51158699554904e-07 -1.29359217235188e-07 -2.27744642483133e-07 7.04065989970205e-08;6.74886341820129e-08 -1.02009407061935e-08 -3.30790296448812e-08 1.64959795655031e-08 1.40641779998855e-08;1.31706886235108e-09 -1.06243701278671e-09 -2.85573799673944e-09 3.72566568681289e-09 2.48402582003925e-09;-3.68427463251097e-11 -1.90028122983781e-10 -3.98586561768697e-11 1.14458831693287e-11 -2.27722300377854e-12;-7.90029729611056e-11 3.81213646526419e-11 4.63303426711788e-11 1.52294835905903e-11 -2.99094751490726e-12;-2.36146602045017e-11 1.03852674709985e-11 -4.47242126307100e-12 5.30884113537806e-12 1.68499023262969e-12;-3.30107358134527e-13 -4.73989085379655e-13 5.17199549822684e-13 2.34951744478255e-13 2.05931351608192e-13;0.000430215687511780 -0.000132831373000014 -3.41830835017045e-05 4.70312161436033e-06 -3.84807179340006e-05;1.66861163032403e-05 -8.10092908523550e-06 8.20658107437905e-06 6.12399025026683e-06 -1.85536495631911e-06;1.53552093641337e-06 2.19486495660361e-06 -1.07253805120137e-06 -4.72141767909137e-07 4.00744581573216e-07;2.56647305130757e-07 -8.07492046592274e-08 -2.05858469296168e-07 1.09784168930599e-07 -7.76823030181225e-08;1.77744008115031e-08 1.64134677817420e-08 4.86163044879020e-09 1.13334251800856e-08 -7.17260621115426e-09;1.61133063219326e-09 -1.85414677057024e-09 -2.13798537812651e-09 1.15255123229679e-09 2.24504700129464e-09;1.23344223096739e-10 -1.20385012169848e-10 -2.18038256346433e-12 3.23033120628279e-11 8.01179568213400e-11;-6.55745274387847e-12 1.22127104697198e-11 5.83805016355883e-12 -8.31201582509817e-12 1.90985373872656e-12;-2.89199983667265e-12 5.05962500506667e-12 1.28092925110279e-12 5.60353813743813e-13 1.76753731968770e-12;-1.61678729774956e-13 -3.92206170988615e-13 -9.04941327579237e-14 1.89847694200763e-13 4.10008676756463e-14;-1.16808369005656e-13 -9.97464591430510e-14 7.46366550245722e-15 2.53398578153179e-14 1.06510689748906e-14;-0.000113716921384790 -0.000131902722651488 -0.000162844886485788 7.90171538739454e-06 -0.000178768066961413;-2.13146535366500e-06 -3.57818705543597e-05 -1.50825855069298e-05 -2.17909259570022e-05 -8.19332236308581e-06;-2.88001138617357e-06 -2.09957465440793e-06 6.81466526687552e-08 3.58308906974448e-07 -4.18502067223724e-07;-1.10761444317605e-07 6.91773860777929e-08 8.17125372450372e-08 -2.16476237959181e-08 7.59221970502074e-08;-9.56994224818941e-09 6.64104921728432e-09 6.33077902928348e-09 2.85721181743727e-09 -6.39666681678123e-09;4.62558627839842e-10 -1.69014863754621e-09 -2.80260429599733e-10 4.27558937623863e-11 -1.66926133269027e-10;-7.23385132663753e-11 5.51961193545280e-11 3.04070791942335e-11 3.23227055919062e-12 8.47312431934829e-11;-1.61189613765486e-11 1.66868155925172e-11 1.05370341694715e-11 -4.41495859079592e-12 -2.24939051401750e-12;-8.72229568056267e-13 1.88613726203286e-12 1.21711137534390e-14 -1.13342372297867e-12 -6.87151975256052e-13;7.99311988544090e-15 4.46150979586709e-14 7.50406779454998e-14 -3.20385428942275e-14 -1.26543636054393e-14;4.80503817699514e-14 -3.35545623603729e-14 -1.18546423610485e-14 4.19419209985980e-15 -1.73525614436880e-14;-1.20464898830163e-15 -8.80752065000456e-16 -1.22214298993313e-15 1.69928513019657e-15 1.93593051311405e-16;1.68528879784841e-05 3.57144412031081e-05 -1.65999910125077e-05 5.40370336805755e-05 0.000118138122851376;-3.28151779115881e-05 1.04231790790798e-05 -2.80761862890640e-06 2.98996152515593e-06 -2.67641158709985e-06;-2.08664816151978e-06 -1.64463884697475e-06 6.79099429284834e-08 7.23955842946495e-07 -6.86378427465657e-07;-2.88205823027255e-09 2.38319699493291e-09 1.14169347509045e-07 8.12981074994402e-08 -1.56957943666988e-07;-7.09711403570189e-09 6.29470515502988e-09 3.50833306577579e-09 8.31289199649054e-09 -2.14221463168338e-09;-8.11910123910038e-10 3.34047829618955e-10 3.70619377446490e-10 3.30426088213373e-10 4.86297305597865e-11;1.98628160424161e-11 -4.98557831380098e-12 -5.90523187802174e-12 -1.27027116925122e-12 1.49982368570355e-11;2.62289263262748e-12 3.91242360693861e-12 6.56035499387192e-12 -1.17412941089401e-12 -9.40878197853394e-13;-3.37805010124487e-13 5.39454874299593e-13 -2.41569839991525e-13 -2.41572016820792e-13 -3.01983673057198e-13;-1.85034053857964e-13 4.31132161871815e-14 4.13497222026824e-15 -4.60075514595980e-14 -1.92454846400146e-14;2.96113888929854e-15 -1.11688534391626e-14 3.76275373238932e-15 -3.72593295948136e-15 1.98205490249604e-16;1.40074667864629e-15 -5.15564234798333e-16 3.56287382196512e-16 5.07242777691587e-16 -2.30405782826134e-17;2.96822530176851e-16 -4.77029898301223e-17 1.12782285532775e-16 1.58443229778573e-18 8.22141904662969e-17];
bnm_bh = [0 0 0 0 0;0 0 0 0 0;-2.29210587053658e-06 -2.33805004374529e-06 -7.49312880102168e-07 -5.12022747852006e-07 5.88926055066172e-07;0 0 0 0 0;-4.63382754843690e-06 -2.23853015662938e-06 8.14830531656518e-07 1.15453269407116e-06 -4.53555450927571e-07;-6.92432096320778e-07 -2.98734455136141e-07 1.48085153955641e-08 1.37881746148773e-07 -6.92492118460215e-09;0 0 0 0 0;-1.91507979850310e-06 -1.83614825459598e-06 -7.46807436870647e-07 -1.28329122348007e-06 5.04937180063059e-07;-8.07527103916713e-07 2.83997840574570e-08 -6.01890498063025e-08 -2.48339507554546e-08 2.46284627824308e-08;-2.82995069303093e-07 1.38818274596408e-09 3.22731214161408e-09 2.87731153972404e-10 1.53895537278496e-08;0 0 0 0 0;-6.68210270956800e-07 -2.19104833297845e-06 1.30116691657253e-07 4.78445730433450e-07 -4.40344300914051e-07;-2.36946755740436e-07 -1.32730991878204e-07 1.83669593693860e-08 7.90218931983569e-08 -4.70161979232584e-08;1.07746083292179e-07 -4.17088637760330e-09 -1.83296035841109e-09 -5.80243971371211e-09 -2.11682361167439e-09;-5.44712355496109e-08 1.89717032256923e-09 2.27327316287804e-10 7.78400728280038e-10 8.82380487618991e-12;0 0 0 0 0;-5.61707049615673e-08 -1.09066447089585e-06 -2.25742250174119e-07 -8.64367795924377e-07 1.06411275240680e-08;2.41782935157918e-08 -3.65762298303819e-08 -6.93420659586875e-08 -3.97316214341991e-08 -2.08767816486390e-08;6.38293030383436e-08 1.11377936334470e-08 6.91424941454782e-09 1.39887159955004e-09 5.25428749022906e-09;1.09291268489958e-08 1.23935926756516e-10 3.92917259954515e-10 -1.79144682483562e-10 -9.11802874917597e-10;-4.40957607823325e-09 1.45751390560667e-10 1.24641258165301e-10 -6.45810339804674e-11 -8.92894658893326e-12;0 0 0 0 0;1.54754294162102e-08 -1.60154742388847e-06 -4.08425188394881e-07 6.18170290113531e-09 -2.58919765162122e-07;1.37130642286873e-08 -6.67813955828458e-08 -7.01410996605609e-09 3.82732572660461e-08 -2.73381870915135e-08;2.19113155379218e-08 4.11027496396868e-09 6.33816020485226e-09 -1.49242411327524e-09 -6.14224941851705e-10;6.26573021218961e-09 5.17137416480052e-10 -3.49784328298676e-10 1.13578756343208e-10 2.80414613398411e-10;1.65048133258794e-11 1.00047239417239e-10 1.05124654878499e-10 -3.03826002621926e-11 4.57155388334682e-11;6.20221691418381e-11 9.75852610098156e-12 -5.46716005756984e-12 1.31643349569537e-11 3.61618775715470e-12;0 0 0 0 0;-1.03938913012708e-06 -1.78417431315664e-07 2.86040141364439e-07 1.83508599345952e-08 -1.34452220464346e-07;-4.36557481393662e-08 7.49780206868834e-09 -8.62829428674082e-09 5.50577793039009e-09 -9.46897502333254e-09;3.43193738406672e-10 1.13545447306468e-08 1.25242388852214e-09 6.03221501959620e-10 1.57172070361180e-09;-4.73307591021391e-10 1.70855824051391e-10 -2.62470421477037e-11 2.04525835988874e-10 -1.17859695928164e-10;-3.36185995299839e-10 3.19243054562183e-11 1.17589412418126e-10 -1.35478747434514e-12 5.11192214558542e-11;3.19640547592136e-11 2.94297823804643e-12 -1.00651526276990e-11 -1.67028733953153e-12 3.03938833625503e-12;1.68928641118173e-11 -7.90032886682002e-13 -1.40899773539137e-12 7.76937592393354e-13 7.32539820298651e-13;0 0 0 0 0;2.32949756055277e-07 1.46237594908093e-07 -1.07770884952484e-07 1.26824870644476e-07 -2.36345735961108e-08;8.89572676497766e-08 7.24810004121931e-08 2.67583556180119e-08 2.48434796111361e-08 -3.55004782858686e-09;-1.00823909773603e-08 8.84433929029076e-10 -2.55502517594511e-10 -5.48034274059119e-10 -8.50241938494079e-10;1.13259819566467e-09 5.55186945221216e-10 7.63679807785295e-11 -1.70067998092043e-11 1.57081965572493e-10;-2.37748192185353e-10 2.45463764948000e-11 3.23208414802860e-11 -2.72624834520723e-12 8.14449183666500e-12;-1.54977633126025e-11 4.58754903157884e-12 -1.25864665839074e-12 2.44139868157872e-12 -1.82827441958193e-12;3.28285563794513e-12 -1.10072329225465e-12 -7.23470501810935e-13 5.85309745620389e-13 4.11317589687125e-13;4.57596974384170e-13 9.84198128213558e-14 3.34503817702830e-14 7.08431086558307e-15 2.79891177268807e-14;0 0 0 0 0;-3.67820719155580e-07 6.98497901205902e-07 1.83397388750300e-07 2.39730262495372e-07 -2.58441984368194e-07;5.17793954077994e-08 5.54614175977835e-08 1.75026214305232e-09 -2.55518450411346e-09 -6.12272723006537e-09;-7.94292648157198e-09 -1.01709107852895e-09 -1.49251241812310e-09 9.32827213605682e-10 -8.24490722043118e-10;1.36410408475679e-11 2.16390220454971e-10 1.24934806872235e-10 -6.82507825145903e-11 -4.01575177719668e-11;-1.41619917600555e-11 -1.54733230409082e-11 1.36792829351538e-11 1.11157862104733e-12 2.08548465892268e-11;-3.56521723755846e-12 4.47877185884557e-12 -6.34096209274637e-16 -1.13010624512348e-12 -2.82018136861041e-13;2.22758955943441e-12 -4.63876465559380e-13 -5.80688019272507e-13 2.45878690598655e-13 1.49997666808106e-13;-6.26833903786958e-14 2.73416335780807e-14 1.91842340758425e-14 1.67405061129010e-14 -2.45268543953704e-17;1.81972870222228e-14 5.43036245069085e-15 1.92476637107321e-15 8.78498602508626e-17 -1.42581647227657e-15;0 0 0 0 0;9.74322164613392e-07 -5.23101820582724e-07 -2.81997898176227e-07 4.54762451707384e-08 -3.34645078118827e-08;-6.75813194549663e-09 3.49744702199583e-08 -5.09170419895883e-09 5.24359476874755e-09 4.96664262534662e-09;4.53858847892396e-10 -1.49347392165963e-09 -2.00939511362154e-09 9.30987163387955e-10 9.74450200826854e-11;-4.92900885858693e-10 5.34223033225688e-12 1.08501839729368e-10 -6.43526142089173e-11 -3.11063319142619e-11;1.38469246386690e-11 -7.91180584906922e-12 2.26641656746936e-13 4.55251515177956e-12 6.05270575117769e-12;4.02247935664225e-12 1.82776657951829e-12 -1.28348801405445e-13 -2.16257301300350e-13 -5.54363979435025e-14;4.15005914461687e-13 -2.00647573581168e-13 -1.67278251942946e-13 1.30332398257985e-13 1.52742363652434e-13;6.36376500056974e-14 1.65794532815776e-14 -3.80832559052662e-15 -6.40262894005341e-16 2.42577181848072e-15;-5.55273521249151e-15 3.69725182221479e-15 2.02114207545759e-15 -4.50870833392161e-16 9.62950493696677e-17;1.00935904205024e-17 6.54751873609395e-17 -1.09138810997186e-16 -8.62396750098759e-17 -3.82788257844306e-17;0 0 0 0 0;4.21958510903678e-07 -8.30678271007705e-08 -3.47006439555247e-07 -3.36442823712421e-08 9.90739768222027e-08;2.64389033612742e-08 2.65825090066479e-09 -1.28895513428522e-08 -7.07182694980098e-10 7.10907165301180e-09;6.31203524153492e-09 -1.67038260990134e-09 1.33104703539822e-09 8.34376495185149e-10 -2.52478613522612e-10;1.18414896299279e-10 -2.57745052288455e-11 2.88295935685818e-11 -3.27782977418354e-11 -1.05705000036156e-11;-4.20826459055091e-12 -6.97430607432268e-12 -3.90660545970607e-12 -3.90449239948755e-13 -4.60384797517466e-13;-9.47668356558200e-13 6.53305025354881e-13 2.63240185434960e-13 1.40129115015734e-13 3.85788887132074e-14;2.23947810407291e-13 7.35262771548253e-15 -3.83348211931292e-14 4.20376514344176e-14 4.26445836468461e-14;-3.88008154470596e-16 2.28561424667750e-15 -8.73599966653373e-16 2.14321147947665e-15 6.38631825071920e-16;-8.62165565535721e-15 1.79742912149810e-15 1.01541125038661e-15 -7.91027655831866e-17 -4.06505132825230e-16;-2.35355054392189e-16 -6.13997759731013e-17 -2.73490528665965e-17 2.63895177155121e-17 -4.47531057245187e-18;6.01909706823530e-17 5.35520010856833e-18 -2.15530106132531e-18 -2.46778496746231e-18 -7.09947296442799e-19;0 0 0 0 0;-3.75005956318736e-07 -5.39872297906819e-07 -1.19929654883034e-07 4.52771083775007e-08 1.82790552943564e-07;7.82606642505646e-09 -1.68890832383153e-08 -8.45995188378997e-09 1.42958730598502e-09 3.21075754133531e-09;4.28818421913782e-09 -1.07501469928219e-09 8.84086350297418e-10 9.74171228764155e-10 8.59877149602304e-12;1.28983712172521e-10 -6.96375160373676e-11 -2.13481436408896e-11 1.33516375568179e-11 -1.65864626508258e-11;-4.48914384622368e-12 9.68953616831263e-13 -1.61372463422897e-12 -2.09683563440448e-12 -1.90096826314068e-12;-1.12626619779175e-13 3.34903159106509e-14 -1.21721528343657e-13 7.46246339290354e-14 3.68424909859186e-13;5.08294274367790e-14 2.83036159977090e-14 1.48074873486387e-14 -9.59633528834945e-15 -1.26231060951100e-14;-4.01464098583541e-16 1.97047929526674e-15 -5.29967950447497e-16 -3.59120406619931e-16 1.69690933982683e-16;-1.73919209873841e-15 7.52792462841274e-16 3.65589287101147e-16 -7.79247612043812e-17 -8.24599670368999e-17;-4.61555616150128e-17 4.94529746019753e-19 -1.09858157212270e-17 3.95550811124928e-18 3.23972399884100e-18;-2.27040686655766e-17 -3.27855689001215e-18 -3.30649011116861e-19 9.08748546536849e-19 8.92197599890994e-19;5.67241944733762e-18 3.84449400209976e-19 1.77668058015537e-19 2.00432838283455e-20 -2.00801461564767e-19];
bnm_bw = [0 0 0 0 0;0 0 0 0 0;-9.56715196386889e-06 -3.68040633020420e-08 1.27846786489883e-07 1.32525487755973e-06 1.53075361125066e-06;0 0 0 0 0;-7.17682617983607e-06 2.89994188119445e-06 -2.97763578173405e-07 8.95742089134942e-07 3.44416325304006e-07;-8.02661132285210e-07 3.66738692077244e-07 -3.02880965723280e-07 3.54144282036103e-07 -1.68873066391463e-07;0 0 0 0 0;-2.89640569283461e-06 -7.83566373343614e-07 -8.36667214682577e-07 -7.41891843549121e-07 -9.23922655636489e-08;-1.06144662284862e-06 1.57709930505924e-07 1.04203025714319e-07 1.20783300488461e-07 -1.38726055821134e-07;-4.16549018672265e-07 -1.35220897698872e-07 -6.40269964829901e-08 1.63258283210837e-08 -2.57958025095959e-08;0 0 0 0 0;3.52324885892419e-06 -2.26705543513814e-07 1.53835589488292e-06 -3.75263061267433e-07 3.69384057396017e-07;-2.06569149157664e-07 -9.36260183227175e-08 -3.55985284353048e-08 -9.13671163891094e-08 6.93156256562600e-09;1.32437594740782e-07 4.44349887272663e-08 -3.38192451721674e-08 -3.97263855781102e-08 -1.93087822995800e-09;-1.29595244818942e-07 -1.40852985547683e-08 1.42587592939760e-09 7.05779876554001e-09 -1.00996269264535e-08;0 0 0 0 0;4.06960756215938e-06 -1.97898540226986e-06 7.21905857553588e-08 -1.19908881538755e-06 -5.67561861536903e-08;6.53369660286999e-08 -2.42818687866392e-07 -1.66203004559493e-08 -2.41512414151897e-08 4.45426333411018e-08;1.44650670663281e-07 8.50666367433859e-09 -4.61165612004307e-09 4.88527987491045e-09 1.06277326713172e-08;1.86770937103513e-08 -6.44197940288930e-10 -7.60456736846174e-09 -9.97186468682689e-10 8.73229752697716e-10;-1.00206566229113e-08 1.33934372663121e-09 1.41691503439220e-09 8.72352590578753e-10 -8.04561626629829e-10;0 0 0 0 0;3.07161843116618e-06 1.82962085656470e-06 1.87728623016069e-07 7.10611617623261e-07 2.26499092250481e-07;4.50766403064905e-08 -1.67752393078256e-07 2.47844723639070e-08 -3.56484348424869e-09 -1.56634836636584e-08;3.77011651881090e-08 -7.23045828480496e-09 5.22995988863761e-09 -1.03740320341306e-09 4.57839777217789e-09;8.09495635883121e-09 -3.01977244420529e-10 -2.30104544933093e-09 3.63658580939428e-10 4.39320811714867e-10;9.37087629961269e-11 1.00780920426635e-09 1.28140539913350e-10 -6.65795285522138e-12 4.71732796198631e-11;-8.88504487069155e-11 -1.63253810435461e-10 7.22669710644299e-11 5.64715132584527e-11 -1.08949308197617e-12;0 0 0 0 0;-2.64054293284174e-07 -2.37611606117256e-06 -1.83671059706264e-06 -3.12199354841993e-07 -1.05598289276114e-07;7.41706968747147e-08 -1.64359098062646e-08 -3.09750224040234e-08 -9.68640079410317e-09 -7.90399057863403e-08;-1.00254376564271e-08 1.12528248631191e-08 -2.67841549174100e-09 -2.69481819323647e-09 1.56550607475331e-09;-2.18568129350729e-09 6.26422056977450e-10 1.95007291427316e-09 3.14226463591125e-10 -3.62000388344482e-10;-9.30451291747549e-10 5.62175549482704e-11 1.01022849902012e-10 5.18675856498499e-11 5.37561696283235e-11;5.33151334468794e-11 1.07571307336725e-10 -1.31714567944652e-11 -4.17524405900018e-11 -2.16737797893502e-12;4.69916869001309e-11 -4.34516364859583e-12 -6.61054225868897e-12 -5.75845818545368e-12 -2.32180293529175e-12;0 0 0 0 0;-3.50305843086926e-06 1.76085131953403e-06 8.16661224478572e-07 4.09111042640801e-07 -9.85414469804995e-08;1.44670876127274e-07 -1.41331228923029e-08 -3.06530152369269e-08 -1.46732098927996e-08 -2.30660839364244e-08;-2.00043052422933e-08 1.72145861031776e-09 2.13714615094209e-09 1.02982676689194e-09 -1.64945224692217e-10;1.23552540016991e-09 1.42028470911613e-09 8.79622616627508e-10 -7.44465600265154e-10 -7.17124672589442e-11;-6.67749524914644e-10 -5.77722874934050e-11 3.40077806879472e-11 4.26176076541840e-11 8.23189659748212e-11;-4.62771648935992e-11 -7.24005305716782e-13 1.18233730497485e-12 5.18156973532267e-12 -1.53329687155297e-12;4.75581699468619e-12 -3.79782291469732e-12 1.33077109836853e-12 -1.02426020107120e-12 3.10385019249130e-13;1.66486090578792e-12 1.08573672403649e-12 1.26268044166279e-13 -1.23509297742757e-13 -1.81842007284038e-13;0 0 0 0 0;9.93870680202303e-08 -1.85264736035628e-06 -5.58942734710854e-07 -5.54183448316270e-07 -3.95581289689398e-08;7.88329069002365e-08 2.04810091451078e-08 3.74588851000076e-09 3.42429296613803e-08 -2.00840228416712e-08;-5.93700447329696e-10 -6.57499436973459e-10 -6.90560448220751e-09 3.56586371051089e-09 7.33310245621566e-11;-6.38101662363634e-11 4.23668020216529e-10 -2.43764895979202e-10 -9.31466610703172e-11 -3.17491457845975e-10;1.50943725382470e-11 -6.11641188685078e-11 -4.37018785685645e-11 -2.32871158949602e-11 4.19757251950526e-11;-1.18165328825853e-11 -9.91299557532438e-13 6.40908678055865e-14 2.41049422936434e-12 -8.20746054454953e-14;6.01892101914838e-12 -8.78487122873450e-13 -1.58887481332294e-12 -3.13556902469604e-13 5.14523727801645e-14;-1.50791729401891e-13 -1.45234807159695e-13 1.65302377570887e-13 -5.77094211651483e-15 9.22218953528393e-14;-1.85618902787381e-14 5.64333811864051e-14 -9.94311377945570e-15 -2.40992156199999e-15 -2.19196760659665e-14;0 0 0 0 0;-8.16252352075899e-08 1.61725487723444e-06 9.55522506715921e-07 4.02436267433511e-07 -2.80682052597712e-07;7.68684790328630e-09 -5.00940723761353e-09 -2.43640127974386e-08 -2.59119930503129e-08 3.35015169182094e-08;7.97903115186673e-09 3.73803883416618e-09 3.27888334636662e-09 1.37481300578804e-09 -1.10677168734482e-10;-1.67853012769912e-09 -1.61405252173139e-10 -1.98841576520056e-10 -1.46591506832192e-11 9.35710487804660e-11;4.08807084343221e-11 -3.74514169689568e-11 -3.03638493323910e-11 -5.02332555734577e-12 -8.03417498408344e-12;6.48922619024579e-12 1.96166891023817e-12 -1.96968755122868e-12 -5.20970156382361e-12 -1.62656885103402e-12;1.28603518902875e-12 -4.88146958435109e-13 -3.37034886991840e-13 1.37393696103000e-14 4.41398325716943e-14;1.48670014793021e-13 4.41636026364555e-14 2.06210477976005e-14 -3.43717583585390e-14 -1.21693704024213e-14;-1.67624180330244e-14 6.59317111144238e-15 2.57238525440646e-15 -3.21568425020512e-17 5.29659568026553e-15;7.85453466393227e-16 6.91252183915939e-16 -1.20540764178454e-15 -3.85803892583301e-16 3.46606994632006e-16;0 0 0 0 0;2.86710087625579e-06 -1.68179842305865e-06 -8.48306772016870e-07 -7.08798062479598e-07 -1.27469453733635e-07;2.11824305734993e-09 2.02274279084379e-08 1.61862253091554e-08 3.25597167111807e-08 3.40868964045822e-09;1.21757111431438e-08 1.68405530472906e-09 1.55379338018638e-09 -3.81467795805531e-10 2.53316405545058e-09;-9.98413758659768e-11 5.38382145421318e-10 3.92629628330704e-10 -1.43067134097778e-10 3.74959329667113e-12;-1.57270407028909e-11 -9.02797202317592e-12 8.45997059887690e-12 4.71474382524218e-12 5.41880986596427e-12;-1.20658618702054e-12 7.12940685593433e-13 1.02148613026937e-12 1.63063852348169e-13 1.74048793197708e-13;3.80559390991789e-13 1.19678271353485e-13 9.72859455604188e-14 5.42642400031729e-14 8.18796710714586e-14;-4.69629218656902e-14 5.59889038686206e-15 2.05363292795059e-15 5.38599403288686e-15 -2.68929559474202e-15;-1.88759348081742e-14 5.20975954705924e-15 -4.43585653096395e-16 5.57436617793556e-16 -3.95922805817677e-16;-9.80871456373282e-16 2.50857658461759e-17 -1.24253000050963e-16 6.00857065211394e-17 3.53799635311500e-18;2.49370713054872e-16 -1.49119714269816e-17 -3.12276052640583e-17 -2.42001662334001e-17 -1.69766504318143e-17;0 0 0 0 0;-1.69222102455713e-06 1.64277906173064e-06 5.28855114364096e-07 4.28159853268650e-07 -1.57362445882665e-07;1.67656782413678e-08 -3.77746114074055e-08 -2.21564555842165e-08 -3.37071806992217e-08 1.47454008739800e-08;1.06080499491408e-08 3.21990403709678e-09 3.87301757435359e-09 2.92241827834347e-10 -1.86619473655742e-11;1.62399669665839e-10 3.51322865845172e-10 2.67086377702958e-11 -1.31596563625491e-10 3.14164569507034e-11;-2.02180016657259e-11 2.03305178342732e-11 6.34969032565839e-12 5.99522296668787e-12 -4.46275273451008e-12;-9.88409290158885e-13 -1.47692750858224e-13 3.14655550730530e-13 -2.41857189187879e-13 4.47727504501486e-13;1.71430777754854e-13 1.73950835042486e-13 5.92323956541558e-14 8.06625710171825e-15 2.33252485755634e-14;-1.74184545690134e-15 -8.18003353124179e-16 -6.62369006497819e-16 4.16303374396147e-15 7.06513748014024e-15;-6.02936238677014e-15 1.89241084885229e-15 1.99097881944270e-17 -6.99974290696640e-16 -2.69504942597709e-17;-4.65632962602379e-16 3.70281995445114e-18 -9.04232973763345e-17 2.20847370761932e-17 7.62909453726566e-17;-6.25921477907943e-17 -2.10532795609842e-17 -1.03808073867183e-17 1.15091380049019e-18 4.66794445408388e-19;9.39427013576903e-18 9.17044662931859e-19 2.04132745117549e-18 -1.72364063154625e-19 -1.18098896532163e-18];
bnm_ch = [0 0 0 0 0;0 0 0 0 0;3.44092035729033e-05 -1.21876825440561e-05 -1.87490665238967e-05 -2.60980336247863e-05 4.31639313264615e-06;0 0 0 0 0;-2.60125613000133e-05 1.70570295762269e-05 3.08331896996832e-05 1.66256596588688e-05 -1.07841055501996e-05;8.74011641844073e-06 -2.25874169896607e-06 6.50985196673747e-07 1.30424765493752e-06 -1.85081244549542e-07;0 0 0 0 0;3.77496505484964e-05 -1.08198973553337e-05 -1.67717574544937e-05 -3.22476096673598e-05 1.12281888201134e-05;-7.68623378647958e-07 -4.01400837153063e-06 -2.16390246700835e-06 -1.76912959937924e-06 -1.12740084951955e-06;-2.37092815818895e-06 -9.52317223759653e-07 -2.22722065579131e-07 -6.25157619772530e-08 1.86582003894639e-08;0 0 0 0 0;-6.10254317785872e-05 -2.51815503068494e-05 2.01046207874667e-05 7.21107723367308e-06 -1.30692058660457e-05;-9.60655417241537e-06 -7.31381721742373e-06 -2.52767927589636e-06 9.09039973214621e-07 -6.76454911344246e-07;-2.25743206384908e-08 2.33058746737575e-07 2.24746779293445e-07 6.78551351968876e-08 1.25076011387284e-07;-2.25744112770133e-07 -1.44429560891636e-07 -2.96810417448652e-08 -5.93858519742856e-08 -2.43210229455420e-08;0 0 0 0 0;7.45721015256308e-06 -3.81396821676410e-05 -1.41086198468687e-05 -2.28514517574713e-05 7.28638705683277e-06;-5.77517778169692e-06 -3.93061211403839e-06 -2.17369763310752e-06 -1.48060935583664e-07 -2.74200485662814e-07;4.52962035878238e-07 9.80990375495214e-07 4.67492045269286e-07 -8.31032252212116e-09 1.69426023427740e-07;7.20536791795515e-10 2.75612253452141e-09 2.47772119382536e-09 4.30621825021233e-09 -2.86498479499428e-08;-2.46253956492716e-08 -3.10300833499669e-09 8.06559148724445e-09 2.98197408430123e-10 6.32503656532846e-09;0 0 0 0 0;-6.01147094179306e-05 -3.16631758509869e-05 4.10038115100010e-06 3.55215057231403e-07 -2.23606515237408e-06;-2.85937516921923e-06 -3.67775706610630e-06 -5.06445540401637e-07 8.21776759711184e-07 -5.98690271725558e-07;7.77122595418965e-07 3.60896376754085e-07 3.88610487893381e-07 -4.39533892679537e-08 -6.26882227849174e-08;1.05759993661891e-07 2.58009912408833e-08 -1.51356049060972e-08 -1.13335813107412e-09 5.37470857850370e-10;7.99831506181984e-09 1.67423735327465e-09 2.94736760548677e-09 -1.56727133704788e-09 8.46186800849124e-10;3.07727104043851e-09 3.93584215798484e-10 3.86721562770643e-11 1.72181091277391e-10 -2.16915737920145e-10;0 0 0 0 0;-1.16335389078126e-05 -1.39864676661484e-05 2.52546278407717e-06 -8.79152625440188e-06 -8.97665132187974e-06;-3.95874550504316e-06 -1.17976262528730e-07 7.03189926369300e-07 3.38907065351535e-07 -3.67714052493558e-07;2.29082449370440e-07 5.72961531093329e-07 4.21969662578894e-08 1.24112958141431e-08 9.56404486571888e-08;1.44631865298671e-09 6.19368473895584e-09 1.67110424041236e-09 2.57979463602951e-09 -6.90806907510366e-09;1.77235802019153e-09 -8.14388846228970e-10 4.50421956523579e-09 5.67452314909707e-10 2.47610443675560e-09;4.85932343880617e-10 2.24864117422804e-10 -2.22534534468511e-10 -7.96395824973477e-11 3.12587399902493e-12;-3.20173937255409e-11 -1.29872402028088e-11 -4.24092901203818e-11 2.66570185704416e-11 -5.25164954403909e-12;0 0 0 0 0;-1.36010179191872e-05 1.77873053642413e-05 4.80988546657119e-06 3.46859608161212e-06 -1.73247520896541e-06;2.00020483116258e-06 2.43393064079673e-06 1.21478843695862e-06 1.95582820041644e-07 -3.11847995109088e-07;-8.13287218979310e-09 1.05206830238665e-08 6.54040136224164e-09 -1.96402660575990e-08 -1.40379796070732e-08;4.01291020310740e-08 2.92634301047947e-08 6.04179709273169e-09 8.61849065020545e-10 5.98065429697245e-09;-1.06149335032911e-09 -4.39748495862323e-10 8.83040310269353e-10 3.49392227277679e-10 8.57722299002622e-10;-1.25049888909390e-11 2.05203288281631e-10 1.37817670505319e-11 6.82057794430145e-11 -9.41515631694254e-11;7.47196022644130e-12 -2.51369898528782e-11 -2.12196687809200e-11 1.55282119505201e-11 9.99224438231805e-12;-7.90534019004874e-13 3.55824506982589e-12 8.00835777767281e-13 8.73460019069655e-13 1.34176126600106e-12;0 0 0 0 0;3.12855262465316e-05 1.31629386003608e-05 2.65598119437581e-06 8.68923340949135e-06 -7.51164082949678e-06;1.56870792650533e-06 1.89227301685370e-06 4.15620385341985e-07 -2.74253787880603e-07 -4.28826210119200e-07;-9.99176994565587e-08 -1.10785129426286e-07 -1.10318125091182e-07 6.22726507350764e-09 -3.39214566386250e-08;1.24872975018433e-08 1.10663206077249e-08 5.40658975901469e-09 -2.79119137105115e-09 -2.47500096192502e-09;1.11518917154060e-10 -4.21965763244849e-10 3.26786005211229e-10 1.93488254914545e-10 7.00774679999972e-10;1.50889220040757e-10 1.03130002661366e-10 -3.09481760816903e-11 -4.47656630703759e-11 -7.36245021803800e-12;-1.91144562110285e-12 -1.11355583995978e-11 -1.76207323352556e-11 8.15289793192265e-12 3.45078925412654e-12;-2.73248710476019e-12 -1.65089342283056e-13 -2.20125355220819e-13 5.32589191504356e-13 5.70008982140874e-13;8.06636928368811e-13 1.30893069976672e-13 9.72079137767479e-14 3.87410156264322e-14 -5.56410013263563e-14;0 0 0 0 0;2.02454485403216e-05 -9.77720471118669e-06 -4.35467548126223e-06 2.19599868869063e-06 -3.26670819043690e-06;-3.21839256310540e-08 8.38760368015005e-07 -5.08058835724060e-07 4.16177282491396e-08 1.53842592762120e-07;-1.57377633165313e-07 -7.86803586842404e-08 -7.40444711426898e-08 3.15259864117954e-08 5.60536231567172e-09;-3.26080428920229e-10 -3.14576780695439e-09 8.46796096612981e-10 -2.59329379174262e-09 -8.01054756588382e-10;-4.58725236153576e-11 -6.87847958546571e-11 8.18226480126754e-12 1.81082075625897e-10 1.74510532938256e-10;7.60233505328792e-11 4.76463939581321e-11 -2.47198455442033e-11 -8.83439688929965e-12 5.93967446277316e-13;-8.92919292558887e-12 -4.38524572312029e-12 -4.02709146060896e-12 4.84344426425295e-12 5.12869042781520e-12;1.91518361809952e-12 3.06846255371817e-13 -2.44830265306345e-13 7.86297493099244e-14 2.72347805801980e-13;9.09936624159538e-14 7.20650818861447e-15 2.45383991578283e-14 -4.79580974186462e-15 3.64604724046944e-14;-4.63611142770709e-14 1.73908246420636e-15 -4.41651410674801e-15 -6.61409045306922e-16 -1.60016049099639e-15;0 0 0 0 0;6.17105245892845e-06 -1.04342983738457e-05 -1.72711741097994e-05 -8.16815967888426e-07 3.42789959967593e-06;-2.44014060833825e-07 2.06991837444652e-07 -3.85805819475679e-07 1.67162359832166e-08 4.15139610402483e-07;8.18199006804020e-08 -3.20013409049159e-08 5.94000906771151e-08 2.24122167188946e-08 -1.33796186160409e-08;7.66269294674338e-11 -6.07862178874828e-10 4.95795757186248e-10 -3.07589245481422e-10 3.44456287710689e-10;-1.84076250254929e-10 -1.30985312312781e-10 -1.52547325533276e-10 -2.51000125929512e-11 -1.93924012590455e-11;-2.93307452197665e-11 2.88627386757582e-11 5.58812021182217e-12 -1.68692874069187e-13 1.80464313900575e-12;-9.59053874473003e-13 6.04803122874761e-13 -9.80015608958536e-13 1.70530372034214e-12 1.70458664160775e-12;2.80169588226043e-13 9.09573148053551e-14 2.16449186617004e-14 1.15550091496353e-13 4.97772796761321e-14;-3.04524400761371e-14 3.42845631349694e-14 2.44230630602064e-14 5.76017546103056e-16 -9.74409465961093e-15;5.98765340844291e-15 -2.63942474859535e-15 -1.80204805804437e-15 -1.84981819321183e-16 -5.85073392163660e-16;-2.37069441910133e-15 2.87429226086856e-16 -1.67055963193389e-16 2.72110684914090e-18 8.46646962667892e-17;0 0 0 0 0;-2.71386164105722e-05 -1.41834938338454e-05 -2.00777928859929e-07 5.94329804681196e-07 8.61856994375586e-06;-3.93656495458664e-08 -6.36432821807576e-07 -2.47887475106438e-07 -2.64906446204966e-08 1.10689794197004e-07;5.25319489188562e-08 9.00866357158695e-09 5.00693379572512e-08 2.47269011056404e-08 -7.27648556194598e-09;1.87207107149043e-09 -1.46428282396138e-09 -2.71812237167257e-10 8.44902265891466e-10 -5.62683870906027e-10;-1.08295119666184e-10 4.75553388543793e-11 -5.49429386495686e-11 -6.60907871731611e-11 -5.97347322824822e-11;-4.95118306815571e-12 5.31083735234970e-13 -1.93679746327378e-12 -1.61770521840510e-12 1.23276727202510e-11;6.68582682909900e-13 7.38288575160449e-13 5.47630483499201e-13 -1.00770258118914e-13 -1.65564928475981e-13;5.80963409268471e-14 6.93474288078737e-14 6.60728092794315e-15 -5.21029056725202e-15 -1.11283532854883e-16;-4.10567742688903e-15 1.62252646805882e-14 1.00774699865989e-14 -2.44793214897877e-16 -1.59283906414563e-15;1.84669506619904e-17 8.28473337813919e-17 -1.53400662078899e-16 -5.01060672199689e-17 -2.20727935766132e-16;2.65355116203636e-16 -3.70233146147684e-17 3.52689394451586e-18 -8.62215942516328e-18 9.26909361974526e-18;9.94266950643135e-17 4.17028699663441e-18 -7.65153491125819e-21 -5.62131270981041e-18 -3.03732817297438e-18];
bnm_cw = [0 0 0 0 0;0 0 0 0 0;-0.000209104872912563 -1.41530274973540e-05 3.00318745764815e-05 -1.82864291318284e-05 -7.62965409959238e-06;0 0 0 0 0;-0.000186336519900275 0.000191256553935638 7.28356195304996e-05 3.59637869639906e-05 -2.53927226167388e-05;0.000108195343799485 -6.97050977217619e-05 -6.68037133871099e-05 2.30387653190503e-05 -1.22735483925784e-05;0 0 0 0 0;0.000119941091277039 -7.70547844186875e-05 -8.15376297964528e-05 1.06005789545203e-05 2.31177232268720e-05;-1.77494760217164e-05 -1.37061385686605e-05 -1.74805936475816e-05 -6.91745900867532e-07 -7.10231790947787e-06;-1.47564103733219e-05 2.08890785485260e-06 3.19876879447867e-06 9.43984664503715e-07 -4.90480527577521e-06;0 0 0 0 0;4.93300138389457e-05 -6.77641298460617e-05 -3.25043347246397e-05 8.33226714911921e-06 8.11499972792905e-06;-2.80449863471272e-05 -1.04367606414606e-05 1.64473584641163e-07 -3.57420965807816e-06 2.95887156564038e-06;1.88835280111533e-06 5.69125761193702e-07 -2.22757382799409e-06 -1.96699131032252e-07 -2.91861219283659e-07;-4.69918971436680e-06 -7.00778948636735e-07 2.97544157334673e-09 3.86100512544410e-07 2.30939653701027e-07;0 0 0 0 0;1.77050610394149e-05 -3.18353071311574e-05 3.04232260950316e-05 -6.26821316488169e-05 -1.75094810002378e-06;9.25605901565775e-06 -8.25179123302247e-06 6.74032752408358e-06 3.22192289084524e-06 6.09414500075259e-06;4.28233825242200e-06 2.10470570087927e-07 -4.75050074985668e-07 -4.89382663470592e-07 8.75232347469207e-07;8.50393520366934e-07 1.58764911467186e-07 -2.16267638321210e-07 -7.43341300487416e-10 1.75131729813230e-07;-2.87064111623119e-07 4.50393893102830e-08 6.63315044416690e-08 7.61199387418853e-08 -6.05694385243652e-09;0 0 0 0 0;-1.95692079507947e-05 5.15486098887851e-05 3.00852761598173e-05 1.21485028343416e-05 -6.72450521493428e-06;5.34496867088158e-06 3.90973451680699e-06 3.70148924718425e-06 5.73731499938212e-08 5.52258220288780e-07;3.39950838185315e-07 -5.63443976772634e-07 4.52082211980595e-07 -2.57094645806243e-07 -6.84885762924729e-08;2.15793276880684e-07 2.05911354090873e-07 1.33747872341142e-08 -2.07997626478952e-08 -3.69812938736019e-08;2.11952749403224e-09 4.04317822544732e-08 2.40972024883650e-09 8.56289126938059e-09 2.31035283490200e-08;-2.08402298813248e-09 -8.50243600879112e-09 2.60895410117768e-09 -6.69156841738591e-10 -5.16280278087006e-09;0 0 0 0 0;0.000124901291436683 -5.70770326719086e-05 -8.44887248105015e-05 -3.11442665354698e-05 -1.12982893252046e-05;-8.38934444233944e-06 1.56860091415414e-06 -1.77704563531825e-06 -5.70219068898717e-08 -4.30377735031244e-06;3.72965318017681e-07 6.98175439446187e-07 1.75760544807919e-08 1.59731284857151e-07 3.62363848767891e-07;-2.32148850787091e-07 -4.21888751852973e-08 8.35926113952108e-08 -2.24572480575674e-08 -6.92114100904503e-08;-2.92635642210745e-09 3.38086229163415e-09 4.72186694662901e-09 -8.32354437305758e-11 4.19673890995627e-09;-1.26452887692900e-09 1.91309690886864e-09 1.54755631983655e-09 -1.09865169400249e-09 1.83645326319994e-10;9.92539437011905e-10 -2.96318203488300e-10 1.17466020823486e-10 -5.00185957995526e-10 -8.54777591408537e-11;0 0 0 0 0;-0.000182885335404854 7.27424724520089e-05 3.05286278023427e-05 2.55324463432562e-05 -6.39859510763234e-06;-5.21449265232557e-06 -6.70572386081398e-06 -3.95473351292738e-06 -6.41023334372861e-07 -3.11616331059009e-06;2.37090789071727e-07 3.58427517014705e-07 2.55709192777007e-07 8.44593804408541e-08 9.27243162355359e-09;7.24370898432057e-08 -7.43945120337710e-09 8.61751911975683e-10 -2.34651212610623e-08 2.94052921681456e-09;-1.22127317934425e-08 -3.89758984276768e-09 4.12890383904924e-11 2.06528068002723e-09 1.73488696972270e-09;-5.44137406907620e-10 -4.81034553189921e-10 -2.56101759039694e-11 3.21880564410154e-10 -2.70195343165250e-11;1.08394225300546e-10 -7.99525492688661e-11 1.73850287030654e-10 -8.06390014426271e-11 -7.63143364291160e-13;-3.41446959267441e-11 2.72675729042792e-11 5.69674704865345e-12 -3.38402998344892e-12 -2.96732381931007e-12;0 0 0 0 0;2.91161315987250e-05 -7.24641166590735e-05 -8.58323519857884e-06 -1.14037444255820e-05 1.32244819451517e-05;1.24266748259826e-06 -4.13127038469802e-06 -8.47496394492885e-07 5.48722958754267e-07 -1.98288551821205e-06;-1.70671245196917e-08 1.36891127083540e-08 -2.80901972249870e-07 -5.45369793946222e-09 -9.58796303763498e-08;1.14115335901746e-08 2.79308166429178e-08 -1.71144803132413e-08 4.86116243565380e-09 -8.13061459952280e-09;-1.19144311035824e-09 -1.28197815211763e-09 -1.22313592972373e-09 6.23116336753674e-10 2.11527825898689e-09;4.94618645030426e-10 -1.01554483531252e-10 -3.58808808952276e-10 1.23499783028794e-10 -1.21017599361833e-10;1.33959569836451e-10 -1.87140898812283e-11 -3.04265350158941e-11 -1.42907553051431e-11 -1.09873858099638e-11;1.30277419203512e-11 -4.95312627777245e-12 2.23070215544358e-12 1.66450226016423e-12 6.26222944728474e-12;-4.40721204874728e-12 2.99575133064885e-12 -1.54917262009097e-12 8.90015664527060e-14 -1.59135267012937e-12;0 0 0 0 0;-4.17667211323160e-05 1.39005215116294e-05 1.46521361817829e-05 3.23485458024416e-05 -8.57936261085263e-06;9.48491026524450e-07 1.67749735481991e-06 6.80159475477603e-07 -1.34558044496631e-06 1.62108231492249e-06;-2.67545753355631e-07 -3.31848493018159e-08 1.05837219557465e-07 1.55587655479400e-07 -2.84996014386667e-08;-5.15113778734878e-08 8.83630725241303e-09 3.36579455982772e-09 -6.22350102096402e-09 5.03959133095369e-09;2.04635880823035e-11 -1.07923589059151e-09 -6.96482137669712e-10 -4.70238500452793e-10 -6.60277903598297e-10;-2.41897168749189e-11 1.33547763615216e-10 -5.13534673658908e-11 -8.32767177662817e-11 5.72614717082428e-11;7.55170562359940e-12 -1.57123461699055e-11 -1.48874069619124e-11 -7.10529462981252e-13 -7.99006335025107e-12;2.41883156738960e-12 2.97346980183361e-12 1.28719977731450e-12 -2.49240876894143e-12 6.71155595793198e-13;4.16995565336914e-13 -1.71584521275288e-13 -7.23064067359978e-14 2.45405880599037e-13 4.43532934905830e-13;3.56937508828997e-14 2.43012511260300e-14 -7.96090778289326e-14 -1.59548529636358e-14 8.99103763000507e-15;0 0 0 0 0;0.000117579258399489 -4.52648448635772e-05 -2.69130037097862e-05 -3.82266335794366e-05 -4.36549257701084e-06;-1.43270371215502e-06 1.21565440183855e-06 8.53701136074284e-07 1.52709810023665e-06 1.22382663462904e-06;3.06089147519664e-07 9.79084123751975e-08 7.96524661441178e-08 4.54770947973458e-08 2.22842369458882e-07;-9.94254707745127e-09 1.43251376378012e-08 1.93911753685160e-08 -6.52214645690987e-09 -1.97114016452408e-09;-9.20751919828404e-10 -9.44312829629076e-10 7.24196738163952e-11 -6.71801072324561e-11 2.33146774065873e-10;-1.43544298956410e-11 1.78464235318769e-10 7.69950023012326e-11 -4.22390057304453e-12 3.05176324574816e-11;-7.88053753973990e-12 -3.20207793051003e-12 1.01527407317625e-12 6.02788185858449e-12 1.14919530900453e-11;-1.21558899266069e-12 5.31300597882986e-13 3.44023865079264e-13 -6.22598216726224e-14 -5.47031650765402e-14;-4.15627948750943e-13 2.77620907292721e-13 -8.99784134364011e-14 1.07254247320864e-13 6.85990080564196e-14;-3.91837863922901e-14 9.74714976816180e-15 6.79982450963903e-15 -2.41420876658572e-15 -2.20889384455344e-15;9.25912068402776e-15 -4.02621719248224e-15 -2.43952036351187e-15 -1.97006876049866e-15 1.03065621527869e-16;0 0 0 0 0;-0.000103762036940193 4.38145356960292e-05 2.43406920349913e-05 7.89103527673736e-06 -1.66841465339160e-05;-1.18428449371744e-06 -1.30188721737259e-06 -1.88013557116650e-06 -1.01342046295303e-06 9.21813037802502e-07;1.51836068712460e-07 1.11362553803933e-07 1.55375052233052e-07 1.94450910788747e-09 -1.73093755828342e-08;-3.77758211813121e-09 1.23323969583610e-08 1.72510045250302e-09 -1.88609789458597e-09 1.28937597985937e-09;-1.07947760393523e-09 5.26051570105365e-10 -3.67657536332496e-11 3.16110123523840e-10 -3.24273198242170e-10;-2.00385649209820e-12 2.54703869682390e-11 4.08563622440851e-12 -4.83350348928636e-11 -3.98153443845079e-13;2.73094467727215e-12 5.08900664114903e-12 -7.66669089075134e-13 2.50015592643012e-12 4.29763262853853e-12;6.53946487537890e-13 -2.24958413781008e-13 6.74638861781238e-15 3.28537647613903e-14 2.54199700290116e-13;-1.09122051193505e-13 8.36362392931501e-14 -3.90750153912300e-14 -5.44915910741950e-14 2.43816947219217e-14;-1.41882561550134e-14 1.00455397812713e-14 2.63347255121581e-15 1.53043256823601e-15 2.49081021428095e-15;-1.17256193152654e-15 1.05648985031971e-16 1.31778372453016e-16 1.44815198666577e-16 -3.72532768618480e-16;2.66203457773766e-16 -7.67224608659658e-17 3.51487351031864e-18 4.10287131339291e-17 -6.72171711728514e-17];


%READ THE RESPECTIVE LINES FROM THE MATRICES
anm_bh_A0 = anm_bh(:,1);   anm_bh_A1 = anm_bh(:,2);   anm_bh_B1 = anm_bh(:,3);   anm_bh_A2 = anm_bh(:,4);   anm_bh_B2 = anm_bh(:,5);
anm_bw_A0 = anm_bw(:,1);   anm_bw_A1 = anm_bw(:,2);   anm_bw_B1 = anm_bw(:,3);   anm_bw_A2 = anm_bw(:,4);   anm_bw_B2 = anm_bw(:,5);
anm_ch_A0 = anm_ch(:,1);   anm_ch_A1 = anm_ch(:,2);   anm_ch_B1 = anm_ch(:,3);   anm_ch_A2 = anm_ch(:,4);   anm_ch_B2 = anm_ch(:,5);
anm_cw_A0 = anm_cw(:,1);   anm_cw_A1 = anm_cw(:,2);   anm_cw_B1 = anm_cw(:,3);   anm_cw_A2 = anm_cw(:,4);   anm_cw_B2 = anm_cw(:,5);
bnm_bh_A0 = bnm_bh(:,1);   bnm_bh_A1 = bnm_bh(:,2);   bnm_bh_B1 = bnm_bh(:,3);   bnm_bh_A2 = bnm_bh(:,4);   bnm_bh_B2 = bnm_bh(:,5);
bnm_bw_A0 = bnm_bw(:,1);   bnm_bw_A1 = bnm_bw(:,2);   bnm_bw_B1 = bnm_bw(:,3);   bnm_bw_A2 = bnm_bw(:,4);   bnm_bw_B2 = bnm_bw(:,5);
bnm_ch_A0 = bnm_ch(:,1);   bnm_ch_A1 = bnm_ch(:,2);   bnm_ch_B1 = bnm_ch(:,3);   bnm_ch_A2 = bnm_ch(:,4);   bnm_ch_B2 = bnm_ch(:,5);
bnm_cw_A0 = bnm_cw(:,1);   bnm_cw_A1 = bnm_cw(:,2);   bnm_cw_B1 = bnm_cw(:,3);   bnm_cw_A2 = bnm_cw(:,4);   bnm_cw_B2 = bnm_cw(:,5);


%ASSIGNMENT
el = satELEV;
polDist = pi/2 - lat;
    
%a.) CALCULATE LEGENDRE POLYNOMIALS
 
%DEGREE n AND ORDER m
nmax = 12;

%UNIT VECTOR
x = sin(polDist)*cos(lon);
y = sin(polDist)*sin(lon);
z = cos(polDist);

%LEGENDRE POLYNOMIALS
V(1,1) = 1;
W(1,1) = 0;
V(2,1) = z * V(1,1);
W(2,1) = 0;

for n = 2:nmax
    V(n+1,1) = ((2*n-1) * z * V(n,1) - (n-1) * V(n-1,1)) / n;
    W(n+1,1) = 0;
end

for m = 1:nmax
    V(m+1,m+1) = (2*m-1) * (x*V(m,m) - y*W(m,m));
    W(m+1,m+1) = (2*m-1) * (x*W(m,m) + y*V(m,m));
    if (m < nmax) 
        V(m+2,m+1) = (2*m+1) * z * V(m+1,m+1);
        W(m+2,m+1) = (2*m+1) * z * W(m+1,m+1);
    end
    for n = m+2:nmax
        V(n+1,m+1) = ((2*n-1)*z*V(n,m+1) - (n+m-1)*V(n-1,m+1)) / (n-m);
        W(n+1,m+1) = ((2*n-1)*z*W(n,m+1) - (n+m-1)*W(n-1,m+1)) / (n-m);
    end
end

%b.)DETERMINE THE COEFFICIENTS bh, bw, ch and cw 
%INITIALIZE
bh_A0 = 0; bh_A1 = 0; bh_B1 = 0; bh_A2 = 0; bh_B2 = 0;
bw_A0 = 0; bw_A1 = 0; bw_B1 = 0; bw_A2 = 0; bw_B2 = 0;
ch_A0 = 0; ch_A1 = 0; ch_B1 = 0; ch_A2 = 0; ch_B2 = 0;
cw_A0 = 0; cw_A1 = 0; cw_B1 = 0; cw_A2 = 0; cw_B2 = 0;
i = 0;

for n = 0:nmax
    for m = 0:n
        
        i = i+1;
        
        bh_A0 = bh_A0 + (anm_bh_A0(i)*V(n+1,m+1) + bnm_bh_A0(i)*W(n+1,m+1));
        bh_A1 = bh_A1 + (anm_bh_A1(i)*V(n+1,m+1) + bnm_bh_A1(i)*W(n+1,m+1));
        bh_B1 = bh_B1 + (anm_bh_B1(i)*V(n+1,m+1) + bnm_bh_B1(i)*W(n+1,m+1));
        bh_A2 = bh_A2 + (anm_bh_A2(i)*V(n+1,m+1) + bnm_bh_A2(i)*W(n+1,m+1));
        bh_B2 = bh_B2 + (anm_bh_B2(i)*V(n+1,m+1) + bnm_bh_B2(i)*W(n+1,m+1));
        
        bw_A0 = bw_A0 + (anm_bw_A0(i)*V(n+1,m+1) + bnm_bw_A0(i)*W(n+1,m+1));
        bw_A1 = bw_A1 + (anm_bw_A1(i)*V(n+1,m+1) + bnm_bw_A1(i)*W(n+1,m+1));
        bw_B1 = bw_B1 + (anm_bw_B1(i)*V(n+1,m+1) + bnm_bw_B1(i)*W(n+1,m+1));
        bw_A2 = bw_A2 + (anm_bw_A2(i)*V(n+1,m+1) + bnm_bw_A2(i)*W(n+1,m+1));
        bw_B2 = bw_B2 + (anm_bw_B2(i)*V(n+1,m+1) + bnm_bw_B2(i)*W(n+1,m+1));
        
        ch_A0 = ch_A0 + (anm_ch_A0(i)*V(n+1,m+1) + bnm_ch_A0(i)*W(n+1,m+1));
        ch_A1 = ch_A1 + (anm_ch_A1(i)*V(n+1,m+1) + bnm_ch_A1(i)*W(n+1,m+1));
        ch_B1 = ch_B1 + (anm_ch_B1(i)*V(n+1,m+1) + bnm_ch_B1(i)*W(n+1,m+1));
        ch_A2 = ch_A2 + (anm_ch_A2(i)*V(n+1,m+1) + bnm_ch_A2(i)*W(n+1,m+1));
        ch_B2 = ch_B2 + (anm_ch_B2(i)*V(n+1,m+1) + bnm_ch_B2(i)*W(n+1,m+1));
        
        cw_A0 = cw_A0 + (anm_cw_A0(i)*V(n+1,m+1) + bnm_cw_A0(i)*W(n+1,m+1));
        cw_A1 = cw_A1 + (anm_cw_A1(i)*V(n+1,m+1) + bnm_cw_A1(i)*W(n+1,m+1));
        cw_B1 = cw_B1 + (anm_cw_B1(i)*V(n+1,m+1) + bnm_cw_B1(i)*W(n+1,m+1));
        cw_A2 = cw_A2 + (anm_cw_A2(i)*V(n+1,m+1) + bnm_cw_A2(i)*W(n+1,m+1));
        cw_B2 = cw_B2 + (anm_cw_B2(i)*V(n+1,m+1) + bnm_cw_B2(i)*W(n+1,m+1));
        
    end
end

%ADDING THE SEASONAL AMPLITUDES FOR THE SPECIFIED DOY TO THE MEAN VALUES
bh = bh_A0 + bh_A1*cos(doy/365.25*2*pi) + bh_B1*sin(doy/365.25*2*pi) + bh_A2*cos(doy/365.25*4*pi) + bh_B2*sin(doy/365.25*4*pi);
bw = bw_A0 + bw_A1*cos(doy/365.25*2*pi) + bw_B1*sin(doy/365.25*2*pi) + bw_A2*cos(doy/365.25*4*pi) + bw_B2*sin(doy/365.25*4*pi);
ch = ch_A0 + ch_A1*cos(doy/365.25*2*pi) + ch_B1*sin(doy/365.25*2*pi) + ch_A2*cos(doy/365.25*4*pi) + ch_B2*sin(doy/365.25*4*pi);
cw = cw_A0 + cw_A1*cos(doy/365.25*2*pi) + cw_B1*sin(doy/365.25*2*pi) + cw_A2*cos(doy/365.25*4*pi) + cw_B2*sin(doy/365.25*4*pi);

%CALCULATING THE HYDROSTATIC AND WET MAPPING FACTORS
vmf3h = (1+(ah/(1+bh/(1+ch))))   /   (sin(el)+(ah/(sin(el)+bh/(sin(el)+ch))));
vmf3w = (1+(aw/(1+bw/(1+cw))))   /   (sin(el)+(aw/(sin(el)+bw/(sin(el)+cw))));

% HEIGHT CORRECTION FOR THE HYDROSTATIC PART [Niell, 1996]
a_ht = 2.53d-5;
b_ht = 5.49d-3;
c_ht = 1.14d-3;
h_ell_km     = h_ell/1000;   % convert height to km
ht_corr_coef = 1/sin(el)   -   (1+(a_ht/(1+b_ht/(1+c_ht))))  /  (sin(el)+(a_ht/(sin(el)+b_ht/(sin(el)+c_ht))));
ht_corr      = ht_corr_coef * h_ell_km;
vmf3h          = vmf3h + ht_corr;

%******************************************END OF vmf3.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%C.3.SUBROUTINE TO COMPUTE VIENNA MAPPING FUNCTION 1(VMF1) COEFFICIENTS
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [ah1,aw1] = VMF1coe(UTCtime,latD,lonD,hgt,VMFgrids,grid_h,...
                             grid_w,grid_resh,grid_resw,Timevar_h,...
                                             Timevar_w,MFh_model,MFw_model)

%**************************************************************************
%***DESCRIPTION:
%              This subroutine determines the Tropospheric Mapping         *
%              Function using the Vienna Mapping Functions 1(VMF1)         * 
%              coefficients(ah,aw) from various VMF1 Models(i.e. either... *
%              from VMF1 grid file(s) or from GPT2w models for Specific site* 
%              given receiver/station position in [lat Lon hgt], receiver  *
%              reception time in the utc time format[Yr Mon Day Hr Min Sec]*
%              ,GPT2w grid values,grid resolution  and time variation      *
%              indicator(0 or 1)                                           *
%--------------------------------------------------------------------------*                                     
%USAGE:                                                                    *                                                           *
%      [ah1,aw1] = VMF1coe(UTCtime,latD,lonD,hgt,VMFgrids,grid_h,...       *
%                             grid_w,grid_resh,grid_resw,Timevar_h,...     *
%                                            Timevar_w,MFh_model,MFw_model)*       
%--------------------------------------------------------------------------*
%****INPUTs:                                                               *
%1.    UTCtime : UTC time in [Year,Month,Day,Hour,Minute,Seconds]          *
%2.       latD : station ellipsoidal latitude in [degrees]                 *
%3.       lonD : station ellipsoidal longitude in [degrees]                *
%4.        hgt : station ellipsoidal height in [meters]                    *
%5.VMFgrids    : Structure array of VMF grid files which can be :          *
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
%--------------------------------------------------------------------------*
%6.      grid   : grid values in cells extracted from the grid file        *
%    grid_res   : grid resolution (°) (possible: 1 or 5)                   *
%7    Timevar   : case 1: no time variation but static quantities          *
%                 case 0: with time variation (annual and semiannual terms)*
%8.  MF_model   :  mapping function model type(i.e.GPT2 OR GPT2w)          *
%--------------------------------------------------------------------------*
%***OUTPUTs:                                                               *
%          ah1  : Vienna Hydrostatic mapping function 1 coefficient        *
%          aw1  : Vienna Wet mapping function  1 coefficient               *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-                                      *
%WRITTEN by:  Samuel Osah,Msc Geomatic Engineering ,2016                   *    
%             Email: osahsamuel@yahoo.ca                                   *
%             Tel:+233 (0)246137410/+233 (0)509438484                      *        
%==========================================================================+
%**************************************************************************+
%**************************************************************************+
% setappdata(0,'VMF_grid_hw',VMF_grid)
% setappdata(0,'VMF_grid_mf',VMF_grid)
%***CONVERT Lat and Lon IN DEGREES TO RADIAN 
latR = (latD./180).*pi; %[radian]
lonR = (lonD./180).*pi ;%[radian]    

%RETRIEVE SOME SAVED DATA FROM TropModel or TropoModel subroutines
%GET DRY & WET MODELs 
dryModel = getappdata(0,'dryModel');%get selected Hydrostatic model
wetModel = getappdata(0,'wetModel');%get selected Wet model

%GET VMF grid FILE FOUND/NOT FOUND INDICATOR
VMF_grid_found = getappdata(0,'VMF_grid_found');

%*****************1ST PREAMBLE
%SEARCH VMF1 GRID FILEs IN DIRECTORY & COMPUTE MF COEFFICIENTs
%USING "SearchReadVMFgrid.m" FUNCTION
%NOTE:
% "SearchReadVMFgrid.m" will search and read the following grid file type : 
%1.VMFG_yrmonthday.H00 or VMFG_yrmonthday.H06 or VMFG_yrmonthday.H12 or
%  VMFG_yrmonthday.H18 eg:[VMFG_20180101.H00, VMFG_20180101.H06,... 
%                          VMFG_20180101. H12,VMFG_20180101.H18]

%NOTE: 
%     GENERALLY, VMF1 COEFFICIENTS ARE EXTRACTED OR COMPUTED FROM VMF1 GRID FILES.
%     WHERE VMF1 GRID FILES ARE NOT AVAILABLE, WE THEN RESORT TO THOSE FROM THE GPT2w or GPT2 MODELS

%IF HYDROSTATIC & WET MAPPING FUNCTION(MF) MODELS ARE THE SAME
if any([strncmpi(MFh_model,MFw_model,12),strfind(MFh_model,MFw_model),any([strncmpi(MFh_model,'VMF gridded ZHD',15),...
        strncmpi(MFw_model,'VMF gridded ZWD',15)])])
    
    %IF MF MODEL IS ANY OF THE VIENNA MAPPING FUNCTION 1(VMF1) MODELS
   if any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),strncmpi(MFh_model,'VMF gridded ZHD',15),...
           strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),strncmpi(MFw_model,'VMF gridded ZWD',15)])
    
       %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines WERE THE VMF GRID MODELS
      if any([any([strncmpi(dryModel,'VMF gridded ZHD',15),strfind(dryModel,'VMF gridded ZHD')])
              any([strncmpi(wetModel,'VMF gridded ZWD',15),strfind(wetModel,'VMF gridded ZWD')])])
          
         %FIRST GET VMF COEFFIFFIENTS FROM TropModel or TropoModel subroutines 
         if all([~isempty(getappdata(0,'ah_vmfg1')),~isempty(getappdata(0,'aw_vmfg1'))])
            
            %RETRIEVE STORED VMF COEFFICIENTS 
            ah1 = getappdata(0,'ah_vmfg1');
            aw1 = getappdata(0,'aw_vmfg1');
            VMF_model = 'VMF1';
            
         else %IF STORED DATA [getappdata(0,'ah_vmfg1')& getappdata(0,'aw_vmfg1')] ARE EMPTY([])
               
              if VMF_grid_found == 1%IF VMF GRID FILES ARE FOUND,THEN COMPUTE VMF COE USING "SearchReadVMFgrid"
                  
                 %Call the "SearchReadVMFgrid.m" function TO COMPUTE VMF1 COEFFICIENTS
                 [ah1,aw1,VMF_model] = SearchReadVMFgrid(UTCtime,latD,lonD,hgt,VMFgrids);
                 
              else %IF VMF GRID FILES ARE NOT FOUND
                  %ASSIGN EMPTY([]) MATRIX
                  ah1       = [];
                  aw1       = [];
                  VMF_model = [];
                  
              end
                  
        
         end
         
      else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines WERE NOT THE VMF GRID MODELS 
          
           if VMF_grid_found == 1 %IF VMF GRID FILES ARE FOUND,THEN COMPUTE VMF COE USING "SearchReadVMFgrid"
               
              %Call the "SearchReadVMFgrid.m" function TO COMPUTE VMF1 COEFFICIENTS
              [ah1,aw1,VMF_model] = SearchReadVMFgrid(UTCtime,latD,lonD,hgt,VMFgrids);
              
           else %IF VMF GRID FILES ARE NOT FOUND
               
                %ASSIGN EMPTY([]) MATRIX
                ah1       = [];
                aw1       = [];
                VMF_model = [];
                  
           end 
          
      end %//if any([any([strncmpi(dryModel,'VMF ZHD',7),strfind(dryModel,'VMF ZHD')])
          %          any([strncmpi(wetModel,'VMF ZWD',7),strfind(wetModel,'VMF ZWD')])])

      %IF 1ST PREAMBLE FAILS, GET COEFFICIENTS FROM GPT2w MODELS 
      %EG:.gpt2_5w.grd or gpt2_1w.grd or gpt2_1wA.grd
      
      %IF VMF GRID MODEL IS THE VMF3 MODEL OR VMF GRID FILES WERE NOT AVAILABLE
      if any([strncmpi(VMF_model,'VMF3',4),all([isempty(VMF_model),isempty(ah1),isempty(aw1)])])

          %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel
          %subroutines ARE THE SAME
         if any([strncmpi(dryModel,wetModel,14),strfind(dryModel,wetModel)]) 
          
            if any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                    strncmpi(dryModel,'GPT2 (5° x 5°)',14)])
             
               %RETRIEVE STORED VMF COEFFICIENTS
               ah1 = getappdata(0,'ah_vmf1');
               aw1 = getappdata(0,'aw_vmf1');
               
            else
                 try 
                       %CHECK LENGTH OF GPT GRID VALUES{FOR GPT2W,LENGTH SHOULD BE 10}
                      if any([length(grid_h) == 10,length(grid_w) == 10])
                
                          %Call the "vmfCOE_gpt2w.m" function TO COMPUTE VMF1 COEFFICIENTS
                          [ah1,aw1] = vmfCOE_gpt2w(UTCtime,latR,lonR,hgt,grid_h,grid_resh,Timevar_h);
           
                      elseif length(grid_h) == 8 %IF LENGTH IS EQUALS 8 {FOR GPT2,LENGTH SHOULD BE 8}
                  
                             %Call the "gvmfCOE_gpt2.m" function TO COMPUTE VMF1 COEFFICIENTS
                             [ah1,aw1] = vmfCOE_gpt2(UTCtime,latR,lonR,hgt,grid_h,Timevar_h);     
                      end     
           
                 catch    
                        %Call the "vmfCOE_gpt2.m" function TO COMPUTE VMF1 COEFFICIENTS
                        [ah1,aw1] = vmfCOE_gpt2(UTCtime,latR,lonR,hgt,grid_h,Timevar_h);
                 end          
               
            end  
            
         else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel %subroutines ARE NOT THE SAME 
          
              %IF DRY MODEL IS ANY OF THE GPT2w MODEL 1° OR 5° GRID & WET MODEL IS NOT
              if all([any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                           strncmpi(dryModel,'GPT2 (5° x 5°)',14)]),...
                      ~any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                            strncmpi(wetModel,'GPT2 (5° x 5°)',14)])])     
                
                 %RETRIEVE STORED VMF COEFFICIENTS
                 if any([strncmpi(getappdata(0,'Drymodel'),'GPT2',4),strfind(getappdata(0,'Drymodel'),'GPT2')])
                    
                    ah1 = getappdata(0,'ah_vmf1');
                    aw1 = getappdata(0,'aw_gpt2_vmf1');
                    
                 elseif any([strncmpi(getappdata(0,'Drymodel'),'GPT2w',5),strfind(getappdata(0,'Drymodel'),'GPT2w')])
                        
                        ah1 = getappdata(0,'ah_vmf1');
                        aw1 = getappdata(0,'aw_gpt2w_vmf1');
                        
                 end
                 
              %IF WET MODEL IS ANY OF THE GPT2w MODEL (1° OR 5° GRID) & DRY MODEL IS NOT   
              elseif all([~any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                               strncmpi(dryModel,'GPT2 (5° x 5°)',14)]),...
                           any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                                strncmpi(wetModel,'GPT2 (5° x 5°)',14)])])     
                
                 %RETRIEVE STORED VMF COEFFICIENTS
                 if any([strncmpi(getappdata(0,'Wetmodel'),'GPT2',4),strfind(getappdata(0,'Wetmodel'),'GPT2')])
                    
                    ah1 = getappdata(0,'ah_gpt2_vmf1');
                    aw1 = getappdata(0,'aw_vmf1');
                    
                 elseif any([strncmpi(getappdata(0,'Wetmodel'),'GPT2w',5),strfind(getappdata(0,'Wetmodel'),'GPT2w')])
                     
                        ah1 = getappdata(0,'ah_gpt2w_vmf1');
                        aw1 = getappdata(0,'aw_vmf1'); 
                        
                 end
                 
              else %IF DRY OR WET TROPO MODEL IS NOT ANY OF THE GPT2w MODEL (1° OR 5° GRID)
                  
                   try 
                       %CHECK LENGTH OF GPT GRID VALUES{FOR GPT2W,LENGTH SHOULD BE 10}
                      if any([length(grid_h) == 10,length(grid_w) == 10])
                
                          %Call the "vmfCOE_gpt2w.m" function TO COMPUTE VMF1 COEFFICIENTS
                          [ah1,aw1] = vmfCOE_gpt2w(UTCtime,latR,lonR,hgt,grid_h,grid_resh,Timevar_h);
           
                      elseif length(grid_h) == 8 %IF LENGTH IS EQUALS 8 {FOR GPT2,LENGTH SHOULD BE 8}
                  
                             %Call the "gvmfCOE_gpt2.m" function TO COMPUTE VMF1 COEFFICIENTS
                             [ah1,aw1] = vmfCOE_gpt2(UTCtime,latR,lonR,hgt,grid_h,Timevar_h);     
                      end    
           
                   catch   
                        %Call the "vmfCOE_gpt2.m" function TO COMPUTE VMF1 COEFFICIENTS
                        [ah1,aw1] = vmfCOE_gpt2(UTCtime,latR,lonR,hgt,grid_h,Timevar_h);
                   end        
                  
              end  %//if all([any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                      %            strncmpi(dryModel,'GPT2 (5° x 5°)',14)]),...
                      %      ~any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                      %            strncmpi(wetModel,'GPT2 (5° x 5°)',14)])])  
              
         end %//if any([strncmpi(drymodel,wetmodel,14),strfind(drymodel,wetmodel)]) 
         
      end  %//if any([strncmpi(VMF_model,'VMF3',4),all([isempty(VMF_model),isempty(ah1),isempty(aw1)])])
 

%IF MF MODEL IS RATHER ANY OF THE GPT2w MODEL(1° OR 5°) AND NOT ANY OF VMF1 MODELS
elseif any([strncmpi(MFh_model,'GPT2w (1° x 1°)',14),strncmpi(MFh_model,'GPT2w (5° x 5°)',14)]),...
       any([strncmpi(MFw_model,'GPT2w (1° x 1°)',14),strncmpi(MFw_model,'GPT2w (5° x 5°)',14)])
       
       %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE THE SAME
       if any([strncmpi(dryModel,wetModel,14),strfind(dryModel,wetModel)]) 
          
            %IF DRY TROPO MODEL IS ANY OF THE GPT2w MODEL(1° OR 5°)
            if any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14)])
             
               %RETRIEVE STORED VMF COEFFICIENTS 
               ah1 = getappdata(0,'ah_vmf1');
               aw1 = getappdata(0,'aw_vmf1');  
            end 
               
       else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE NOT THE SAME
             
             %IF DRY MODEL IS ANY OF THE GPT2w MODEL (1° OR 5° GRID) & WET MODEL IS NOT
             if all([any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14)]),...
                     ~any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14)])])     
                
                %RETRIEVE STORED VMF COEFFICIENTS
                ah1 = getappdata(0,'ah_vmf1');
                aw1 = getappdata(0,'aw_gpt2w_vmf1');
                
             %IF WET MODEL IS ANY OF THE GPT2w MODEL (1° OR 5° GRID) & DRY MODEL IS NOT    
             elseif all([~any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14)]),...
                          any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14)])])     
                
                   %RETRIEVE STORED VMF COEFFICIENTS   
                   ah1 = getappdata(0,'ah_gpt2w_vmf1');
                   aw1 = getappdata(0,'aw_vmf1');
                    
             else %IF DRY OR WET TROPO MODEL IS NOT ANY OF THE GPT2w MODEL (1° OR 5° GRID)  
                  %Call the "vmfCOE_gpt2w.m" function TO COMPUTE VMF1 COEFFICIENTS
                  [ah1,aw1] = vmfCOE_gpt2w(UTCtime,latR,lonR,hgt,grid_h,grid_resh,Timevar_h);
                 
             end  
              
       end  
 
%IF MF MODEL IS THE GPT2 MODEL(5° GRID) AND NOT ANY OF VMF1 OR GPT2w MODELS       
elseif any([strncmpi(MFh_model,'GPT2 (5° x 5°)',14),strncmpi(MFw_model,'GPT2 (5° x 5°)',14)])
        
       %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE THE SAME
       if any([strncmpi(dryModel,wetModel,14),strfind(dryModel,wetModel)]) 
          
          %IF DRY TROPO MODEL IS THE GPT2 MODEL(5° GRID) 
          if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2 (5° x 5°)',14)])
             
             %RETRIEVE STORED VMF COEFFICIENTS
             ah1 = getappdata(0,'ah_vmf1');
             aw1 = getappdata(0,'aw_vmf1');
               
          end  
               
       else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE NOT THE SAME
           
            %IF DRY MODEL IS THE GPT2 MODEL(5° GRID) & WET MODEL IS NOT
           if all([strncmpi(dryModel,'GPT2 (5° x 5°)',14),~strncmpi(wetModel,'GPT2 (5° x 5°)',14)])
               
              %RETRIEVE STORED VMF COEFFICIENTS
              ah1 = getappdata(0,'ah_vmf1');
              aw1 = getappdata(0,'aw_gpt2_vmf1'); 
                
            %IF WET MODEL IS THE GPT2 MODEL(5° GRID) & DRY MODEL IS NOT    
           elseif all([~strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2 (5° x 5°)',14)])
               
                  %RETRIEVE STORED VMF COEFFICIENTS
                  ah1 = getappdata(0,'ah_gpt2_vmf1');
                  aw1 = getappdata(0,'aw_vmf1'); 
                  
           else %IF DRY OR WET TROPO MODEL IS NOT THE GPT2 MODEL(5° GRID) 
               
                %Call the "vmfCOE_gpt2.m" function & COMPUTE VMF1 COEFFICIENTS
                [ah1,aw1] = vmfCOE_gpt2(UTCtime,latR,lonR,hgt,grid_h,Timevar_h);
                
           end
       
       end
       
   end %//if any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),strncmpi(MFh_model,'VMF ZHD',7),...
       %          strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),strncmpi(MFw_model,'VMF ZWD',7)])
    
else %IF HYDROSTATIC & WET MAPPING FUNCTION(MF) MODELS ARE NOT THE SAME
  
    %IF HYDROSTATIC MF MODEL IS ANY OF THE VMF1 OR GPT MODELS & WET MF MODEL IS NOT
    if all([any([strncmpi(MFh_model,'VMF',3),strfind(MFh_model,'VMF'),strncmpi(MFh_model,'GPT2',4),strfind(MFh_model,'GPT2')]),...
            ~any([strncmpi(MFw_model,'VMF',3),strfind(MFw_model,'VMF'),strncmpi(wetModel,'GPT2',4),strfind(wetModel,'GPT2')]),...
             any([strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12)])])
      
       %IF HYDROSTATIC MF MODEL IS ANY OF THE VMF1 MODELS
       if any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),strncmpi(MFh_model,'VMF g',5)])
           
          %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines WERE THE VMF GRID MODELS 
          if any([any([strncmpi(dryModel,'VMF gridded ZHD',15),strfind(dryModel,'VMF gridded ZHD')]),...
                  any([strncmpi(wetModel,'VMF gridded ZWD',15),strfind(wetModel,'VMF gridded ZWD')])])
              
             %FIRST GET VMF COEFFIFFIENTS FROM TropModel or TropoModel subroutines
             if all([~isempty(getappdata(0,'ah_vmfg1')),~isempty(getappdata(0,'aw_vmfg1'))])
              
                %RETRIEVE STORED VMF COEFFICIENTS 
                ah1 = getappdata(0,'ah_vmfg1');
                aw1 = getappdata(0,'aw_vmfg1');
                VMF_model = 'VMF1';
            
             else %IF STORED DATA [getappdata(0,'ah_vmfg1')& getappdata(0,'aw_vmfg1')] ARE EMPTY([])
               
                 if VMF_grid_found == 1 %IF VMF GRID FILES ARE FOUND,THEN COMPUTE VMF COE USING "SearchReadVMFgrid"
               
                    %Call the "SearchReadVMFgrid.m" function TO COMPUTE VMF1 COEFFICIENTS
                    [ah1,aw1,VMF_model] = SearchReadVMFgrid(UTCtime,latD,lonD,hgt,VMFgrids);
              
                 else %IF VMF GRID FILES ARE NOT FOUND
               
                     %ASSIGN EMPTY([]) MATRIX 
                     ah1       = [];
                     aw1       = [];
                     VMF_model = [];
                  
                 end  
        
             end 
         
          else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE NOT THE VMF GRID MODELS
               
               if VMF_grid_found == 1 %IF VMF GRID FILES ARE FOUND,THEN COMPUTE VMF COE USING "SearchReadVMFgrid"
               
                  %Call the "SearchReadVMFgrid.m" function TO COMPUTE VMF1 COEFFICIENTS
                  [ah1,aw1,VMF_model] = SearchReadVMFgrid(UTCtime,latD,lonD,hgt,VMFgrids);
              
               else %IF VMF GRID FILES ARE NOT FOUND
               
                    %ASSIGN EMPTY([]) MATRIX  
                    ah1       = [];
                    aw1       = [];
                    VMF_model = [];
                  
               end   
          
          end
          
          %IF 1ST PREAMBLE FAILS, GET COEFFICIENTS FROM GPT2w MODELS 
          %EG:.gpt2_5w.grd or gpt2_1w.grd or gpt2_1wA.grd
      
          %IF VMF GRID MODEL IS THE VMF3 MODEL OR VMF GRID FILES WERE NOT AVAILABLE
          if any([strncmpi(VMF_model,'VMF3',4),all([isempty(VMF_model),isempty(ah1),isempty(aw1)])])

             %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE THE SAME
             if any([strncmpi(dryModel,wetModel,14),strfind(dryModel,wetModel)]) 
          
                if any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                        strncmpi(dryModel,'GPT2 (5° x 5°)',14)])
             
                   %RETRIEVE STORED VMF COEFFICIENTS
                   ah1 = getappdata(0,'ah_vmf1');
                   aw1 = getappdata(0,'aw_vmf1');
                end  
            
             else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel %subroutines ARE NOT THE SAME 
          
                  %IF DRY MODEL IS ANY OF THE GPT2w MODEL 1° OR 5° GRID & WET MODEL IS NOT
                  if all([any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                               strncmpi(dryModel,'GPT2 (5° x 5°)',14)]),...
                          ~any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                                strncmpi(wetModel,'GPT2 (5° x 5°)',14)])])     
                
                     %RETRIEVE STORED VMF COEFFICIENTS
                     if any([strncmpi(getappdata(0,'Drymodel'),'GPT2',4),strfind(getappdata(0,'Drymodel'),'GPT2')])
                         
                        ah1 = getappdata(0,'ah_vmf1');
                        aw1 = getappdata(0,'aw_gpt2_vmf1');
                        
                     elseif any([strncmpi(getappdata(0,'Drymodel'),'GPT2w',5),strfind(getappdata(0,'Drymodel'),'GPT2w')])
                            
                            ah1 = getappdata(0,'ah_vmf1');
                            aw1 = getappdata(0,'aw_gpt2w_vmf1');     
                     end
                 
                  %IF WET MODEL IS ANY OF THE GPT2w MODEL (1° OR 5° GRID) & DRY MODEL IS NOT   
                  elseif all([~any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                                    strncmpi(dryModel,'GPT2 (5° x 5°)',14)]),...
                               any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                                    strncmpi(wetModel,'GPT2 (5° x 5°)',14)])])     
                
                         %RETRIEVE STORED VMF COEFFICIENTS 
                         if any([strncmpi(getappdata(0,'Wetmodel'),'GPT2',4),strfind(getappdata(0,'Wetmodel'),'GPT2')])
                         
                            ah1 = getappdata(0,'ah_gpt2_vmf1');
                            aw1 = getappdata(0,'aw_vmf1');
                            
                         elseif any([strncmpi(getappdata(0,'Wetmodel'),'GPT2w',5),strfind(getappdata(0,'Wetmodel'),'GPT2w')])
                 
                                ah1 = getappdata(0,'ah_gpt2w_vmf1');
                                aw1 = getappdata(0,'aw_vmf1');  
                         end
                         
                  else %IF DRY OR WET TROPO MODEL IS NOT ANY OF THE GPT2w MODEL (1° OR 5° GRID)
                  
                       try 
                          %CHECK LENGTH OF GPT GRID VALUES{FOR GPT2W,LENGTH SHOULD BE 10}
                          if length(grid_h) == 10
                
                             %Call the "vmfCOE_gpt2w.m" function TO COMPUTE VMF1 COEFFICIENTS
                             [ah1,aw1] = vmfCOE_gpt2w(UTCtime,latR,lonR,hgt,grid_h,grid_resh,Timevar_h);
           
                          elseif length(grid_h) == 8 %IF LENGTH IS EQUALS 8 {FOR GPT2,LENGTH SHOULD BE 8}
                  
                                 %Call the "gvmfCOE_gpt2.m" function TO COMPUTE VMF1 COEFFICIENTS
                                 [ah1,aw1] = vmfCOE_gpt2(UTCtime,latR,lonR,hgt,grid_h,Timevar_h);     
                          end     
           
                       catch    
                            %Call the "vmfCOE_gpt2.m" function
                            [ah1,aw1] = vmfCOE_gpt2(UTCtime,latR,lonR,hgt,grid_h,Timevar_h);
                       end         
                  
                  end %//if all([any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                      %               strncmpi(dryModel,'GPT2 (5° x 5°)',14)]),...
                      %         ~any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                      %               strncmpi(wetModel,'GPT2 (5° x 5°)',14)])])  
              
             end %//if any([strncmpi(drymodel,wetmodel,14),strfind(drymodel,wetmodel)]) 
         
          end %//if any([strncmpi(VMF_model,'VMF3',4),all([isempty(VMF_model),isempty(ah1),isempty(aw1)])])
        
       %IF MF MODEL IS RATHER ANY OF THE GPT2w MODEL(1° OR 5°) AND NOT ANY OF VMF1 MODELS
       elseif any([strncmpi(MFh_model,'GPT2w (1° x 1°)',14),strncmpi(MFh_model,'GPT2w (5° x 5°)',14)])
              
             %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE THE SAME
             if any([strncmpi(dryModel,wetModel,14),strfind(dryModel,wetModel)]) 
          
                %IF DRY TROPO MODEL IS ANY OF THE GPT2w MODEL(1° OR 5°)
                if any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14)])
             
                   %RETRIEVE STORED VMF COEFFICIENTS 
                   ah1 = getappdata(0,'ah_vmf1');
                   aw1 = getappdata(0,'aw_vmf1');  
                end   
               
             else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE NOT THE SAME
             
                  %IF DRY MODEL IS ANY OF THE GPT2w MODEL (1° OR 5° GRID) & WET MODEL IS NOT
                  if all([any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14)]),...
                          ~any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14)])])     
                
                     %RETRIEVE STORED VMF COEFFICIENTS
                     ah1 = getappdata(0,'ah_vmf1');
                     aw1 = getappdata(0,'aw_gpt2w_vmf1');
                
                     %IF WET MODEL IS ANY OF THE GPT2w MODEL (1° OR 5° GRID) & DRY MODEL IS NOT    
                  elseif all([~any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14)]),...
                               any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14)])])     
                
                         %RETRIEVE STORED VMF COEFFICIENTS   
                         ah1 = getappdata(0,'ah_gpt2w_vmf1');
                         aw1 = getappdata(0,'aw_vmf1');
                    
                  else %IF DRY OR WET TROPO MODEL IS NOT ANY OF THE GPT2w MODEL (1° OR 5° GRID)  
                       %Call the "vmfCOE_gpt2w.m" function
                      [ah1,aw1] = vmfCOE_gpt2w(UTCtime,latR,lonR,hgt,grid_h,grid_resh,Timevar_h);
                 
                  end   
              
             end %//if any([strncmpi(drymodel,wetmodel,14),strfind(drymodel,wetmodel)])   
 
       %IF MF MODEL IS THE GPT2 MODEL(5° GRID) AND NOT ANY OF VMF1 OR GPT2w MODELS       
       elseif any([strncmpi(MFh_model,'GPT2 (5° x 5°)',14),strncmpi(MFw_model,'GPT2 (5° x 5°)',14)])
        
              %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE THE SAME
              if any([strncmpi(dryModel,wetModel,14),strfind(dryModel,wetModel)]) 
          
                 %IF DRY TROPO MODEL IS THE GPT2 MODEL(5° GRID) 
                 if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2 (5° x 5°)',14)])
             
                    %RETRIEVE STORED VMF COEFFICIENTS
                    ah1 = getappdata(0,'ah_vmf1');
                    aw1 = getappdata(0,'aw_vmf1');
               
                 end   
               
              else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE NOT THE SAME
           
                   %IF DRY MODEL IS THE GPT2 MODEL(5° GRID) & WET MODEL IS NOT
                   if all([strncmpi(dryModel,'GPT2 (5° x 5°)',14),~strncmpi(wetModel,'GPT2 (5° x 5°)',14)])
               
                      %RETRIEVE STORED VMF COEFFICIENTS
                      ah1 = getappdata(0,'ah_vmf1');
                      aw1 = getappdata(0,'aw_gpt2_vmf1'); 
                
                   %IF WET MODEL IS THE GPT2 MODEL(5° GRID) & DRY MODEL IS NOT    
                   elseif all([~strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2 (5° x 5°)',14)])
               
                          %RETRIEVE STORED VMF COEFFICIENTS
                          ah1 = getappdata(0,'ah_gpt2_vmf1');
                          aw1 = getappdata(0,'aw_vmf1'); 
                  
                   else %IF DRY OR WET TROPO MODEL IS NOT THE GPT2 MODEL(5° GRID) 
               
                        %Call the "vmfCOE_gpt2.m" function & COMPUTE VMF1 COEFFICIENTS
                        [ah1,aw1] = vmfCOE_gpt2(UTCtime,latR,lonR,hgt,grid_h,Timevar_h);
                
                   end 
       
              end %//if any([strncmpi(drymodel,wetmodel,14),strfind(drymodel,wetmodel)])
       
       end  %//if any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),strncmpi(MFh_model,'VMF ZHD',7),...
            %          strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),strncmpi(MFw_model,'VMF ZWD',7)]) 
        
            
    %IF HYDROSTATIC MF MODEL IS NOT ANY OF THE VMF1 OR GPT MODELS BUT WET MF MODEL IS
    elseif all([any([~strncmpi(MFh_model,'VMF',3),~strfind(MFh_model,'VMF'),~strncmpi(MFh_model,'GPT2',4),~strfind(MFh_model,'GPT2'),...
                     any([~strncmpi(MFh_model,'VMF1(1° x 1°)',12),~strncmpi(MFh_model,'VMF1(5° x 5°)',12)])]),...
                any([strncmpi(MFw_model,'VMF',3),strfind(MFw_model,'VMF'),strncmpi(wetModel,'GPT2',4),strfind(wetModel,'GPT2')])]) 
           
           %IF WET MF MODEL IS ANY OF THE VMF1 MODELS
           if any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),strncmpi(MFw_model,'VMF g',5)])
           
              %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines WERE THE VMF GRID MODELS 
              if any([any([strncmpi(dryModel,'VMF gridded ZHD',15),strfind(dryModel,'VMF gridded ZHD')]),...
                      any([strncmpi(wetModel,'VMF gridded ZWD',15),strfind(wetModel,'VMF gridded ZWD')])])
              
                 %FIRST GET VMF COEFFIFFIENTS FROM TropModel or TropoModel subroutines
                 if all([~isempty(getappdata(0,'ah_vmfg1')),~isempty(getappdata(0,'aw_vmfg1'))])
              
                    %RETRIEVE STORED VMF COEFFICIENTS 
                    ah1 = getappdata(0,'ah_vmfg1');
                    aw1 = getappdata(0,'aw_vmfg1');
                    VMF_model = 'VMF1';
            
                 else %IF STORED DATA [getappdata(0,'ah_vmfg1')& getappdata(0,'aw_vmfg1')] ARE EMPTY([])
               
                     if VMF_grid_found == 1 %IF VMF GRID FILES ARE FOUND,THEN COMPUTE VMF COE USING "SearchReadVMFgrid"
               
                        %Call the "SearchReadVMFgrid.m" function TO COMPUTE VMF1 COEFFICIENTS
                        [ah1,aw1,VMF_model] = SearchReadVMFgrid(UTCtime,latD,lonD,hgt,VMFgrids);
              
                     else %IF VMF GRID FILES ARE NOT FOUND
               
                         %ASSIGN EMPTY([]) MATRIX   
                         ah1       = [];
                         aw1       = [];
                         VMF_model = [];
                  
                     end    
        
                 end  
         
              else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE NOT THE VMF GRID MODELS
               
                   if VMF_grid_found == 1 %IF VMF GRID FILES ARE FOUND,THEN COMPUTE VMF COE USING "SearchReadVMFgrid"
               
                      %Call the "SearchReadVMFgrid.m" function TO COMPUTE VMF1 COEFFICIENTS
                      [ah1,aw1,VMF_model] = SearchReadVMFgrid(UTCtime,latD,lonD,hgt,VMFgrids);
              
                   else %IF VMF GRID FILES ARE NOT FOUND
               
                       %ASSIGN EMPTY([]) MATRIX   
                       ah1       = [];
                       aw1       = [];
                       VMF_model = [];
                  
                   end    

              end 
          
          %IF 1ST PREAMBLE FAILS, GET COEFFICIENTS FROM GPT2w MODELS 
          %EG:.gpt2_5w.grd or gpt2_1w.grd or gpt2_1wA.grd
      
          %IF VMF GRID MODEL IS THE VMF3 MODEL OR VMF GRID FILES WERE NOT AVAILABLE
          if any([strncmpi(VMF_model,'VMF3',4),all([isempty(VMF_model),isempty(ah1),isempty(aw1)])])

             %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE THE SAME
             if any([strncmpi(dryModel,wetModel,14),strfind(dryModel,wetModel)]) 
          
                if any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                        strncmpi(dryModel,'GPT2 (5° x 5°)',14)])
             
                   %RETRIEVE STORED VMF COEFFICIENTS
                   ah1 = getappdata(0,'ah_vmf1');
                   aw1 = getappdata(0,'aw_vmf1');
                   
                else
                   
                     try 
                        %CHECK LENGTH OF GPT GRID VALUES{FOR GPT2W,LENGTH SHOULD BE 10}
                        if length(grid_w) == 10
                
                           %Call the "vmfCOE_gpt2w.m" function
                           [ah1,aw1] = vmfCOE_gpt2w(UTCtime,latR,lonR,hgt,grid_w,grid_resw,Timevar_w);
           
                        elseif length(grid_w) == 8 %IF LENGTH IS EQUALS 8 {FOR GPT2,LENGTH SHOULD BE 8}
                  
                                 %Call the "gvmfCOE_gpt2.m" function TO COMPUTE VMF1 COEFFICIENTS
                                 [ah1,aw1] = vmfCOE_gpt2(UTCtime,latR,lonR,hgt,grid_w,Timevar_w);     
                        end      
           
                     catch      
                           %Call the "vmfCOE_gpt2.m" function TO COMPUTE VMF1 COEFFICIENTS
                           [ah1,aw1] = vmfCOE_gpt2(UTCtime,latR,lonR,hgt,grid_w,Timevar_w);
                     end           
                       
                end  
            
             else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel %subroutines ARE NOT THE SAME 
          
                  %IF DRY MODEL IS ANY OF THE GPT2w MODEL 1° OR 5° GRID & WET MODEL IS NOT
                  if all([any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                               strncmpi(dryModel,'GPT2 (5° x 5°)',14)]),...
                          ~any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                                strncmpi(wetModel,'GPT2 (5° x 5°)',14)])])     
                
                     %RETRIEVE STORED VMF COEFFICIENTS
                     if any([strncmpi(getappdata(0,'Drymodel'),'GPT2',4),strfind(getappdata(0,'Drymodel'),'GPT2')])
                         
                        ah1 = getappdata(0,'ah_vmf1');
                        aw1 = getappdata(0,'aw_gpt2_vmf1');
                        
                     elseif any([strncmpi(getappdata(0,'Drymodel'),'GPT2w',5),strfind(getappdata(0,'Drymodel'),'GPT2w')])
                            
                            ah1 = getappdata(0,'ah_vmf1');
                            aw1 = getappdata(0,'aw_gpt2w_vmf1');     
                     end
                 
                  %IF WET MODEL IS ANY OF THE GPT2w MODEL (1° OR 5° GRID) & DRY MODEL IS NOT   
                  elseif all([~any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                                    strncmpi(dryModel,'GPT2 (5° x 5°)',14)]),...
                               any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                                    strncmpi(wetModel,'GPT2 (5° x 5°)',14)])])     
                
                         %RETRIEVE STORED VMF COEFFICIENTS           
                         if any([strncmpi(getappdata(0,'Wetmodel'),'GPT2',4),strfind(getappdata(0,'Wetmodel'),'GPT2')])
                         
                            ah1 = getappdata(0,'ah_gpt2_vmf1');
                            aw1 = getappdata(0,'aw_vmf1');
                            
                         elseif any([strncmpi(getappdata(0,'Wetmodel'),'GPT2w',5),strfind(getappdata(0,'Wetmodel'),'GPT2w')])
                 
                                ah1 = getappdata(0,'ah_gpt2w_vmf1');
                                aw1 = getappdata(0,'aw_vmf1');  
                         end
                 
                  else %IF DRY OR WET TROPO MODEL IS NOT ANY OF THE GPT2w MODEL (1° OR 5° GRID)
                  
                       try 
                          %CHECK LENGTH OF GPT GRID VALUES{FOR GPT2W,LENGTH SHOULD BE 10}
                          if length(grid_w) == 10
                
                             %Call the "vmfCOE_gpt2w.m" function
                             [ah1,aw1] = vmfCOE_gpt2w(UTCtime,latR,lonR,hgt,grid_w,grid_resw,Timevar_w);
           
                          elseif length(grid_w) == 8 %IF LENGTH IS EQUALS 8 {FOR GPT2,LENGTH SHOULD BE 8}
                  
                                 %Call the "gvmfCOE_gpt2.m" function TO COMPUTE VMF1 COEFFICIENTS
                                 [ah1,aw1] = vmfCOE_gpt2(UTCtime,latR,lonR,hgt,grid_w,Timevar_w);     
                          end     
           
                       catch    
                            %Call the "vmfCOE_gpt2.m" function TO COMPUTE VMF1 COEFFICIENTS
                            [ah1,aw1] = vmfCOE_gpt2(UTCtime,latR,lonR,hgt,grid_w,Timevar_w);
                       end         
                  
                  end %//if all([any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                      %               strncmpi(dryModel,'GPT2 (5° x 5°)',14)]),...
                      %          ~any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                      %                strncmpi(wetModel,'GPT2 (5° x 5°)',14)])])  
              
             end %//if any([strncmpi(drymodel,wetmodel,14),strfind(drymodel,wetmodel)]) 
         
          end %//if any([strncmpi(VMF_model,'VMF3',4),all([isempty(VMF_model),isempty(ah1),isempty(aw1)])])
        
       %IF WET MF MODEL IS RATHER ANY OF THE GPT2w MODEL(1° OR 5°) AND NOT ANY OF VMF1 MODELS
       elseif any([strncmpi(MFw_model,'GPT2w (1° x 1°)',14),strncmpi(MFw_model,'GPT2w (5° x 5°)',14)])
              
             %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE THE SAME
             if any([strncmpi(dryModel,wetModel,14),strfind(dryModel,wetModel)]) 
          
                %IF DRY TROPO MODEL IS ANY OF THE GPT2w MODEL(1° OR 5°)
                if any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14)])
             
                   %RETRIEVE STORED VMF COEFFICIENTS 
                   ah1 = getappdata(0,'ah_vmf1');
                   aw1 = getappdata(0,'aw_vmf1');  
                end   
               
             else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE NOT THE SAME
             
                  %IF DRY MODEL IS ANY OF THE GPT2w MODEL (1° OR 5° GRID) & WET MODEL IS NOT
                  if all([any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14)]),...
                          ~any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14)])])     
                
                     %RETRIEVE STORED VMF COEFFICIENTS
                     ah1 = getappdata(0,'ah_vmf1');
                     aw1 = getappdata(0,'aw_gpt2w_vmf1');
                
                     %IF WET MODEL IS ANY OF THE GPT2w MODEL (1° OR 5° GRID) & DRY MODEL IS NOT    
                  elseif all([~any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14)]),...
                               any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14)])])     
                
                         %RETRIEVE STORED VMF COEFFICIENTS   
                         ah1 = getappdata(0,'ah_gpt2w_vmf1');
                         aw1 = getappdata(0,'aw_vmf1');
                    
                  else %IF DRY OR WET TROPO MODEL IS NOT ANY OF THE GPT2w MODEL (1° OR 5° GRID)  
                       %Call the "vmfCOE_gpt2w.m" function TO COMPUTE VMF1 COEFFICIENTS
                      [ah1,aw1] = vmfCOE_gpt2w(UTCtime,latR,lonR,hgt,grid_w,grid_resw,Timevar_w);
                 
                  end   
              
             end %//if any([strncmpi(drymodel,wetmodel,14),strfind(drymodel,wetmodel)])   
 
       %IF WET MF MODEL IS THE GPT2 MODEL(5° GRID) AND NOT ANY OF THE VMF1 OR GPT2w MODELS       
       elseif strncmpi(MFw_model,'GPT2 (5° x 5°)',14)
        
              %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE THE SAME
              if any([strncmpi(dryModel,wetModel,14),strfind(dryModel,wetModel)]) 
          
                 %IF DRY TROPO MODEL IS THE GPT2 MODEL(5° GRID) 
                 if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2 (5° x 5°)',14)])
             
                    %RETRIEVE STORED VMF COEFFICIENTS
                    ah1 = getappdata(0,'ah_vmf1');
                    aw1 = getappdata(0,'aw_vmf1');
               
                 end   
               
              else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE NOT THE SAME
           
                   %IF DRY MODEL IS THE GPT2 MODEL(5° GRID) & WET MODEL IS NOT
                   if all([strncmpi(dryModel,'GPT2 (5° x 5°)',14),~strncmpi(wetModel,'GPT2 (5° x 5°)',14)])
               
                      %RETRIEVE STORED VMF COEFFICIENTS
                      ah1 = getappdata(0,'ah_vmf1');
                      aw1 = getappdata(0,'aw_gpt2_vmf1'); 
                
                   %IF WET MODEL IS THE GPT2 MODEL(5° GRID) & DRY MODEL IS NOT    
                   elseif all([~strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2 (5° x 5°)',14)])
               
                          %RETRIEVE STORED VMF COEFFICIENTS
                          ah1 = getappdata(0,'ah_gpt2_vmf1');
                          aw1 = getappdata(0,'aw_vmf1'); 
                  
                   else %IF DRY OR WET TROPO MODEL IS NOT THE GPT2 MODEL(5° GRID) 
               
                        %Call the "vmfCOE_gpt2.m" function & COMPUTE VMF1 COEFFICIENTS
                        [ah1,aw1] = vmfCOE_gpt2(UTCtime,latR,lonR,hgt,grid_w,Timevar_w);
                
                   end 
       
              end %//if any([strncmpi(drymodel,wetmodel,14),strfind(drymodel,wetmodel)])
       
           end   %//if any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),strncmpi(MFw_model,'VMF Z',5)]) 
         
         
    end %//if all([any([strncmpi(MFh_model,'VMF',3),strfind(MFh_model,'VMF'),strncmpi(MFh_model,'GPT2',4),strfind(MFh_model,'GPT2')]),...
        %          any([~strncmpi(MFw_model,'VMF',3),~strfind(MFw_model,'VMF'),~strncmpi(Wetmodel,'GPT2',4),~strfind(Wetmodel,'GPT2')])])
            
            
end %//if any([strncmpi(MFh_model,MFw_model,12),strfind(MFh_model,MFw_model),any([strncmpi(MFh_model,'VMF ZHD',7),...
    %          strncmpi(MFw_model,'VMF ZWD',7)])])
        
%******************************************END OF VMF1coe.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%C.4.SUBROUTINE TO COMPUTE VIENNA MAPPING FUNCTION 3(VMF3) COEFFICIENTS
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [ah3,aw3] = VMF3coe(UTCtime,latD,lonD,hgt,VMFgrids,grid_h,...
                             grid_w,grid_resh,grid_resw,Timevar_h,...
                                             Timevar_w,MFh_model,MFw_model)

%**************************************************************************
%***DESCRIPTION:
%              This subroutine determines the Tropospheric Mapping         *
%              Function using the Vienna Mapping Functions 3(VMF3)         * 
%              coefficients(ah,aw) from various VMF1 Models(i.e. either... *
%              from VMF1 grid file(s) or from GPT3 models for Specific site* 
%              given receiver/station position in [lat Lon hgt], receiver  *
%              reception time in the utc time format[Yr Mon Day Hr Min Sec]*
%              ,GPT3 grid values,grid resolution  and time variation      *
%              indicator(0 or 1)                                           *
%                                                                          *                                       
%--------------------------------------------------------------------------*                                     
%USAGE:                                                                    *                                                           *
%      [ah3,aw3] = VMF3coe(UTCtime,latD,lonD,hgt,VMFgrids,grid_h,...       *
%                             grid_w,grid_resh,grid_resw,Timevar_h,...     *
%                                            Timevar_w,MFh_model,MFw_model)*       
%--------------------------------------------------------------------------*
%****INPUTs:                                                               *
%1.    UTCtime : UTC time in [Year,Month,Day,Hour,Minute,Seconds]          *
%2.       latD : station ellipsoidal latitude in [degrees]                 *
%3.       lonD : station ellipsoidal longitude in [degrees]                *
%4.        hgt : station ellipsoidal height in [meters]                    *
%5.VMFgrids    : Structure array of VMF grid files which can be :          *
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
%--------------------------------------------------------------------------*
%6.      grid   : grid values in cells extracted from the grid file        *
%    grid_res   : grid resolution (°) (possible: 1 or 5)                   *
%7    Timevar   : case 1: no time variation but static quantities          *
%                 case 0: with time variation (annual and semiannual terms)*
%8.  MF_model   :  mapping function model type(i.e.GPT2 OR GPT2w)          *
%--------------------------------------------------------------------------*
%***OUTPUTs:                                                               *
%          ah3  : Vienna Hydrostatic mapping function 3 coefficient        *
%          aw3  : Vienna Wet mapping function  3 coefficient               *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%WRITTEN by:  Samuel Osah,Msc Geomatic Engineering ,2016                   *    
%             Email: osahsamuel@yahoo.ca                                   *
%             Tel:+233 (0)246137410/+233 (0)509438484                      *        
%==========================================================================+
%**************************************************************************+
%**************************************************************************+

%***CONVERT Lat and Lon IN DEGREES TO RADIAN 
latR = (latD./180).*pi; %[radian]
lonR = (lonD./180).*pi ;%[radian]    

%RETRIEVE SOME SAVED DATA FROM TropModel or TropoModel subroutines
%GET DRY & WET MODELs 
dryModel = getappdata(0,'dryModel');%get selected Hydrostatic model
wetModel = getappdata(0,'wetModel');%get selected Wet model

%GET VMF grid FILE FOUND/NOT FOUND INDICATOR
VMF_grid_found = getappdata(0,'VMF_grid_found');
                 
%*****************1ST PREAMBLE
%SEARCH VMF1 GRID FILEs IN DIRECTORY & COMPUTE MF COEFFICIENTs
%USING "SearchReadVMFgrid.m" FUNCTION
%NOTE:
% "SearchReadVMFgrid.m" will search and read the following grid file type : 
%1.VMF3_yrmonthday.H00 or VMF3_yrmonthday.H06 or VMF3_yrmonthday.H12 or
%  VMF3_yrmonthday.H18 eg:[VMF3_20180101.H00, VMF3_20180101.H06,... 
%                          VMF3_20180101. H12,VMF3_20180101.H18]

%NOTE: 
%     GENERALLY, VMF3 COEFFICIENTS ARE EXTRACTED OR COMPUTED FROM VMF3 GRID FILES.
%     WHERE VMF3 GRID FILES ARE NOT AVAILABLE, WE THEN RESORT TO THOSE FROM THE GPT3(1° or 5°) MODELS

%IF HYDROSTATIC & WET MAPPING FUNCTION(MF) MODELS ARE THE SAME
if any([strncmpi(MFh_model,MFw_model,12),strfind(MFh_model,MFw_model),any([strncmpi(MFh_model,'VMF gridded ZHD',15),...
        strncmpi(MFw_model,'VMF gridded ZWD',15)])])
    
    %IF MF MODEL IS ANY OF THE VIENNA MAPPING FUNCTION 1(VMF1) MODELS
   if any([strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12),strncmpi(MFh_model,'VMF gridded ZHD',15),...
           strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12),strncmpi(MFw_model,'VMF gridded ZWD',15)])
    
       %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines WERE THE VMF GRID MODELS
      if any([any([strncmpi(dryModel,'VMF gridded ZHD',15),strfind(dryModel,'VMF gridded ZHD')])
              any([strncmpi(wetModel,'VMF gridded ZWD',15),strfind(wetModel,'VMF gridded ZWD')])])
          
         %FIRST GET VMF3 COEFFIFFIENTS FROM TropModel or TropoModel subroutines 
         if all([~isempty(getappdata(0,'ah_vmfg3')),~isempty(getappdata(0,'aw_vmfg3'))])
            
            %RETRIEVE STORED VMF3 COEFFICIENTS 
            ah3 = getappdata(0,'ah_vmfg3');
            aw3 = getappdata(0,'aw_vmfg3');
            VMF_model = 'VMF3';
            
         else %IF STORED DATA [getappdata(0,'ah_vmfg3')& getappdata(0,'aw_vmfg3')] ARE EMPTY([])
               
              if VMF_grid_found == 1 %IF VMF GRID FILES ARE FOUND,THEN COMPUTE VMF COE USING "SearchReadVMFgrid"
               
                 %Call the "SearchReadVMFgrid.m" function TO COMPUTE VMF3 COEFFICIENTS
                 [ah3,aw3,VMF_model] = SearchReadVMFgrid(UTCtime,latD,lonD,hgt,VMFgrids);
              
              else %IF VMF GRID FILES ARE NOT FOUND
               
                   %ASSIGN EMPTY([]) MATRIX     
                   ah3       = [];
                   aw3       = [];
                   VMF_model = [];
                  
              end %//if VMF_grid_found == 1    
        
         end
         
      else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines WERE NOT THE VMF GRID MODELS 
          
           if VMF_grid_found == 1 %IF VMF GRID FILES ARE FOUND,THEN COMPUTE VMF COE USING "SearchReadVMFgrid"
               
              %Call the "SearchReadVMFgrid.m" function TO COMPUTE VMF3 COEFFICIENTS
              [ah3,aw3,VMF_model] = SearchReadVMFgrid(UTCtime,latD,lonD,hgt,VMFgrids);
              
           else %IF VMF GRID FILES ARE NOT FOUND
               
                %ASSIGN EMPTY([]) MATRIX     
                ah3       = [];
                aw3       = [];
                VMF_model = [];
                  
           end  %//if VMF_grid_found == 1
          
      end

      %IF 1ST PREAMBLE FAILS, GET COEFFICIENTS FROM GPT3 MODELS 
      %EG:.gpt3_5.grd or gpt3_1.grd 
      
      %IF VMF GRID MODEL IS THE VMF1 MODEL OR VMF3 GRID FILES WERE NOT AVAILABLE
      if any([strncmpi(VMF_model,'VMF1',4),all([isempty(VMF_model),isempty(ah3),isempty(aw3)])])

          %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel
          %subroutines ARE THE SAME
         if any([strncmpi(dryModel,wetModel,14),strfind(dryModel,wetModel)]) 
          
            if any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)])
             
               %RETRIEVE STORED VMF3 COEFFICIENTS
               ah3 = getappdata(0,'ah_vmf3');
               aw3 = getappdata(0,'aw_vmf3');
               
            end 
            
         else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel %subroutines ARE NOT THE SAME 
          
              %IF DRY MODEL IS ANY OF THE GPT2w MODEL 1° OR 5° GRID & WET MODEL IS NOT
              if all([any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),...
                      ~any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])])     
                
                 %RETRIEVE STORED VMF3 COEFFICIENTS
                 if any([strncmpi(getappdata(0,'Drymodel'),'GPT3',4),strfind(getappdata(0,'Drymodel'),'GPT3')])
                    
                    ah3 = getappdata(0,'ah_vmf3');
                    aw3 = getappdata(0,'aw_gpt3_vmf3');
                    
                 end 
                 
              %IF WET MODEL IS ANY OF THE GPT2w MODEL (1° OR 5° GRID) & DRY MODEL IS NOT   
              elseif all([~any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),...
                           any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])])     
                
                 %RETRIEVE STORED VMF3 COEFFICIENTS
                 if any([strncmpi(getappdata(0,'Wetmodel'),'GPT3',4),strfind(getappdata(0,'Wetmodel'),'GPT3')])
                    
                    ah3 = getappdata(0,'ah_gpt3_vmf3');
                    aw3 = getappdata(0,'aw_vmf3');
                          
                 end
                 
              else %IF DRY OR WET TROPO MODEL IS NOT ANY OF THE GPT2w MODEL (1° OR 5° GRID)
                  
                    %CHECK LENGTH OF GPT GRID VALUES{FOR GPT3,LENGTH SHOULD BE 14}
                    if any([length(grid_h) == 14,length(grid_w) == 14])
                
                       %Call the "vmfCOE_gpt3.m" function TO COMPUTE VMF3 COEFFICIENTS
                       [ah3,aw3] = vmfCOE_gpt3(UTCtime,latR,lonR,hgt,grid_h,grid_resh,Timevar_h);
           
                    end       
                  
              end  %//if all([any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),...
                   %          ~any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])])
              
         end %//if any([strncmpi(drymodel,wetmodel,14),strfind(drymodel,wetmodel)]) 
         
      end  %//if any([strncmpi(VMF_model,'VMF1',4),all([isempty(VMF_model),isempty(ah3),isempty(aw3)])])
 

%IF MF MODEL IS RATHER ANY OF THE GPT3 MODEL(1° OR 5°) AND NOT ANY OF VMF3 MODELS
elseif any([strncmpi(MFh_model,'GPT3 (1° x 1°)',14),strncmpi(MFh_model,'GPT3 (5° x 5°)',14)]),...
       any([strncmpi(MFw_model,'GPT3 (1° x 1°)',14),strncmpi(MFw_model,'GPT3 (5° x 5°)',14)])
       
       %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE THE SAME
       if any([strncmpi(dryModel,wetModel,14),strfind(dryModel,wetModel)]) 
          
            %IF DRY TROPO MODEL IS ANY OF THE GPT3 MODEL(1° OR 5°)
            if any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)])
             
               %RETRIEVE STORED VMF3 COEFFICIENTS 
               ah3 = getappdata(0,'ah_vmf3');
               aw3 = getappdata(0,'aw_vmf3');  
            end 
               
       else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE NOT THE SAME
             
             %IF DRY MODEL IS ANY OF THE GPT3 MODEL (1° OR 5° GRID) & WET MODEL IS NOT
             if all([any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),...
                     ~any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])])     
                
                %RETRIEVE STORED VMF3 COEFFICIENTS
                ah3 = getappdata(0,'ah_vmf3');
                aw3 = getappdata(0,'aw_gpt3_vmf3');
                
             %IF WET MODEL IS ANY OF THE GPT3 MODEL (1° OR 5° GRID) & DRY MODEL IS NOT    
             elseif all([~any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),...
                          any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])])     
                
                   %RETRIEVE STORED VMF3 COEFFICIENTS   
                   ah3 = getappdata(0,'ah_gpt3_vmf3');
                   aw3 = getappdata(0,'aw_vmf3');
                    
             else %IF DRY OR WET TROPO MODEL IS NOT ANY OF THE GPT3 MODEL (1° OR 5° GRID)  
                  %Call the "vmfCOE_gpt3.m" function TO COMPUTE VMF3 COEFFICIENTS
                  [ah3,aw3] = vmfCOE_gpt3(UTCtime,latR,lonR,hgt,grid_h,grid_resh,Timevar_h);
                 
             end %//if all([any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),...
                 %         ~any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])])
                     
              
       end %//if any([strncmpi(drymodel,wetmodel,14),strfind(drymodel,wetmodel)])  
 
  
   end %//if any([strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12),strncmpi(MFh_model,'VMF ZHD',7),...
       %          strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12),strncmpi(MFw_model,'VMF ZWD',7)])
    
else %IF HYDROSTATIC & WET MAPPING FUNCTION(MF) MODELS ARE NOT THE SAME
    
    %IF HYDROSTATIC MF MODEL IS ANY OF THE VMF1 OR GPT3 MODELS & WET MF MODEL IS NOT
    if all([any([strncmpi(MFh_model,'VMF',3),strfind(MFh_model,'VMF'),strncmpi(MFh_model,'GPT3',4),strfind(MFh_model,'GPT3')]),...
            any([~strncmpi(MFw_model,'VMF',3),~strfind(MFw_model,'VMF'),~strncmpi(wetModel,'GPT3',4),~strfind(wetModel,'GPT3')])])
       
       %IF HYDROSTATIC MF MODEL IS ANY OF THE VMF3 MODELS
       if any([strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12),strncmpi(MFh_model,'VMF g',5)])
           
          %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines WERE THE VMF GRID MODELS 
          if any([any([strncmpi(dryModel,'VMF gridded ZHD',15),strfind(dryModel,'VMF gridded ZHD')]),...
                  any([strncmpi(wetModel,'VMF gridded ZWD',15),strfind(wetModel,'VMF gridded ZWD')])])
              
             %FIRST GET VMF3 COEFFIFFIENTS FROM TropModel or TropoModel subroutines
             if all([~isempty(getappdata(0,'ah_vmfg3')),~isempty(getappdata(0,'aw_vmfg3'))])
              
                %RETRIEVE STORED VMF3 COEFFICIENTS 
                ah3 = getappdata(0,'ah_vmfg3');
                aw3 = getappdata(0,'aw_vmfg3');
                VMF_model = 'VMF3';
            
             else %IF STORED DATA [getappdata(0,'ah_vmfg3')& getappdata(0,'aw_vmfg3')] ARE EMPTY([])
               
                  if VMF_grid_found == 1 %IF VMF GRID FILES ARE FOUND,THEN COMPUTE VMF COE USING "SearchReadVMFgrid"
               
                    %Call the "SearchReadVMFgrid.m" function TO COMPUTE VMF3 COEFFICIENTS
                    [ah3,aw3,VMF_model] = SearchReadVMFgrid(UTCtime,latD,lonD,hgt,VMFgrids);
              
                  else %IF VMF GRID FILES ARE NOT FOUND
               
                       %ASSIGN EMPTY([]) MATRIX     
                       ah3       = [];
                       aw3       = [];
                       VMF_model = [];
                  
                  end %//if VMF_grid_found == 1
        
             end 
         
          else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE NOT THE VMF GRID MODELS
               
               if VMF_grid_found == 1 %IF VMF GRID FILES ARE FOUND,THEN COMPUTE VMF COE USING "SearchReadVMFgrid"
               
                 %Call the "SearchReadVMFgrid.m" function TO COMPUTE VMF3 COEFFICIENTS
                 [ah3,aw3,VMF_model] = SearchReadVMFgrid(UTCtime,latD,lonD,hgt,VMFgrids);
              
               else %IF VMF GRID FILES ARE NOT FOUND
               
                    %ASSIGN EMPTY([]) MATRIX    
                    ah3       = [];
                    aw3       = [];
                    VMF_model = [];
                  
               end %//if VMF_grid_found == 1
          
          end
          
          %IF 1ST PREAMBLE FAILS, GET COEFFICIENTS FROM GPT3 MODELS 
          %EG:.gpt3_5.grd or gpt3_1.grd 
      
          %IF VMF GRID MODEL IS THE VMF1 MODEL OR VMF3 GRID FILES WERE NOT AVAILABLE
          if any([strncmpi(VMF_model,'VMF1',4),all([isempty(VMF_model),isempty(ah3),isempty(aw3)])])

             %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE THE SAME
             if any([strncmpi(dryModel,wetModel,14),strfind(dryModel,wetModel)]) 
          
                if any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)])
             
                   %RETRIEVE STORED VMF3 COEFFICIENTS
                   ah3 = getappdata(0,'ah_vmf3');
                   aw3 = getappdata(0,'aw_vmf3');
                end  
            
             else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel %subroutines ARE NOT THE SAME 
          
                  %IF DRY MODEL IS ANY OF THE GPT3 MODEL 1° OR 5° GRID & WET MODEL IS NOT
                  if all([any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),...
                         ~any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])])     
                
                     %RETRIEVE STORED VMF3 COEFFICIENTS
                     if any([strncmpi(getappdata(0,'Drymodel'),'GPT3',4),strfind(getappdata(0,'Drymodel'),'GPT3')])
                         
                        ah3 = getappdata(0,'ah_vmf3');
                        aw3 = getappdata(0,'aw_gpt3_vmf3');
                           
                     end
                 
                  %IF WET MODEL IS ANY OF THE GPT3 MODEL (1° OR 5° GRID) & DRY MODEL IS NOT   
                  elseif all([~any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),...
                               any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])])     
                
                         %RETRIEVE STORED VMF3 COEFFICIENTS 
                         if any([strncmpi(getappdata(0,'Wetmodel'),'GPT3',4),strfind(getappdata(0,'Wetmodel'),'GPT3')])
                         
                            ah3 = getappdata(0,'ah_gpt3_vmf3');
                            aw3 = getappdata(0,'aw_vmf3');
                            
                         end
                         
                  else %IF DRY OR WET TROPO MODEL IS NOT ANY OF THE GPT3 MODEL (1° OR 5° GRID)
                  
                       %CHECK LENGTH OF GPT3 GRID VALUES{FOR GPT3,LENGTH SHOULD BE 14}
                       if length(grid_h) == 14
                
                          %Call the "vmfCOE_gpt3.m" function TO COMPUTE VMF3 COEFFICIENTS
                          [ah3,aw3] = vmfCOE_gpt3(UTCtime,latR,lonR,hgt,grid_h,grid_resh,Timevar_h);
           
                       end      
             
                  end %//if if all([any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),...
                      %            ~any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])])
                         
              
             end %//if any([strncmpi(drymodel,wetmodel,14),strfind(drymodel,wetmodel)]) 
         
          end %//if any([strncmpi(VMF_model,'VMF1',4),all([isempty(VMF_model),isempty(ah3),isempty(aw3)])])
        
       %IF MF MODEL IS RATHER ANY OF THE GPT3 MODEL(1° OR 5°) AND NOT ANY OF VMF3 MODELS
       elseif any([strncmpi(MFh_model,'GPT3 (1° x 1°)',14),strncmpi(MFh_model,'GPT3 (5° x 5°)',14)])
              
             %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE THE SAME
             if any([strncmpi(dryModel,wetModel,14),strfind(dryModel,wetModel)]) 
          
                %IF DRY TROPO MODEL IS ANY OF THE GPT3 MODEL(1° OR 5°)
                if any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)])
             
                   %RETRIEVE STORED VMF3 COEFFICIENTS 
                   ah3 = getappdata(0,'ah_vmf3');
                   aw3 = getappdata(0,'aw_vmf3');  
                end   
               
             else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE NOT THE SAME
             
                  %IF DRY MODEL IS ANY OF THE GPT3 MODEL (1° OR 5° GRID) & WET MODEL IS NOT
                  if all([any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),...
                         ~any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])])     
                
                     %RETRIEVE STORED VMF3 COEFFICIENTS
                     ah3 = getappdata(0,'ah_vmf3');
                     aw3 = getappdata(0,'aw_gpt3_vmf3');
                
                  else %IF DRY OR WET TROPO MODEL IS NOT ANY OF THE GPT3 MODEL (1° OR 5° GRID)  
                       %Call the "vmfCOE_gpt3.m" function TO COMPUTE VMF3 COEFFICIENTS
                       [ah3,aw3] = vmfCOE_gpt3(UTCtime,latR,lonR,hgt,grid_h,grid_resh,Timevar_h);
                 
                  end   
              
             end %//if any([strncmpi(drymodel,wetmodel,14),strfind(drymodel,wetmodel)])   
 
       
       end  %//if any([strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12),strncmpi(MFh_model,'VMF ZHD',7),...
            %          strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12),strncmpi(MFw_model,'VMF ZWD',7)]) 
        
            
    %IF HYDROSTATIC MF MODEL IS NOT ANY OF THE VMF3 OR GPT MODELS BUT WET MF MODEL IS
    elseif all([any([~strncmpi(MFh_model,'VMF',3),~strfind(MFh_model,'VMF'),~strncmpi(MFh_model,'GPT3',4),~strfind(MFh_model,'GPT3')]),...
                any([strncmpi(MFw_model,'VMF',3),strfind(MFw_model,'VMF'),strncmpi(wetModel,'GPT3',4),strfind(wetModel,'GPT3')])]) 
            
           %IF WET MF MODEL IS ANY OF THE VMF3 MODELS
           if any([strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12),strncmpi(MFw_model,'VMF g',5)])
           
              %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines WERE THE VMF GRID MODELS 
              if any([any([strncmpi(dryModel,'VMF gridded ZHD',15),strfind(dryModel,'VMF gridded ZHD')]),...
                      any([strncmpi(wetModel,'VMF gridded ZWD',15),strfind(wetModel,'VMF gridded ZWD')])])
              
                 %FIRST GET VMF COEFFIFFIENTS FROM TropModel or TropoModel subroutines
                 if all([~isempty(getappdata(0,'ah_vmfg3')),~isempty(getappdata(0,'aw_vmfg3'))])
              
                    %RETRIEVE STORED VMF COEFFICIENTS 
                    ah3 = getappdata(0,'ah_vmfg3');
                    aw3 = getappdata(0,'aw_vmfg3');
                    VMF_model = 'VMF3';
            
                 else %IF STORED DATA [getappdata(0,'ah_vmfg3')& getappdata(0,'aw_vmfg3')] ARE EMPTY([])
               
                      if VMF_grid_found == 1 %IF VMF GRID FILES ARE FOUND,THEN COMPUTE VMF COE USING "SearchReadVMFgrid"
               
                         %Call the "SearchReadVMFgrid.m" function TO COMPUTE VMF3 COEFFICIENTS
                         [ah3,aw3,VMF_model] = SearchReadVMFgrid(UTCtime,latD,lonD,hgt,VMFgrids);
              
                      else %IF VMF GRID FILES ARE NOT FOUND
               
                           %ASSIGN EMPTY([]) MATRIX          
                           ah3       = [];
                           aw3       = [];
                           VMF_model = [];
                  
                      end %//if VMF_grid_found == 1
        
                 end  
         
              else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE NOT THE VMF GRID MODELS
               
                   if VMF_grid_found == 1 %IF VMF GRID FILES ARE FOUND,THEN COMPUTE VMF COE USING "SearchReadVMFgrid"
               
                      %Call the "SearchReadVMFgrid.m" function TO COMPUTE VMF3 COEFFICIENTS
                      [ah3,aw3,VMF_model] = SearchReadVMFgrid(UTCtime,latD,lonD,hgt,VMFgrids);
              
                   else %IF VMF GRID FILES ARE NOT FOUND
               
                        %ASSIGN EMPTY([]) MATRIX          
                        ah3       = [];
                        aw3       = [];
                        VMF_model = [];
                  
                   end %//if VMF_grid_found == 1

              end 
          
          %IF 1ST PREAMBLE FAILS, GET COEFFICIENTS FROM GPT3 MODELS 
          %EG:.gpt3_5.grd or gpt3_1.grd 
      
          %IF VMF GRID MODEL IS THE VMF1 MODEL OR VMF3 GRID FILES WERE NOT AVAILABLE
          if any([strncmpi(VMF_model,'VMF1',4),all([isempty(VMF_model),isempty(ah3),isempty(aw3)])])

             %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE THE SAME
             if any([strncmpi(dryModel,wetModel,14),strfind(dryModel,wetModel)]) 
          
                if any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)])
             
                   %RETRIEVE STORED VMF COEFFICIENTS
                   ah3 = getappdata(0,'ah_vmf3');
                   aw3 = getappdata(0,'aw_vmf3');
                end  
            
             else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel %subroutines ARE NOT THE SAME 
          
                  %IF DRY MODEL IS ANY OF THE GPT2w MODEL 1° OR 5° GRID & WET MODEL IS NOT
                  if all([any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),...
                          ~any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])])     
                
                     %RETRIEVE STORED VMF COEFFICIENTS
                     if any([strncmpi(getappdata(0,'Drymodel'),'GPT3',4),strfind(getappdata(0,'Drymodel'),'GPT3')])
                         
                        ah3 = getappdata(0,'ah_vmf3');
                        aw3 = getappdata(0,'aw_gpt3_vmf3');
                        
                     end 
                 
                  %IF WET MODEL IS ANY OF THE GPT3 MODEL (1° OR 5° GRID) & DRY MODEL IS NOT   
                  elseif all([~any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),...
                               any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])])     
                
                         %RETRIEVE STORED VMF3 COEFFICIENTS           
                         if any([strncmpi(getappdata(0,'Wetmodel'),'GPT3',4),strfind(getappdata(0,'Wetmodel'),'GPT3')])
                         
                            ah3 = getappdata(0,'ah_gpt3_vmf3');
                            aw3 = getappdata(0,'aw_vmf3');
                            
                         end
                 
                  else %IF DRY OR WET TROPO MODEL IS NOT ANY OF THE GPT3 MODEL (1° OR 5° GRID)
                  
                       %CHECK LENGTH OF GPT3 GRID VALUES{FOR GPT3,LENGTH SHOULD BE 14}
                        if length(grid_w) == 14
                
                           %Call the "vmfCOE_gpt3.m" function TO COMPUTE VMF3 COEFFICIENTS
                           [ah3,aw3] = vmfCOE_gpt3(UTCtime,latR,lonR,hgt,grid_w,grid_resw,Timevar_w);
           
                        end      
             
                  
                  end %//if all([any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),...
                      %         ~any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])])     
              
             end %//if any([strncmpi(drymodel,wetmodel,14),strfind(drymodel,wetmodel)]) 
         
          end %//if any([strncmpi(VMF_model,'VMF1',4),all([isempty(VMF_model),isempty(ah3),isempty(aw3)])])
        
       %IF WET MF MODEL IS RATHER ANY OF THE GPT3 MODEL(1° OR 5°) AND NOT ANY OF VMF3 MODELS
       elseif any([strncmpi(MFw_model,'GPT3 (1° x 1°)',14),strncmpi(MFw_model,'GPT3 (5° x 5°)',14)])
              
             %CHECK IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE THE SAME
             if any([strncmpi(dryModel,wetModel,14),strfind(dryModel,wetModel)]) 
          
                %IF DRY TROPO MODEL IS ANY OF THE GPT2w MODEL(1° OR 5°)
                if any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)])
             
                   %RETRIEVE STORED VMF COEFFICIENTS 
                   ah3 = getappdata(0,'ah_vmf3');
                   aw3 = getappdata(0,'aw_vmf3');  
                end   
               
             else %IF DRY & WET TROPO MODELS FROM THE TropModel or TropoModel subroutines ARE NOT THE SAME
             
                  %IF DRY MODEL IS ANY OF THE GPT3 MODEL (1° OR 5° GRID) & WET MODEL IS NOT
                  if all([any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),...
                         ~any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])])     
                
                     %RETRIEVE STORED VMF COEFFICIENTS
                     ah3 = getappdata(0,'ah_vmf3');
                     aw3 = getappdata(0,'aw_gpt3_vmf3');
                
                     %IF WET MODEL IS ANY OF THE GPT3 MODEL (1° OR 5° GRID) & DRY MODEL IS NOT    
                  elseif all([~any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]),...
                               any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])])     
                
                         %RETRIEVE STORED VMF COEFFICIENTS   
                         ah3 = getappdata(0,'ah_gpt3_vmf3');
                         aw3 = getappdata(0,'aw_vmf3');
                    
                  else %IF DRY OR WET TROPO MODEL IS NOT ANY OF THE GPT3 MODEL (1° OR 5° GRID)  
                       %Call the "vmfCOE_gpt3.m" function TO COMPUTE VMF3 COEFFICIENTS
                      [ah3,aw3] = vmfCOE_gpt3(UTCtime,latR,lonR,hgt,grid_w,grid_resw,Timevar_w);
                 
                  end    
              
             end %//if any([strncmpi(drymodel,wetmodel,14),strfind(drymodel,wetmodel)])   
 
       
           end   %//if any([strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12),strncmpi(MFw_model,'VMF Z',5)]) 
         
         
    end %//if all([any([strncmpi(MFh_model,'VMF',3),strfind(MFh_model,'VMF'),strncmpi(MFh_model,'GPT3',4),strfind(MFh_model,'GPT3')]),...
        %          any([~strncmpi(MFw_model,'VMF',3),~strfind(MFw_model,'VMF'),~strncmpi(Wetmodel,'GPT3',4),~strfind(Wetmodel,'GPT3')])])
            
            
end %//if any([strncmpi(MFh_model,MFw_model,12),strfind(MFh_model,MFw_model),any([strncmpi(MFh_model,'VMF ZHD',7),...
    %          strncmpi(MFw_model,'VMF ZWD',7)])])
        
%******************************************END OF VMF3coe.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  
%C.5.SUBROUTINE TO COMPUTE VMF COEFFICIENTS FROM GPT2w 1° & 5° GRID MODEL
function [ah,aw] = vmfCOE_gpt2w(UTCtime,lat,lon,hell,grid,grid_res,Timevar)

%**************************************************************************
%DESCRIPTION:
%           This subroutine determines the hydrostatic and wet mapping ... * 
%           function coefficients ah and aw, for specific sites near the   * 
%           Earth surface. It is based on a 1° or 5°  external grid file   *                                                *
%          ('gpt2_5w.grd' or('gpt2_1w.grd') or('gpt2_5.grd')               *

%USAGE:                                                                    *
%      [ah,aw] = vmfCOE_gpt2w(UTCtime,lat,lon,hell,grid,grid_res,Timevar)  *

%INPUTs:                                                                   *
%1.   UTCtime:  UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *                                                                
%2.       lat:  Station geodetic latitude in radians [-pi/2:+pi/2] (n x 1) *
%3.       lon:  Station geodetic longitude in radians [-pi:pi] or [0:2pi]  *
%4.      hell:  ellipsoidal height in [meters] (n x 1)                     *
%5.      grid: grid values in cells extracted from the grid file           *
%    grid_res: grid resolution (°) (possible: 1 or 5)
%6    Timevar:  case 1: no time variation but static quantities            *
%               case 0: with time variation (annual and semiannual terms)  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% OUTPUTs:                                                                 *
%1.      ah: hydrostatic mapping function coefficient at zero height (VMF1)* 
%2.      aw: wet mapping function coefficient(VMF1)(vector of length nstat)*

%Original codes by Böhm et al 2014                                         *
%Modified by: Samuel Osah,Msc Geomatic Engineering ,2016                   *    
%             Email: osahsamuel@yahoo.ca                                   *
%             Tel:+233 (0)246137410/+233 (0)509438484                      *        
%==========================================================================+
%**************************************************************************+
%**************************************************************************+

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
      ah=zeros(nrow_time,1);
      aw=deal(ah);%create copy of ah in aw
   
      for i=1:nrow_time
       
          if grid_res == 1 %if grid resolution is 1°
              
             %Call the "GPT2w_1_x_1.m" Function
             [ah(i,1),aw(i,1)] = GPT2w_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid,Timevar); 
          
          elseif grid_res == 5 %if grid resolution is 5°
              
                 %Call the "GPT2w_5_x_5.m" Function
                [ah(i,1),aw(i,1)] = GPT2w_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid,Timevar);
        
          end 
      end
      
   else 
       %*****INITIALIZE OUTPUTs 
      ah=zeros(nrow_pos,nrow_time);
      aw=deal(ah);%create copy of ah in aw
    
     for i=1:nrow_time %LOOP OVER TIME
        for j=1:nrow_pos %LOOP OVER POSITIONS
            
            %****CHECK GRID RESOLUTION
            
            if grid_res == 1 %if grid resolution is 1°
                
               %Call the the 1 degree grid version; "GPT2w_1_x_1.m"  
               [ah(j,i),aw(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);
           
            elseif grid_res == 5 %if grid resolution is 5°
                
                   %Call the the 5 degree grid version; "GPT2w_5_x_5.m"  
                   [ah(j,i),aw(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);                 
            end

        end
     end
    
   end %//if isequal(Identical,1)
    
else
    %*****INITIALIZE OUTPUTs 
    ah=zeros(nrow_pos,nrow_time);
    aw=deal(ah);%create copy of ah in aw
    
    for i=1:nrow_time %LOOP OVER TIME
        for j=1:nrow_pos %LOOP OVER POSITIONS
            
            %****CHECK GRID RESOLUTION
            
            if grid_res == 1 %if grid resolution is 1° 
                
               %Call the the 1 degree grid version; "GPT2w_1_x_1.m"  
               [ah(j,i),aw(j,i)] = GPT2w_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);
           
            elseif grid_res == 5 %if grid resolution is 5° 
                
                   %Call the the 5 degree grid version; "GPT2w_5_x_5.m"  
                   [ah(j,i),aw(j,i)] = GPT2w_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);
                   
            end

           
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF vmfCOE_gtp2w.m  %%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%C.6.SUBROUTINE TO COMPUTE VMF COEFFICIENTS FROM GPT2w 1° & 5° GRID MODEL
function [ah,aw] = vmfCOE_gpt2(UTCtime,lat,lon,hell,grid,Timevar)

%**************************************************************************
%DESCRIPTION:
%           This subroutine determines the hydrostatic and wet mapping ... * 
%           function coefficients ah and aw, for specific sites near the   * 
%           Earth surface. It is based on a 1° or 5°  external grid file   *                                                *
%          ('gpt2_5w.grd' or('gpt2_1w.grd') or('gpt2_5.grd')               *

%USAGE:                                                                    *
%      [ah,aw] = vmfCOE_gpt2w(UTCtime,lat,lon,hell,grid,grid_res,Timevar)  *

%INPUTs:                                                                   *
%1.   UTCtime:  UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *                                                                
%2.       lat:  Station geodetic latitude in radians [-pi/2:+pi/2] (n x 1) *
%3.       lon:  Station geodetic longitude in radians [-pi:pi] or [0:2pi]  *
%4.      hell:  ellipsoidal height in [meters] (n x 1)                     *
%5.      grid: grid values in cells extracted from the grid file           *
%6    Timevar:  case 1: no time variation but static quantities            *
%               case 0: with time variation (annual and semiannual terms)  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% OUTPUTs:                                                                 *
%1.      ah: hydrostatic mapping function coefficient at zero height (VMF1)* 
%2.      aw: wet mapping function coefficient(VMF1)(vector of length nstat)*

%Original codes by Böhm et al 2014                                         *
%Modified by: Samuel Osah,Msc Geomatic Engineering ,2016                   *    
%             Email: osahsamuel@yahoo.ca                                   *
%             Tel:+233 (0)246137410/+233 (0)509438484                      *        
%==========================================================================+
%**************************************************************************+
%**************************************************************************+

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
      ah=zeros(nrow_time,1);
      aw=deal(ah);%create copy of ah in aw
   
      for i=1:nrow_time
       
          %Call the "GPT2_5_x_5.m" Function
          [ah(i,1),aw(i,1)] = GPT2_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid,Timevar);
      end
                                        
   else 
       %*****INITIALIZE OUTPUTs 
      ah=zeros(nrow_pos,nrow_time);
      aw=deal(ah);%create copy of ah in aw
    
     for i=1:nrow_time %LOOP OVER TIME
         
        for j=1:nrow_pos %LOOP OVER POSITIONS
            
            %Call the the 5 degree grid version; "GPT2_5_x_5.m"  
            [ah(j,i),aw(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);                 
        end 
     end 
    
   end  %//if isequal(Identical,1)
    
else
    %*****INITIALIZE OUTPUTs 
    ah=zeros(nrow_pos,nrow_time);
    aw=deal(ah);%create copy of ah in aw
    
    for i=1:nrow_time %LOOP OVER TIME
        
        for j=1:nrow_pos %LOOP OVER POSITIONS
            
            %Call the the 5 degree grid version; "GPT2_5_x_5.m"  
            [ah(j,i),aw(j,i)] = GPT2_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);        
        end   
    end 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF vmfCOE_gtp2.m  %%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%C.7.SUBROUTINE TO COMPUTE VMF COEFFICIENTS FROM GPT2w 5° GRID MODEL
function [ah,aw] = GPT2w_5_x_5(UTCtime,lat,lon,hell,grid,Timevar)

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
%      [ah,aw] = GPT2w_5(UTCtime,lat,lon,hell,Timevar) *
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
%1.      ah: hydrostatic mapping function coefficient at zero height (VMF1)* 
%2.      aw: wet mapping function coefficient (VMF1) (vector of length nstat)
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


%C.8.SUBROUTINE TO COMPUTE VMF COEFFICIENTS FROM GPT2w 1° GRID MODEL
function [ah,aw] = GPT2w_1_x_1(UTCtime,lat,lon,hell,grid,Timevar)

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
%      [ah,aw] = GPT2w_1_x_1(UTCtime,lat,lon,hell,Timevar)                 *
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
%1.      ah: hydrostatic mapping function coefficient at zero height (VMF1)* 
%2.      aw: wet mapping function coefficient (VMF1) (vector of length nstat)
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


%C.9.**SUBROUTINE TO COMPUTE VMF COEFFICIENTS FROM GPT3 1° & 5° GRID MODEL
function [ah,aw] = vmfCOE_gpt3(UTCtime,lat,lon,hell,grid,grid_res,Timevar)

%**************************************************************************
%DESCRIPTION:
%           This subroutine determines the hydrostatic and wet mapping ... * 
%           function coefficients ah and aw, for specific sites near the   * 
%           Earth surface. It is based on a 1° or 5°  external grid file   *                                                *
%           ('gpt3_5.grd' or('gpt3_1.grd')                                  *

%USAGE:                                                                    *
%      [ah,aw] = vmfCOE_gpt3(UTCtime,lat,lon,hell,grid,grid_res,Timevar)  *

%INPUTs:                                                                   *
%1.   UTCtime:  UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *                                                                
%2.       lat:  Station geodetic latitude in radians [-pi/2:+pi/2] (n x 1) *
%3.       lon:  Station geodetic longitude in radians [-pi:pi] or [0:2pi]  *
%4.      hell:  ellipsoidal height in [meters] (n x 1)                     *
%5.      grid: grid values in cells extracted from the grid file           *
%    grid_res: grid resolution (°) (possible: 1 or 5)
%6    Timevar:  case 1: no time variation but static quantities            *
%               case 0: with time variation (annual and semiannual terms)  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% OUTPUTs:                                                                 *
%1.      ah: hydrostatic mapping function coefficient at zero height (VMF1)* 
%2.      aw: wet mapping function coefficient(VMF1)(vector of length nstat)*

%Original codes by Böhm et al 2014                                         *
%Modified by: Samuel Osah,Msc Geomatic Engineering ,2016                   *    
%             Email: osahsamuel@yahoo.ca                                   *
%             Tel:+233 (0)246137410/+233 (0)509438484                      *        
%==========================================================================+
%**************************************************************************+
%**************************************************************************+

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
      ah=zeros(nrow_time,1);
      aw=deal(ah);%create copy of ah in aw
   
      for i=1:nrow_time
       
          if grid_res == 1 %if grid resolution is 1°
              
             %Call the "GPT2w_1_x_1.m" Function
             [ah(i,1),aw(i,1)] = GPT3_1_x_1(UTCtime(i,:),lat(i),lon(i),hell(i),grid,Timevar); 
          
          elseif grid_res == 5 %if grid resolution is 5°
              
                 %Call the "GPT2w_5_x_5.m" Function
                [ah(i,1),aw(i,1)] = GPT3_5_x_5(UTCtime(i,:),lat(i),lon(i),hell(i),grid,Timevar);
        
          end 
      end
      
   else 
       %*****INITIALIZE OUTPUTs 
      ah=zeros(nrow_pos,nrow_time);
      aw=deal(ah);%create copy of ah in aw
    
     for i=1:nrow_time %LOOP OVER TIME
        for j=1:nrow_pos %LOOP OVER POSITIONS
            
            %****CHECK GRID RESOLUTION
            
            if grid_res == 1 %if grid resolution is 1°
                
               %Call the the 1 degree grid version; "GPT2w_1_x_1.m"  
               [ah(j,i),aw(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);
           
            elseif grid_res == 5 %if grid resolution is 5°
                
                   %Call the the 5 degree grid version; "GPT2w_5_x_5.m"  
                   [ah(j,i),aw(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);
                   
            end

           
        end
     end
    
   end %//if isequal(Identical,1)
    
else
    %*****INITIALIZE OUTPUTs 
    ah=zeros(nrow_pos,nrow_time);
    aw=deal(ah);%create copy of ah in aw
    
    for i=1:nrow_time %LOOP OVER TIME
        for j=1:nrow_pos %LOOP OVER POSITIONS
            
            %****CHECK GRID RESOLUTION
            
            if grid_res == 1 %if grid resolution is 1° 
                
               %Call the the 1 degree grid version; "GPT2w_1_x_1.m"  
               [ah(j,i),aw(j,i)] = GPT3_1_x_1(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);
           
            elseif grid_res == 5 %if grid resolution is 5° 
                
                   %Call the the 5 degree grid version; "GPT2w_5_x_5.m"  
                   [ah(j,i),aw(j,i)] = GPT3_5_x_5(UTCtime(i,:),lat(j),lon(j),hell(j),grid,Timevar);  
            end
 
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF vmfCOE_gtp3.m  %%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%C.10.SUBROUTINE TO COMPUTE VMF COEFFICIENTS FROM GPT3 5° GRID MODEL
function [ah,aw] = GPT3_5_x_5(UTCtime,lat,lon,h_ell,grid,Timevar)

%**************************************************************************
%DESCRIPTION:
%            This subroutine determines pressure, temperature, temperature *
%            lapse rate, mean temperature of the water vapor, water vapor  * 
%            pressure, hydrostatic and wet mapping function coefficients   * 
%            ah and aw, water vapour decrease factor, geoid undulation, and* 
%            empirical tropospheric gradients for specific sites near the  *
%            Earth surface. It is based on a  5 x 5 degree external grid   *
%            file ('gpt3_5.grd') with mean  values as well as sine and     * 
%            cosine amplitudes for the annual and semiannual variation of  *
%            the coefficients.
%USAGE:
%      [ah,aw] = GPT3_5_x_5(UTCtime,lat,lon,hell,grid,Timevar)             *
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
%1.      ah: hydrostatic mapping function coefficient at zero height (VMF3)* 
%2.      aw: wet mapping function coefficient (VMF3) (vector of length nstat)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%NOTE:                                                                     *     
%    The hydrostatic mapping function coefficients have to be used with the*
%     height dependent Vienna Mapping Function 3  because the              *
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GPT3_5 x 5.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%C.11.******SUBROUTINE TO COMPUTE VMF COEFFICIENTS FROM GPT3 1° GRID MODEL
function [ah,aw] = GPT3_1_x_1(UTCtime,lat,lon,h_ell,grid,Timevar)

%**************************************************************************
%DESCRIPTION:
%            This subroutine determines pressure, temperature, temperature *
%            lapse rate, mean temperature of the water vapor, water vapor  * 
%            pressure, hydrostatic and wet mapping function coefficients   * 
%            ah and aw, water vapour decrease factor, geoid undulation, and* 
%            empirical tropospheric gradients for specific sites near the  *
%            Earth surface. It is based on a  1 x 1 degree external grid   *
%            file ('gpt3_1.grd') with mean  values as well as sine and     * 
%            cosine amplitudes for the annual and semiannual variation of  *
%            the coefficients.
%USAGE:
%      [ah,aw] = GPT3_1_x_1(UTCtime,lat,lon,hell,grid,Timevar) 
%INPUTs:                                                                   *
%1.   UTCtime:  UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *                                                                
%2.       lat:  Station geodetic latitude in radians [-pi/2:+pi/2] (n x 1) *
%3.       lon:  Station geodetic longitude in radians [-pi:pi] or [0:2pi]  *
%4.      hell:  ellipsoidal height in [meters] (n x 1)                     *
%5.      grid: grid values in cells extracted from the grid file      *
%6.   Timevar:  case 1: no time variation but static quantities            *
%               case 0: with time variation (annual and semiannual terms)  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% OUTPUTs:
%1.      ah: hydrostatic mapping function coefficient at zero height (VMF3)* 
%2.      aw: wet mapping function coefficient (VMF3) (vector of length nstat)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%NOTE:                                                                     *     
%    The hydrostatic mapping function coefficients have to be used with the*
%     height dependent Vienna Mapping Function 3  because the              *
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GPT3_1_x_1.m  %%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%C.12.******SUBROUTINE TO COMPUTE VMF COEFFICIENTS FROM GPT2 5° GRID MODEL
function [ah,aw] = GPT2_5_x_5(UTCtime,lat,lon,hell,grid,Timevar)

%**************************************************************************
%DESCRIPTION:                                                              *
%            This subroutine determines pressure, temperature, temperature *
%            lapse rate,water vapor pressure, hydrostatic and wet mapping  * 
%            function coefficients ah and aw, water vapour decrease        * 
%            factor and geoid undulation for specific sites near the Earth * 
%            surface. It is based on a 5 x 5 degree external grid file     *
%            ('gpt2_5.grd') with mean  values as well as sine and cosine   *
%             amplitudes for the annual  and semiannual variation of the   * 
%             coefficients.                                                *
%USAGE:                                                                    *
%      [ah,aw] = GPT2(UTCtime,lat,lon,hell,grid,Timevar)                   *
%INPUTs:                                                                   *
%1.   UTCtime:  UTC time in [Year,Month,Day,Hour,Minute,Seconds]           *                                                                
%2.       lat:  Station geodetic latitude in radians [-pi/2:+pi/2] (n x 1) *
%3.       lon:  Station geodetic longitude in radians [-pi:pi] or [0:2pi]  *
%4.      hell:  ellipsoidal height in [meters] (n x 1)                     *
%5.      grid:  grid values in cells extracted from the grid file          *
%6    Timevar:  case 1: no time variation but static quantities            *
%               case 0: with time variation (annual and semiannual terms)  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% OUTPUTs:
%1.      ah: hydrostatic mapping function coefficient at zero height (VMF1)* 
%2.      aw: wet mapping function coefficient (VMF1) (vector of length nstat)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GPT2.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%************END OF VIENNA MAPPING FUNCTIONS(VMF1 / VMF3)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



%A_B_C.SUBROUTINE TO COMPUTE JULIAN,MODIFIED DAY(JD,MJD) &  DAY OF YEAR(DoY)
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
 
%A_B_C_1.*********SUBROUTINE TO COMPUTE DAY OF YEAR(DoY)
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

%A_B_C_2*********SUBROUTINE TO CONVERT TWO DIGITs TO YEAR 4 DIGITS
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
%D.SUBROUTINE TO COMPUTE BLACK & EISNER MAPPING FUNCTION
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

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%E.SUBROUTINE TO COMPUTE BLACK MAPPING FUNCTION
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [MFh,MFw]=Black_MF(satELEV,Temperature)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            "Black_MF" is a subroutine that Computes the wet and dry      *
%             Tropospheric mapping function Using Black mapping function   *
%USAGE:                                                                    *                                                           *
%      [MFh,MFw]=Black_MF(satELEV,Temperature)                             *                                                                   *
%INPUT:                                                                    *
%      satELEV : Satellite Elevation in [degrees]                          *
%  Temperature : Surface temperature in [kelvin]                           *
%OUTPUT:                                                                   *
%1.     MFh ==> Hydrostatic/dry mapping function                           *
%2.     MFw ==> Non Hydrostatic/wet mapping function                       *
%==========================================================================
%REFERENCE:                                                                *
%1.        Black, H. D. (1978). An easily implemented algorithm for the    *
%                               tropospheric range correction.             * 
%3.        Seeber, G. (2003). Satellite Geodesy: Foundations, Methods &    *
%                             Applications. 2nd completely revised  Ed     *
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************

%*****Assignment
T = Temperature;

%1.*****INITIALIZING OUTPUT VARIABLEs
satEL=zeros(size(satELEV,1),size(satELEV,2));%Assign zeros of nxm to SHD
[MFh,MFw,Lc]=deal(satEL);%copy the contents of satEl to all the requested outputs
%2.**************COMPUTE MAPPING FUNCTIONs(MFh & MFw)
Re=6378137;%Distance from the center of Earth to the observation station

%COMPUTE MEAN DRY TROPOSPHERE HEIGHT(Hd)
Hd = 148.98.*(T-4.12);%[meters]

%COMPUTE MEAN WET TROPOSPHERE HEIGHT(Hw)
Hw = 13000;%[meters]

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
      %******COMPUTE SCALE FACTOR (Lc)                                                               
      try
         Lc(j,i)=0.167-(0.076+0.00015.*(T(i)-273.16)).*exp(-0.3.*(satEL(j,i).*180/pi));         
      catch
           Lc=0.85;%provision of constant Lc value      
      end                    
          if Lc(j,i)==1               
             MFh(j,i)=csc(satEL(j,i));
             MFw(j,i)=csc(satEL(j,i));             
          else
              MFh(j,i)=(1-((cos(satEL(j,i)))./(1+(1-Lc(j,i)).*(Hd(i)./Re))).^2).^(-0.5);%Hydrostatic Mapping Function
              MFw(j,i)=(1-((cos(satEL(j,i)))./(1+(1-Lc(j,i)).*(Hw./Re))).^2).^(-0.5);%Wet Mapping Function
          
          end %\\if Lc(j,i)==1

   end   %for j=1:nrow

end  %for i=1:ncol

%%%%%%%%%%%%%%%%%%%%%%%%END OF MF_Black.m %%%%%%%%%%%%%%%%%%%%%%% 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%F.SUBROUTINE TO COMPUTE Ifadis MAPPING FUNCTION
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [MFh,MFw]=Ifadis_MF(satELEV,Temperature,Pressure,WVP)
                                                               
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            "Ifadis_MF" is a subroutine that Computes the wet and dry     *
%             Tropospheric mapping function Using Ifadi mapping function   *

%USAGE:                                                                    *
%      General:                                                            *
%              [MFh,MFw]=Ifadis_MF(satELEV,Temperature,Pressure,WVP)       *

%INPUT:                                                                    *
%      satELEV : Satellite Elevation in [degrees]                          *
%  Temperature : Surface temperature in [kelvin]                           *
%     Pressure : Atmospheric Pressure in millibars(mbar /hPa)              *
%          WVP : Water Vapor Pressure(e) in millibars(mbar /hPa)           *
%OUTPUT:                                                                   *
%1.     MFh ==> Hydrostatic/dry mapping function                           *
%2.     MFw ==> Non Hydrostatic/wet mapping function                       *
%==========================================================================
%REFERENCE:                                                                *
%1.         Liu, Z., & Chen, B. (2015). A Comprehensive Evaluation and 
%           Analysis of the Performance of Multiple Tropospheric Models in 
%           China Region. IEEE Transactions on Geoscience and Remote Sensing, 
%           2015(2015), pp.1-17.
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%*****Assignment
T  = Temperature;
P  = Pressure;
es = WVP;

%1.*****INITIALIZING OUTPUT VARIABLEs
satEL=zeros(size(satELEV,1),size(satELEV,2));%Assign zeros of nxm to SHD
[MFh,MFw]=deal(satEL);%copy the contents of satEl to all the requested outputs

%1.1*****INITIALIZE MAPPING COEFFICIENTs
ah=zeros(size(satELEV,2),1);%Find number of stations 
[bh,ch,aw,bw,cw,Ts]=deal(ah);%copy the contents of ah to other coefficients

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
       
       %COMPUTE DRY MAPPING COEFFICIENTs(ah,bh,ch)
       Ts(i)=T(i)-273.16; %Convert temperture kelvin(K) to degree Celcius(C) 
       ah(i,1)=0.1237e-2+(0.1316e-6.*(P(i)-1000))+(0.1378e-5.*(Ts(i)-15))+(0.8057e-5.*sqrt(es(i)));
       bh(i,1)=0.3333e-2+(0.1946e-6.*(P(i)-1000))+(0.1040e-6.*(Ts(i)-15))+(0.1747e-6.*sqrt(es(i)));
       ch(i,1)=0.078;
       
       %COMPUTE WET MAPPING COEFFICIENTs(aw,bw,cw)
       aw(i,1)=0.5236e-2+(0.2471e-6.*(P(i)-1000))-(0.1724e-6.*(Ts(i)-15))+(0.1328e-4.*sqrt(es(i)));
       bw(i,1)=0.1705e-2+(0.7384e-6.*(P(i)-1000))+(0.3767e-6.*(Ts(i)-15))+(0.2147e-4.*sqrt(es(i)));
       cw(i,1)=0.05917;
       
      %NOW THE MAPPING FUNCTIONs(MF)
      MFh(j,i)=1./(sin(satEL(j,i))+(ah(i)./(sin(satEL(j,i))+(bh(i)./(sin(satEL(j,i))+ch(i))))));%Hydrostatic Mapping Function
      MFw(j,i)=1./(sin(satEL(j,i))+(aw(i)./(sin(satEL(j,i))+(bw(i)./(sin(satEL(j,i))+cw(i))))));%Wet Mapping Function 

   end   %for j=1:nrow

end  %for i=1:ncol
%%%%%%%%%%%%%%%%%%%%%%%%END OF Ifadis_MF.m %%%%%%%%%%%%%%%%%%%%%%% 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%G.SUBROUTINE TO COMPUTE Herring MAPPING FUNCTION
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [MFh,MFw]=Herring_MF(lat,h,satELEV,Temperature)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            "Herring_MF" is a subroutine that Computes the wet and dry    *
%            Tropospheric mapping function Using Herring mapping function  *
%USAGE:                                                                    *                                                          
%      [MFh,MFw]=Herring_MF(satELEV,Temperature)                           *                                                                   
%INPUT:                                                                    *
%          lat : station geodetic latitude in [radian]                     *
%            h : station height in [meters]
%      satELEV : Satellite Elevation in [degrees]                          *
%  Temperature : Surface temperature in [kelvin]                           *
%OUTPUT:                                                                   *
%1.     MFh ==> Hydrostatic/dry mapping function                           *
%2.     MFw ==> Non Hydrostatic/wet mapping function                       *
%==========================================================================
%REFERENCE:                                                                *
%1.        Liu, Z., & Chen, B. (2015). A Comprehensive Evaluation and      *
%          Analysis of the Performance of Multiple Tropospheric Models in  *
%          China Region. IEEE Transactions on Geoscience and Remote Sensing, 
%          2015(2015), pp.1-17.                                            *
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************
%*****Assignment
T = Temperature;

%1.*****INITIALIZING OUTPUT VARIABLEs
satEL=zeros(size(satELEV,1),size(satELEV,2));%Assign zeros of nxm to SHD
[MFh,MFw]=deal(satEL);%copy the contents of satEl to all the requested outputs

[denoFrac_dry, denoFrac_wet]=deal(satEL);%copy the contents of satEL 

%1.1*****INITIALIZE MAPPING COEFFICIENTs
ah=zeros(size(satELEV,2),1);%Find number of stations 
[bh,ch,aw,bw,cw,H,Ts]=deal(ah);%copy the contents of ah to other coefficients
[numFrac_dry, numFrac_wet]=deal(ah);%copy the contents of ah 

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
      
       %COMPUTE DRY MAPPING COEFFICIENTs(ah,bh,ch)
       Ts(i)=T(i)-273.16; %Convert temperture kelvin(K) to degree Celcius(C)
       H(i)=h(i)./1000; %Heights converted to kilometer[km]
       
       ah(i,1)=(1.2320+(0.0139.*cos(lat(i)))-(0.0209.*H(i))+(0.00215.*(Ts(i)-10))).*1.0e-3;
       bh(i,1)=(3.1612-(0.1600.*cos(lat(i)))-(0.0331.*H(i))+(0.00206.*(Ts(i)-10))).*1.0e-3;
       ch(i,1)=(171.244-(4.293.*cos(lat(i)))-(0.149.*H(i))-(0.0021.*(Ts(i)-10))).*1.0e-3;
       
       %COMPUTE WET MAPPING COEFFICIENTs(aw,bw,cw)
       aw(i,1)=(0.583-(0.011.*cos(lat(i)))-(0.052.*H(i))+(0.0014.*(Ts(i)-10))).*1.0e-3;
       bw(i,1)=(1.402-(0.102.*cos(lat(i)))-(0.101.*H(i))+(0.0020.*(Ts(i)-10))).*1.0e-3;
       cw(i,1)=(45.85-(1.91.*cos(lat(i)))-(1.29.*H(i))+(0.015.*(Ts(i)-10))).*1.0e-3;
       
      %NOW THE MAPPING FUNCTIONs(MF)
      %1.******HYDROSTATIC
      numFrac_dry(i,1)=(1+(ah(i,1)./(1+(bh(i,1)./(1+ch(i,1))))));
     denoFrac_dry(i,j)=(sin(satEL(j,i))+(ah(i,1)./(sin(satEL(j,i))+(bh(i,1)./(sin(satEL(j,i))+ch(i,1))))));
      
     MFh(j,i) = numFrac_dry(i,1)./denoFrac_dry(i,j); %Hydrostatic Mapping Function
      
      %2.******WET
       numFrac_wet(i,1)=(1+(aw(i,1)./(1+(bw(i,1)./(1+cw(i,1))))));
      denoFrac_wet(i,j)=(sin(satEL(j,i))+(aw(i,1)./(sin(satEL(j,i))+(bw(i,1)./(sin(satEL(j,i))+cw(i,1))))));
      
      MFw(j,i) = numFrac_wet(i,1)./denoFrac_wet(i,j) ;%Wet Mapping Function 

   end   %for j=1:nrow

end  %for i=1:ncol

%%%%%%%%%%%%%%%%%%%%%%%%END OF MF_Herring.m %%%%%%%%%%%%%%%%%%%%%%% 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%H.SUBROUTINE TO COMPUTE Chao MAPPING FUNCTION
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [MFh,MFw]=Chao_MF(satELEV)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            "Chao_MF" is a subroutine that Computes the wet and dry       *
%             Tropospheric mapping function Using Chao mapping function    *

%USAGE:                                                                    *                                                           *
%      [MFh,MFw]=Chao_MF(satELEV)                                          *
%INPUT:                                                                    *
%      satELEV : Satellite Elevation in [degrees]                          *

%OUTPUT:                                                                   *
%1.     MFh ==> Hydrostatic/dry mapping function                           *
%2.     MFw ==> Non Hydrostatic/wet mapping function                       *
%==========================================================================
%REFERENCE:                                                                *
%1.        Chao, C. C. (1971) A new method to predict wet zenith range     *
%          correction from surface measurements. Jet Propulsion Laboratory,* 
%          California Inst.  Technol., Pasadena, CA, USA, Technical Report * 
%          32-1602, pp. 33–37.                                             * 

%2.         Liu, Z., & Chen, B. (2015). A Comprehensive Evaluation and     *
%           Analysis of the Performance of Multiple Tropospheric Models in * 
%           China Region. IEEE Transactions on Geoscience and Remote Sensing, 
%           2015(2015), pp.1-17.                                           *
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************
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
       MFh(j,i)=1./(sin(satEL(j,i))+(0.00143./(tan(satEL(j,i))+0.0445)));%Hydrostatic Mapping Function
       MFw(j,i)=1./(sin(satEL(j,i))+(0.00035./(tan(satEL(j,i))+0.017)));%Wet Mapping Function         
 
   end   %for j=1:nrow

end  %for i=1:ncol

%%%%%%%%%%%%%%%%%%%%%%%%END OF Chao_MF.m %%%%%%%%%%%%%%%%%%%%%%% 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%I.SUBROUTINE TO COMPUTE Moffett MAPPING FUNCTION
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [MFh,MFw]=Moffett_MF(satELEV)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            "Moffett_MF" is a subroutine that Computes the wet and dry    *
%             Tropospheric mapping function Using Moffett Model            *

%USAGE:                                                                    *                                                           *
%      [MFh,MFw]=Moffett_MF(satELEV)                                       *
%INPUT:                                                                    *
%      satELEV : Satellite Elevation in [degrees]                          *

%OUTPUT:                                                                   *
%1.     MFh ==> Hydrostatic/dry mapping function                           *
%2.     MFw ==> Non Hydrostatic/wet mapping function                       *
%==========================================================================
%REFERENCE:                                                                *
%1.       Mendes, V.B. (1999). Modeling the neutral-atmosphere propagation * 
%         delay in radiometric space techniques, Ph.D. Thesis, Department  *
%         of Geodesy and Geomatics Engineering, University of New Brunswick, 
%         Fredericton, New Brunswick, Canada, Technical Report No. 199.    *
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************
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
           el=satEL(j,i).* 180/pi; %convert satEL to degrees
           dry_elev = sqrt(el.^2 + 6.25) * pi / 180 ;  % convert it back to rad
           wet_elev = sqrt(el.^2 + 2.25) * pi / 180  ; % convert it back to rad
           MFh(j,i)=1./sin(dry_elev);%Hydrostatic Mapping Function
           MFw(j,i)=1./wet_elev;%Wet Mapping Function 
       catch
            MFh(j,i)=1./sin(sqrt(satEL(j,i).^2+(6.25* pi / 180)));%Hydrostatic Mapping Function
            MFw(j,i)=1./sin(sqrt(satEL(j,i).^2+(2.25* pi / 180)));%Wet Mapping Function
            
       end
 
   end   %for j=1:nrow

end  %for i=1:ncol

%%%%%%%%%%%%%%%%%%%%%%%%END OF Moffett_MF.m %%%%%%%%%%%%%%%%%%%%%%% 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%J.SUBROUTINE TO COMPUTE Hopfield MAPPING FUNCTION
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [MFh,MFw]=Hopfield_MF(satELEV)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            "Hopfield_MF" is a subroutine that Computes the wet and dry   *
%             Tropospheric mapping function Using Hopfield Model           *

%USAGE:                                                                    *                                                           *
%      [MFh,MFw]=Hopfield_MF(satELEV)                                      *
%INPUT:                                                                    *
%      satELEV : Satellite Elevation in [degrees]                          *

%OUTPUT:                                                                   *
%1.     MFh ==> Hydrostatic/dry mapping function                           *
%2.     MFw ==> Non Hydrostatic/wet mapping function                       *
%==========================================================================
%REFERENCE:                                                                *
%1.--------"GPS Theory and application",edited by B.Parkinson,J.Spilker,   *
%                                                         P.Enge, AIAA,1996*
%2.--------GPS Theory and Practice,Hofmann-Wellenhof, Lichtenegger, and... *
%                            Collins,Fifth, revised edition  pages 106-115.* 
%3.--------"Global Positioning System, Mishra & Enge", pg 172              *
%4.-------- Modeling of Tropospheric Delays Using ANFIS,Wayan Suparta;...  *
%                                                Kemal Maulana Alhasa(2016)*
%5.--------GPS Theory,Algorithms and Applications 2ed,Guochang Xu,June 2007*
%                                                    pg 66 eqns(5.95-5.108)*
%3.--------Seeber, G. (2003). Satellite Geodesy: Foundations, Methods &    *
%                             Applications. 2nd completely revised  Ed
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************
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
           el=satEL(j,i).* 180/pi; %convert satEL to degrees
           dry_elev = sqrt(el.^2 + 6.25) * pi / 180 ;  % convert it back to rad
           wet_elev = sqrt(el.^2 + 2.25) * pi / 180  ; % convert it back to rad
           MFh(j,i)=1./sin(dry_elev);%Hydrostatic Mapping Function
           MFw(j,i)=1./wet_elev;%Wet Mapping Function 
       catch
            MFh(j,i)=1./sin(sqrt(satEL(j,i).^2+(6.25* pi / 180)));%Hydrostatic Mapping Function
            MFw(j,i)=1./sin(sqrt(satEL(j,i).^2+(2.25* pi / 180)));%Wet Mapping Function
            
       end        
 
   end   %for j=1:nrow

end  %for i=1:ncol

%%%%%%%%%%%%%%%%%%%%%%%%END OF Hopfield_MF.m %%%%%%%%%%%%%%%%%%%%%%% 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%K.SUBROUTINE TO COMPUTE Davis MAPPING FUNCTION
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function MFh=Davis_MF(satELEV,Temperature,Pressure,WVP)
                                                               
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            "Davis_MF" is a subroutine that Computes the  dry             *
%             Tropospheric mapping function Using Davis mapping function   *

%USAGE:                                                                    *
%      General:                                                            *
%              MFh=Davis_MF(satELEV,Temperature,Pressure,WVP)              *

%INPUT:                                                                    *
%      satELEV : Satellite Elevation in [degrees]                          *
%  Temperature : Surface temperature in [kelvin]                           *
%     Pressure : Atmospheric Pressure in millibars(mbar /hPa)              *
%          WVP : Water Vapor Pressure(e) in millibars(mbar /hPa)           *
%OUTPUT:                                                                   *
%1.     MFh ==> Hydrostatic/dry mapping function                           *
%==========================================================================
%REFERENCE:                                                                *
%1.        Bai, Z. (2005). Near-real-time GPS sensing of atmospheric water * 
%          vapour,PhD Thesis,Queensland University of Technology.          *

%2.       Schüler, T. (2001). On ground-based GPS tropospheric delay       *
%         estimation. PhD Thesis Univ. der Bundeswehr München.             *
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%*****Assignment
T=Temperature;
P=Pressure;
es=WVP;

%1.*****INITIALIZING OUTPUT VARIABLEs
satEL=zeros(size(satELEV,1),size(satELEV,2));%Assign zeros of nxm to SHD
MFh=deal(satEL);%copy the contents of satEl to MFh

%1.1*****INITIALIZE MAPPING COEFFICIENTs
ah=zeros(size(satELEV,2),1);%Find number of stations 
[bh,ch,Ts]=deal(ah);%copy the contents of ah to other coefficients

beta=zeros(size(satELEV,2),1);
Ht=deal(beta);%copy the contents of beta to Ht

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
      
      Ts(i)=T(i)-273.16; %Convert temperture kelvin(K) to degree Celcius(C) 
      
      %COMPUTE TEMPERATURE LAPSE RATE(beta)
     try
        beta(i,1)=-5.930-0.0359.*Ts(i); %[C/km]
     catch
           beta(i,1)=-6.5; %Default [C/km]
      end
      %COMPUTE HEIGHT OF TROPOPAUSE(Ht)
      try
         Ht(i,1)=7.508+(2.421.*exp(Ts(i)./22.9));%[km]
      catch
           Ht(i,1)=11.3; %Default [km]
      end
      
       %COMPUTE DRY MAPPING COEFFICIENTs(ah,bh,ch)        
       ah(i,1)=0.001185.*(1+(6.071e-5.*(P(i)-1000))-(1.471e-4.*es(i))+(3.072e-3.*(Ts(i)-20))+(1.965e-2.*(beta(i,1)+6.5))-(5.645e-3.*(Ht(i,1)-11.231)));
       bh(i,1)=0.001144.*(1+(1.164e-5.*(P(i)-1000))-(2.795e-4.*es(i))+(3.109e-3.*(Ts(i)-20))+(3.038e-2.*(beta(i,1)+6.5))-(1.217e-3.*(Ht(i,1)-11.231)));
       ch(i,1)=-0.0090;
                   
      %NOW THE MAPPING FUNCTION(MF)
      MFh(j,i)=1./(sin(satEL(j,i))+(ah(i)./(sin(satEL(j,i))+(bh(i)./(sin(satEL(j,i))+ch(i))))));%Hydrostatic Mapping Function
      
   end   %for j=1:nrow

end  %for i=1:ncol
%%%%%%%%%%%%%%%%%%%%%%%%END OF Davis_MF.m %%%%%%%%%%%%%%%%%%%%%%% 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%L.SUBROUTINE TO COMPUTE Saastamoinen MAPPING FUNCTION
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [MFh,MFw]=Saastamoinen_MF(satELEV)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%           "Saastamoinen_MF" is a subroutine that Computes the wet and dry*
%             Tropospheric mapping function Using Saastamoinen Model       *

%USAGE:                                                                    *                                                           *
%      [MFh,MFw]=Saastamoinen_MF(satELEV)                                  *
%INPUT:                                                                    *
%      satELEV : Satellite Elevation in [degrees]                          *

%OUTPUT:                                                                   *
%1.     MFh ==> Hydrostatic/dry mapping function                           *
%2.     MFw ==> Non Hydrostatic/wet mapping function                       *
%==========================================================================
%REFERENCE:                                                                *
%1.--------"GPS Theory and application",edited by B.Parkinson,J.Spilker,   *
%                                                         P.Enge, AIAA,1996*
%2.--------GPS Theory and Practice,Hofmann-Wellenhof, Lichtenegger, and... *
%                            Collins,Fifth, revised edition  pages 106-115.* 
%3.--------"Global Positioning System, Mishra & Enge", pg 172              *
%4.-------- Modeling of Tropospheric Delays Using ANFIS,Wayan Suparta;...  *
%                                                Kemal Maulana Alhasa(2016)*
%5.--------GPS Theory,Algorithms and Applications 2ed,Guochang Xu,June 2007*
%                                                    pg 66 eqns(5.95-5.108)*
%3.--------Seeber, G. (2003). Satellite Geodesy: Foundations, Methods &    *
%                             Applications. 2nd completely revised  Ed
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************
%***************************************************************************
%1.*****INITIALIZING OUTPUT VARIABLEs
satEL=zeros(size(satELEV,1),size(satELEV,2));%Assign zeros of nxm to SHD
[MFh,MFw,ZA]=deal(satEL);%copy the contents of satEl to all the requested outputs

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
       
       ZA(j,i) =(pi/2)-satEL(j,i);%Zenith Angle/Distance in [radian] 
             
       %NOW THE MAPPING FUNCTIONs(MF)
       try          
         MFh(j,i)=1./cos(ZA(j,i));%Hydrostatic Mapping Function
         MFw(j,i)=1./cos(ZA(j,i));%Wet Mapping Function         
      catch 
           MFh(j,i)=1./sin(satEL(j,i));%Hydrostatic Mapping Function
           MFw(j,i)=1./sin(satEL(j,i));%Wet Mapping Function
      end %\\try  

   end   %for j=1:nrow

end  %for i=1:ncol
%%%%%%%%%%%%%%%%%%%%%%%%END OF Saastamoinen_MF.m %%%%%%%%%%%%%%%%%%%%%%% 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%M.SUBROUTINE TO COMPUTE Goad & Goodman MAPPING FUNCTION
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [MFh,MFw]=Goad_Goodman_MF(T,P,es,satELEV)
%DESCRIPTION:                                                              *
%            "Goad_Goodman_MF" is a subroutine that Computes the wet and   *
%            dry Tropospheric mapping function Using  Goad & Goodman (1974)*
%            (Modified Hopfield) Model.                                    *
%SYNTAX:                                                                   *
%       [MFh,MFw]=GG_MF(Dry_ref,Wet_ref,T,P,es,satELEV)                    *
%INPUT:                                                                    *
%            P = atmospheric pressure [mbar/hPa]                           *
%            T = atmospheric temperature in degrees Kelvin [k]             *
%            es = partial pressure of water vapor [mbar/hPa]               *
%      satELEV = Satellite Elevation Angle in [radian]                     *
%      satELEV is in the format (nxm)                                      *
%Where:                                                                    *
%      n = Number of rows representing each satellite                      *
%      m = Number of columns representing each station/receiver            *

%OUTPUT:                                                                   *
%1.     MFh = Hydrostatic/dry mapping function                             *
%2.     MFw = Non Hydrostatic/wet mapping function                         *
%==========================================================================
%REFERENCE:                                                                *
%1.        Goad,C.C.& Goodman, L. (1974) A Modified Tropospheric Refraction* 
%          CorrectionModel. Paper presented at the American Geophysical    *
%          Union Annual Fall Meeting, San rancisco, December 12-17, as     *
%          presented in: "Global Positioning System, Mishra & Enge", pg 172 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
% =========================================================================
%  Written by: OSAH SAMUEL, Msc Geomatic Engineering ,2016                 *    
%       Email: osahsamuel@yahoo.ca                                         *
%         Tel:+233 (0)246137410/+233 (0)509438484                          * 
%==========================================================================
%***************************************************************************
%***************************************************************************
satEL=satELEV;%Assignment
%(1)*****INITIALIZING OUTPUT VARIABLEs
MFh=zeros(size(satEL,1),size(satEL,2));%Assign zeros of nxm to MFh
MFw=deal(MFh);%copy the contents of MFh to MFw

%(2)**COMPUTE DRY & WET TROPO REFRACTIVITY(N) @ MSL

%***USING ESSEN & FROOME REFRACTIVITY MODEL
k1=77.624; %Refractivity constant(K/mbar or K/hpa)
k2 = -12.92; %Refractivity constant(K/mbar or K/hpa)
k3=371900;%Refractivity constant(K^2/mbar or K^2/hpa)
N_dry0 = k1 .* P ./ T ; %Dry refractivity along signal path       
N_wet0 = k2.*(es./ T)+k3.*( es ./ T.^2);%Wet refractivity along signal path

%(3)******COMPUTE HEIGHT(TOP)OF TROPOSPHERE(DRY & WET COMPONENT)
try
   %Using Saastamoinen range Correction model,compute....
   %DRY "TOP" OF TROPOSPHERE
   hDry = (5*(0.002277).*P)./(N_dry0.*1.0e-6);%[meters]

   %*****WET "TOP" OF TROPOSPHERE
   hWet = (5*(0.002277)./(N_wet0.*1.0e-6)).*((1255./T)+0.05).*es;%[meters]

catch 
     %Using Hopfield range Correction model,compute....
     %DRY "TOP" OF TROPOSPHERE
     hDry = 40136 + 148.72 .* (T - 273.16);% [meters]

     %*****WET "TOP" OF TROPOSPHERE
     hWet = 11000;% [meters]
    
end %\\try
      
%(4)**************COMPUTE MAPPING FUNCTIONs(MFh & MFw)

%4.1.FIND NUMBER OF ROWS & COLUMNs IN satELEV
[nrow,ncol]=size(satEL);

%4.2.LOOP OVER USER POSITIONs & SATELLITE ELEVATIONs
for i=1:ncol %Loop over the Number of Receiver positions
    
   for j=1:nrow %Loop over the Number of Satellite Elevations
       
       %4.2.1.CONVERT ELEVATION ANGLEs TO RADIAN  
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
                     
      %COMPUTE TROPOSHERIC POLYNOMIALS(Alpha)/mapping functions
      MFh(j,i)=AlphaValues(satEL(j,i),hDry(i));%Hydrostatic Mapping Function
      MFw(j,i)=AlphaValues(satEL(j,i),hWet(i));%Wet Mapping Function
                  
   end  %\\for j=1:nrow

end  %\\for i=1:ncol
%***********************END OF GG_MF.m ***********************
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%M.1.******SUBROUTINE TO COMPUTE TROPOSHERIC POLYNOMIALS
function dr=AlphaValues(satEL,ht)
% *************************************************************************
% *************************************************************************
%DESCRIPTION:
%             "AlphaValues" Computes Tropospheric Polynomial(Alpha)Values 
%              for both Dry & Wet tropospheric components based on
%              Goad,C.C. & Goodman Tropospheric Model(Modified Hopfield).

%USAGE:dr=AlphaValues(satEL,ht)
%INPUT:
%      satEL-> Satellite Elevation Angle in [radian]
%          ht->Height(TOP)of Troposhere either dry / wet component in meters
%OUTPUT:
%       dr--->Sum of Polynomial(Alpha) values

Re=6378137;%Semi-Major Axis of Earth Ellipsoid (m)
r=sqrt((Re+ht)^2-(Re*cos(satEL))^2)-Re*sin(satEL);%Range to top of dry/wet component
a=-sin(satEL)/ht; b=-cos(satEL)^2/(2*ht*Re);%Auxiliary Constants

try
    alpha=[1,4*a,6*a^2+4*b,4*a*(a^2+3*b),a^4+12*a^2*b+6*b^2,4*a*b*(a^2+3*b),...
      b^2*(6*a^2+4*b),4*a*b^3,b^4];%Aplha values/Polynomials
    dr=sum(alpha.*(r.^(1:9))./(1:9)) ;
catch
      %Compute Aplha values/Polynomials
      rn = zeros(9,1);
      for i = 1:9
         rn(i) = r^(i); 
      end 
      alpha=[1,4*a/2,(6*a^2+4*b)/3,(4*a*(a^2+3*b))/4,(a^4+12*a^2*b+6*b^2)/5,...
           (4*a*b*(a^2+3*b))/6,(b^2*(6*a^2+4*b))/7,(4*a*b^3)/8,(b^4)/9];
       
       dr = sum(alpha*rn);  
end% try
%%%%%%%%%%%%%%%%%%%%%%%%END OF AlphaValues.m %%%%%%%%%%%%%%%%%%%%%%% 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%  --------------------------------------------------------------------
%F.SUBROUTINE TO COMPUTE ZENITH TROPOSPHERIC DELAYS FROM VMF GRID FILES
%  --------------------------------------------------------------------
function [ah,aw,VMF_model] = SearchReadVMFgrid(Time,lat,lon,h,VMFgrids)
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              * 
%            This subroutine reads and extract Vienna Mapping Function(VMF)* 
%            coefficients(ah,aw) from both gridded VMF1& VMF3 files ...    *
%            as available from: 
%            http://vmf.geo.tuwien.ac.at/trop_products/GRID/ for specific  * 
%            sites near the  Earth surface.The VMF1 grid file is based on  *                                                
%            a 2.0° × 2.5° external grid file while the VMF3 grid file     *
%            is based on either a 1° × 1° OR 5° × 5° grid.Both grids(VMF1 &*
%            VMF3) are given for four daily epochs (0h, 6h, 12h, 18h UT).  *
%            The VMF grid file consist two MF coefficients (ah,aw)as well  * 
%            as the Hydrostatic and Wet portions of ZPD(zh,zw).            *                                         *      
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

if any([all([~isfield(VMFgrids,'H00'),~isfield(VMFgrids,'H06'),~isfield(VMFgrids,'H12'),~isfield(VMFgrids,'H18')]),...
        isempty(VMFgrids)])%Determine whether array is empt
    
    
   %Call the "SearchVMFgrids.m" fxn
   [VMF_grid_found,VMFgrids] = SearchVMFgrids() ;%#ok<*ASGLU>
     
end

if any([any([isfield(VMFgrids,'H00'),isfield(VMFgrids,'H06'),isfield(VMFgrids,'H12'),isfield(VMFgrids,'H18')]),...
        ~isempty(VMFgrids)])%Determine whether array is empt
    
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
   
else
    H0files = []; H6files = []; H12files = []; H18files = [];
    
end
        
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
VMFgrid_type = getappdata(0,'VMFgrid_type_mf');%VMF GRID TYPE(VMF1 or VMF3)
VMFgrid_res  = getappdata(0,'VMFgrid_res_mf');%VMF GRID RESOLUTION especially for VMF3(1°x1° OR 5°x5°)

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


%CHECK IF FILES ARE NOT EMPTY([])
%NOTE: 
%     IF IMPORTED VMF GRID FILES DOES NOT MATCH WITH PROCESS DATA TIME
%     THE SYSTEM WILL RETURN EMPTY([]) MATRICES FOR OUTPUT
if all([isempty(H0filesN),isempty(H6filesN),isempty(H12filesN),isempty(H18filesN)])
    
   ah = [];
   aw = [];
   VMF_model = [];
   
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
                  
        %LOOP OVER VMFG FILES
        for k = 1:size(VMFGfileP,1)
               
            VMFgridP = VMFGfileP{k,1};%CHANGE FROM CELL ARRAY TO CHARACTER
              
            %EXTRACT ZENITH DELAYS & CREATE AN ARRAY OF II(POSITION EPOCHs) ROWS,k(VMF EPOCHs) COLUMNs    
            %*****CHECK VMF FILE TYPE(VMF1/VMF3)getTROPdelays
            if any([strcmpi(VMFgrid_type,'VMF1'),strfind(VMFgrid_type,'VMF1'),any(VMFgrid_res ==[ 2,2.5])]) 
                   
               VMF_model = 'VMF1'; %flag to indicate VMF file type
                  
               [ah(i,k),aw(i,k)] = readVMF1grid(date(i,1:3),lat(i),lon(i),h(i),VMFgridP,oro);
                                                                  
            elseif any([strcmpi(VMFgrid_type,'VMF3'),strfind(VMFgrid_type,'VMF3'),any(VMFgrid_res ==[ 1,5])]) 
                   
                   VMF_model = 'VMF3';%flag to indicate VMF file type
                      
                   [ah(i,k),aw(i,k)] = readVMF3grid(date(i,1:3),lat(i),lon(i),h(i),VMFgridP,VMFgrid_res,oro);
                          
            end  %//if any([strcmpi(VMFgrid_type,'VMF1'),strfind(VMFgrid_type,'VMF1'),any(VMFgrid_res ==[ 2,2.5])]) 
                        

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
   ah_mean  = (round(mean(ah,2),5))';
   aw_mean  = (round(mean(aw,2),5))';
    
   %ASSIGNMENT
   ah = ah_mean;
   aw = aw_mean;
        
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
                  
                  [AH(j,k),AW(j,k)] = readVMF1grid(date(j,1:3),lat(i),lon(i),h(i),VMFgridP,oro);
                                                                  
               elseif any([strcmpi(VMFgrid_type,'VMF3'),strfind(VMFgrid_type,'VMF3'),any(VMFgrid_res ==[ 1,5])])
                   
                      VMF_model = 'VMF3';%flag to indicate VMF file type
                      
                      [AH(j,k),AW(j,k)] = readVMF3grid(date(j,1:3),lat(i),lon(i),h(i),VMFgridP,VMFgrid_res,oro); 
                            
               end %//if any([strcmpi(VMFgrid_type,'VMF1'),strfind(VMFgrid_type,'VMF1'),any(VMFgrid_res ==[ 2,2.5])]) 
                   
              %**************CREATE NUMBER OF COLUMNs IN OUTPUT FILEs
              %CREATE A 3D CELL ARRAY FOR THE OUTPUT
              %Create a cell array of No.of observational period/days/date/
              %time(j) as rows,No.of grid epoch[00,06,12,18 UT] files as...
              %columns(k) and the No. of User stations(i) as number of arrays.
              %i.e.Create a j x k x i cell array matrices.
              %------------------------------------------------------------
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
   
    %VMF COEFFICIENTS
    ah_epoch = cell2mat(Ah(:,:,q));
    aw_epoch = cell2mat(Aw(:,:,q));
    
    ZPD_epochs = [ah_epoch aw_epoch];
    
    %*****FIND MEAN OF ALL AVAILABLE EPOCHS[00 06 12 18] FOR DAILY ZPDs
    %I.E. FIND MEAN OF EACH ROW IN ZTD, ZHD & ZWD (MEAN OF VMF EPOCHs ...
    %                                     TO OBTAIN TOTAL DELAY FOR A DAY)
    %NOTE:
    %     THE NUMBER ROWS IN EACH MATRIX REPRESENTS THE NUMBERS OF
    %     DAYS(TIME/DATE) & NUMBER OF COLUMNS ALSO REPRESENTS THE NO. OF
    %     VMF EPOCHS[i.e. 00 06 12 18]
    %----------------------------------------------------------------------
    %ASSIGNMENT
    
    %***************VMF COEFFICIENTS
    %FIND MEAN & REFORMAT COLUMNS INTO DAYS(TIME/DATE)
    ah_mean  = (mean(cell2mat(Ah(:,:,q)),2))';
    aw_mean  = (mean(cell2mat(Aw(:,:,q)),2))'; 
    
    %ASSIGNMENT
    ah(q,:) = ah_mean;
    aw(q,:) = aw_mean;
       
end %//for q = 1 : nrow_Rpos  

ZPD_mfc = [ah aw]; 

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
                  
                  [AH(j,k),AW(j,k)] = readVMF1grid(date(j,1:3),lat(i),lon(i),h(i),VMFgridP,oro);
                                                                  
               elseif any([strcmpi(VMFgrid_type,'VMF3'),strfind(VMFgrid_type,'VMF3'),any(VMFgrid_res ==[ 1,5])])
                   
                      VMF_model = 'VMF3';%flag to indicate VMF file type
                      
                     [AH(j,k),AW(j,k)] = readVMF3grid(date(j,1:3),lat(i),lon(i),h(i),VMFgridP,VMFgrid_res,oro); 
                            
               end %//if any([strcmpi(VMFgrid_type,'VMF1'),strfind(VMFgrid_type,'VMF1'),any(VMFgrid_res ==[ 2,2.5])]) 
                   
              %**************CREATE NUMBER OF COLUMNs IN OUTPUT FILEs
              %CREATE A 3D CELL ARRAY FOR THE OUTPUT
              %Create a cell array of No.of observational period/days/date/
              %time(j) as rows,No.of grid epoch[00,06,12,18 UT] files as...
              %columns(k) and the No. of User stations(i) as number of arrays.
              %i.e.Create a j x k x i cell array matrices.
              %------------------------------------------------------------
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
   
    
    %VMF COEFFICIENTS
    ah_epoch = cell2mat(Ah(:,:,q));
    aw_epoch = cell2mat(Aw(:,:,q));
    
    ZPD_epochs = [ah_epoch aw_epoch];
    
    %*****FIND MEAN OF ALL AVAILABLE EPOCHS[00 06 12 18] FOR DAILY ZPDs
    %I.E. FIND MEAN OF EACH ROW IN ZTD, ZHD & ZWD (MEAN OF VMF EPOCHs ...
    %                                     TO OBTAIN TOTAL DELAY FOR A DAY)
    %NOTE:
    %     THE NUMBER ROWS IN EACH MATRIX REPRESENTS THE NUMBERS OF
    %     DAYS(TIME/DATE) & NUMBER OF COLUMNS ALSO REPRESENTS THE NO. OF
    %     VMF EPOCHS[i.e. 00 06 12 18]
    %----------------------------------------------------------------------
    %ASSIGNMENT
    
    %***************VMF COEFFICIENTS
    %FIND MEAN & REFORMAT COLUMNS INTO DAYS(TIME/DATE)
    ah_mean  = (mean(cell2mat(Ah(:,:,q)),2))';
    aw_mean  = (mean(cell2mat(Aw(:,:,q)),2))'; 
    
    %ASSIGNMENT
    ah(q,:) = ah_mean;
    aw(q,:) = aw_mean;
       
end %//for q = 1 : nrow_Rpos  

ZPD_mfc = [ah aw]; 

end %//if isequal(nrow_time,nrow_pos)

%%%%%%%%%%%%%%%%%%%%%%%%END OF ZPD EXTRACTION %%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%********SUB-ROUTINE TO READ  AND EXTRACT ZENITH TROPOSPHERIC DELAYS

%(F.1)SUBROUTINE TO READ VMF1 GRID FILE & EXTRACT COEFFICIENTs(ZHD,ZWD,ZTD,ah,aw)
function [ah,aw] = readVMF1grid(UTCtime,Lat,Lon,hgt,VMF1_grid,...
                                                            orography_file)                                                     
%**************************************************************************
%DESCRIPTION:
%***********This subroutine determines the Hydrostatic and Wet Mapping ... * 
%           Function(MF)Coefficients ah and aw from the gridded VMF1 files,* 
%           as available from:                                             *
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
%      [ah,aw] = Readvmf_grid(UTCtime,lat,lon,hell,VMF1_grid...
%                                                          ,orography_file)
%INPUTs:
%1.     UTCtime:.........UTC time in [Year,Month,Day,Hour,Minute,Seconds]  *
%2.     Lat: ............station ellipsoidal latitude in [degrees]         *
%3.     Lon: ............station ellipsoidal longitude in [degrees]        *
%4.     hgt: ............station ellipsoidal height in [meters]            *
%5.     VMF1grid_file:...VMF grid file eg:'VMFG_20180101.H00'              *
%6.     orography_file:...ellipsoidal orography. eg:'orography_ell'        *

%OUTPUTs:                                                                  *
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
    
   hcoe = VMF1_data_int_h1(1,8);
   wcoe = VMF1_data_int_h1(1,9);

else
    %BILINEAR INTERPOLATION (INTERPRETED AS TWO 1D LINEAR INTERPOLATIONS 
    %FOR LAT AND LON, BUT PROGRAMMED WITHOUT SUBFUNCTIONS) 
    
    %(a)*******LINEAR INTERPOLATION FOR LONGITUDE
    if ~isequal(VMF1_data_int_h1(1,2), VMF1_data_int_h1(2,2))%if longitude must be interpolated (that is, the point does not have a longitude on the interval [0:2.5:357.5])
        
         ah_lon1 = VMF1_data_int_h1(1,8) + (VMF1_data_int_h1(2,8)-VMF1_data_int_h1(1,8))*(lonDeg-VMF1_data_int_h1(1,2))/(VMF1_data_int_h1(2,2)-VMF1_data_int_h1(1,2));
         ah_lon2 = VMF1_data_int_h1(3,8) + (VMF1_data_int_h1(4,8)-VMF1_data_int_h1(3,8))*(lonDeg-VMF1_data_int_h1(3,2))/(VMF1_data_int_h1(4,2)-VMF1_data_int_h1(3,2));
         aw_lon1 = VMF1_data_int_h1(1,9) + (VMF1_data_int_h1(2,9)-VMF1_data_int_h1(1,9))*(lonDeg-VMF1_data_int_h1(1,2))/(VMF1_data_int_h1(2,2)-VMF1_data_int_h1(1,2));
         aw_lon2 = VMF1_data_int_h1(3,9) + (VMF1_data_int_h1(4,9)-VMF1_data_int_h1(3,9))*(lonDeg-VMF1_data_int_h1(3,2))/(VMF1_data_int_h1(4,2)-VMF1_data_int_h1(3,2));

    else %if the station coincides with the longitude of the grid
        
        ah_lon1 = VMF1_data_int_h1(1,8);
        ah_lon2 = VMF1_data_int_h1(3,8);
        aw_lon1 = VMF1_data_int_h1(1,9);
        aw_lon2 = VMF1_data_int_h1(3,9);
    end
    
    %*****LINEAR INTERPOLATION FOR LATITUDE
    if ~isequal(VMF1_data_int_h1(1,1), VMF1_data_int_h1(3,1))
        
        hcoe = ah_lon1 + (ah_lon2-ah_lon1)*(latDeg-VMF1_data_int_h1(1,1))/(VMF1_data_int_h1(3,1)-VMF1_data_int_h1(1,1));
        wcoe = aw_lon1 + (aw_lon2-aw_lon1)*(latDeg-VMF1_data_int_h1(1,1))/(VMF1_data_int_h1(3,1)-VMF1_data_int_h1(1,1));
   
    else %if the station coincides with the latitude of the grid
        
        hcoe = ah_lon1;%Hydrostatic Coefficient
        wcoe = aw_lon1;%Wet Coefficient
    end
      
end %//if length(unique(index)) == 1 

%*******ASSIGNMENT
ah=hcoe;%Hydrostatic Coefficient
aw=wcoe;%Wet Coefficient


%%%%%%%%%%%%%%%%%%%%%%%%%END OF readVMF1grid.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-     
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%(F.2).SUBROUTINE TO READ VMF3 GRID FILE & EXTRACT COEFFICIENTs(zhd,zwd,ztd,ah,aw)
function [ah,aw] = readVMF3grid(UTCtime,Lat,Lon,hgt,VMF3_grid,...
                                                            gridRESOLUTION,orography_file)                                                     
%**************************************************************************
%DESCRIPTION:
%***********This subroutine determines the Hydrostatic and Wet Mapping ... * 
%           Function(MF)Coefficients ah and aw from the gridded VMF3 files, 
%           as available from:                                             *
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
%      [ah,aw] = Readvmf_grid(UTCtime,lat,lon,hell,VMF1_grid...
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
    
    hcoe = VMF3_data_int_h1(1,8);
    wcoe = VMF3_data_int_h1(1,9);

else
    %BILINEAR INTERPOLATION (INTERPRETED AS TWO 1D LINEAR INTERPOLATIONS 
    %FOR LAT AND LON, BUT PROGRAMMED WITHOUT SUBFUNCTIONS) 
    
    %(a)*******LINEAR INTERPOLATION FOR LONGITUDE
    if ~isequal(VMF3_data_int_h1(1,2), VMF3_data_int_h1(2,2))%if longitude must be interpolated (that is, the point does not have a longitude on the interval [0:2.5:357.5])
         ah_lon1 = VMF3_data_int_h1(1,8) + (VMF3_data_int_h1(2,8)-VMF3_data_int_h1(1,8))*(lonDeg-VMF3_data_int_h1(1,2))/(VMF3_data_int_h1(2,2)-VMF3_data_int_h1(1,2));
         ah_lon2 = VMF3_data_int_h1(3,8) + (VMF3_data_int_h1(4,8)-VMF3_data_int_h1(3,8))*(lonDeg-VMF3_data_int_h1(3,2))/(VMF3_data_int_h1(4,2)-VMF3_data_int_h1(3,2));
         aw_lon1 = VMF3_data_int_h1(1,9) + (VMF3_data_int_h1(2,9)-VMF3_data_int_h1(1,9))*(lonDeg-VMF3_data_int_h1(1,2))/(VMF3_data_int_h1(2,2)-VMF3_data_int_h1(1,2));
         aw_lon2 = VMF3_data_int_h1(3,9) + (VMF3_data_int_h1(4,9)-VMF3_data_int_h1(3,9))*(lonDeg-VMF3_data_int_h1(3,2))/(VMF3_data_int_h1(4,2)-VMF3_data_int_h1(3,2));

    else %if the station coincides with the longitude of the grid
        ah_lon1 = VMF3_data_int_h1(1,8);
        ah_lon2 = VMF3_data_int_h1(3,8);
        aw_lon1 = VMF3_data_int_h1(1,9);
        aw_lon2 = VMF3_data_int_h1(3,9);
        
    end
    
    %*****LINEAR INTERPOLATION FOR LATITUDE
    if ~isequal(VMF3_data_int_h1(1,1), VMF3_data_int_h1(3,1))
        hcoe = ah_lon1 + (ah_lon2-ah_lon1)*(latDeg-VMF3_data_int_h1(1,1))/(VMF3_data_int_h1(3,1)-VMF3_data_int_h1(1,1));
        wcoe = aw_lon1 + (aw_lon2-aw_lon1)*(latDeg-VMF3_data_int_h1(1,1))/(VMF3_data_int_h1(3,1)-VMF3_data_int_h1(1,1));
    else %if the station coincides with the latitude of the grid
         hcoe = ah_lon1;%Hydrostatic Coefficient
         wcoe = aw_lon1;%Wet Coefficient
    end
      
end 
%ASSIGN VALUEs FOR ALL STATIONs
ah  = hcoe;%Hydrostatic Coefficient
aw  = wcoe;%Wet Coefficient

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





