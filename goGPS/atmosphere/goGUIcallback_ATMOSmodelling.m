%**************************************************************************
%DESCRIPTION:                                                              * 
%           "goGUIcallback_ATMOSmodelling" is a Sub-routine that links the *
%           various callbacks for Atmospheric Delay Modelling GUI          *
%           components in goGPS.                                           *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================
% WRITTEN BY:                                                              +
%            OSAH SAMUEL, MSC GEOMATIC ENGINEERING (PhD STUDENT)           +
%            Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com                +
%            Phone:+233(0)246137410 / +233(0)509438484                     + 
%==========================================================================
%**************************************************************************+
%**************************************************************************+
function goGUIcallback_ATMOSmodelling(hObject,handles,gui_goGPS_OpeningFcn)


if nargin == 2
    
  gui_goGPS_OpeningFcn = ' goGPS_OpeningFcn'; 
  
end
      
%LOAD PREVIOUS goGPS ATMOSPHERIC MODELLING GUI COMPONENTS STATE

if nargin == 3
    
   if gui_goGPS_OpeningFcn
    
      %LOAD PREVIOUS GUI STATE
      loadGUIState(handles)
   end 

elseif any(nargin==[1,2])

%                 -------------------------------------
%*****************CALLING OF VARIOUS CALLBACK FUNCTIONS
%                 --------------------------------------
%1.******IONOSPHERIC CORRECTION CALLBACK
if strcmpi(get(hObject,'Tag'),'pop_lIono_correction')
   IONOcorrection_Callback(hObject,handles)
end

%2.******TROPOSPHERIC CORRECTION CALLBACK
%2.1 COMBINE MODELs(RADIO BUTTON(rb)) 
if strcmpi(get(hObject,'Tag'),'rb_Tropo_combine')
    TROPOcorrection_combine_rb_Callback(hObject,handles)
end

%2.2 SEPARATE MODELs(RADIO BUTTON)
if strcmpi(get(hObject,'Tag'),'rb_Tropo_separate')
   TROPOcorrection_separate_rb_Callback(hObject,handles)
end

%2.3 COMBINE MODELs(POPUP MENU) 
if strcmpi(get(hObject,'Tag'),'pop_lTropo_combine')
   TROPOmodels_combine_popm_Callback(hObject, handles)
end


%2.4 HYDROSTATIC MODELs(POPUP MENU) 
if strcmpi(get(hObject,'Tag'),'pop_ITropo_hydrostatic')
   TROPOmodels_hydrostatic_popm_Callback(hObject, handles)
end

%2.5 WET MODELs(POPUP MENU) 
if strcmpi(get(hObject,'Tag'),'pop_ITropo_wet')
   TROPOmodels_wet_popm_Callback(hObject, handles)
end 

%2.6 *****************MAPPING FUNCTIONS
%2.6.1 RADIO BUTTONS
%A.USE MODEL MAPPING FUNCTION
if strcmpi(get(hObject,'Tag'),'rb_model_MF')
   TROPOmodels_mappingfunction_rb_Callback(hObject,handles)
end 

%B.USE OTHER MAPPING MODELS
if strcmpi(get(hObject,'Tag'),'rb_different_MF')
   OTHERmodels_mappingfunction_rb_Callback(hObject,handles)
end 


%2.6.2 POPUP MENUs 
%A.HYDROSTATIC MAPPING FUNCTIONS
if strcmpi(get(hObject,'Tag'),'pop_IHydrostic_MF')
   HYDROSTATIC_mappingfunction_popm_Callback(hObject,handles)
end 

%B.WET MAPPING FUNCTIONS CALLBACK
if strcmpi(get(hObject,'Tag'),'pop_IWet_MF')
   WET_mappingfunction_popm_Callback(hObject,handles)
end 

%2.7 TROPOSPHERIC DELAY ESTIMATION USING PPP
if strcmpi(get(hObject,'Tag'),'flag_tropo')
   PPP_tropoESTIMATION_Callback(hObject,handles)
end 

%3.SOURCE OF METEOROLOGICAL PARAMETERS
%A.POPUP MENU CALLBACK
if strcmpi(get(hObject,'Tag'),'pop_source_metpara')
   METparameters_popm_Callback(hObject,handles)
end 


%B.**************SELECTION OF GRID RESOLUTION(RADIO BUTTON)
%1° GRID RESOLUTION CALLBACK
if strcmpi(get(hObject,'Tag'),'rb_grid_resolution_1')
 GridResolution_1_rb_Callback(hObject,handles)
end

%5° GRID RESOLUTION CALLBACK
if strcmpi(get(hObject,'Tag'),'rb_grid_resolution_5')
 GridResolution_5_rb_Callback(hObject,handles)
end

%4.*****ADD VMD GRIDDED ZENITH DELAYS TO THE OUPT FILE
if strcmpi(get(hObject,'Tag'),'cb_extract_VMF_ZTDs')
   extract_VMF_ZTDs_cb_Callback(hObject,handles)
end

%5.*****RETRIEVE PRECIPITABLE WATER VAPOUR(PWV) CHECKBOX(cb)
if strcmpi(get(hObject,'Tag'),'cb_retrievePWV')
   retrievePWV_cb_Callback(hObject,handles)
end

%5.A*****WEIGHTED MEAN TEMPERATURE(Tm) MODELS POPUP MENU
if strcmpi(get(hObject,'Tag'),'pop_Tm')
   Tm_model_pop_Callback(hObject,handles)
end

%5.B.**************SELECTION OF GRID RESOLUTION(RADIO BUTTON)
%1° GRID RESOLUTION CALLBACK
if strcmpi(get(hObject,'Tag'),'rb_grid_resolution_1_pwv')
   GridResolution_1_rb_pwv_Callback(hObject,handles)
end

%5° GRID RESOLUTION CALLBACK
if strcmpi(get(hObject,'Tag'),'rb_grid_resolution_5_pwv')
   GridResolution_5_rb_pwv_Callback(hObject,handles)
end

%6.*************THE go ! PUSH BUTTON(pb) CALLBACK
if strcmpi(get(hObject,'Tag'),'go_button')
   go_pb_Callback(hObject,handles)
end

end
%--------------------------------------------------------------------------+
%******************************SUB-ROUTINES
%--------------------------------------------------------------------------+

%1.SUB-ROUTINES THAT LINKS THE VARIOUS CALLBACKS FOR ATMOSPHERIC DELAY  
%  MODELLING GUI COMPONETS IN goGPS========================================            

%                        -----------------------
%(1)======================IONOSPHERIC CORRECTION===========================
%                        ----------------------
% --- EXECUTES ON SELECTION CHANGE in pop_lIono_correction.
function IONOcorrection_Callback(hObject,handles)
%--------------------------------------------------------------------------
getconts = cellstr(get(hObject,'String'));%get IIono popup menu contents as cell array
IonoModel=getconts{get(hObject,'Value')};%get selected item from pop_lIono_correction
set(hObject,'TooltipString',strjoin({'Solution using',IonoModel,'model'}))%change Tooltip any time user select Model type
if get(hObject,'Value')==1
   set(hObject,'TooltipString',strjoin({'Solution without','Ionospheric correction'}))%change Tooltip any time user select Model type    
end

%SAVE SELECTED IONO MODEL
setappdata(0,'Iono_model',IonoModel)

%==================================END OF IONOcorrection_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%                          -----------------------
%(2)=======================TROPOSPHERIC DELAY CORRECTION===================
%                          ------------------------

%********PPP TROPOSPHERIC DELAY ESTIMATION
% --- EXECUTES ON BUTTON PRESS IN flag_tropo.
function PPP_tropoESTIMATION_Callback(hObject,handles)

if get(hObject,'value') == 1    
    
   %SAVE goGPS GUI FIGURE HANDLES
   setappdata(0,'handle',handles);

   %OPEN PPP  MAPPING FUNCTION SELECTION DIALOGUE
   %Call the 'MF_ppp.m' fxn
   MF_ppp()
   
end

%           --------------------------------
%2.1********COMBINE MODELs(RADIO BUTTON(rb)) 
%           --------------------------------
% --- EXECUTES ON BUTTON PRESS IN rb_Tropo_combine.
function TROPOcorrection_combine_rb_Callback(hObject,handles)
%--------------------------------------------------------------------------
global tropoModel 

%********ENABLE SOME TOOLS
%1.COMBINE TROPO MODELs
set(hObject,'Value',1,'TooltipString','Solution with Combined Tropospheric Correction [SHD + SWD]')%Tropo_combine Radio Button
set(handles.pop_lTropo_combine,'enable','on')%Tropo combine popup menu
set(handles.text_Tropo_combine,'enable','on')%Tropo combine static text
%SET TOOLTIP STRING
getconts = cellstr(get(handles.pop_lTropo_combine,'String'));%get Tropo_combine popup menu contents as cell array
tropoModel=getconts{get(handles.pop_lTropo_combine,'Value')};%get selected item from pop_lTropo_combine
tropoModelVAL=get(handles.pop_lTropo_combine,'Value');%GET SELECTED VALUE
set(handles.pop_lTropo_combine,'TooltipString',strjoin({'Solution using',tropoModel,'model'}))%change Tooltip any time user select Model type

if isequal(get(handles.pop_lTropo_combine,'value'),1)
   set(handles.pop_lTropo_combine,'TooltipString',strjoin({'Solution without','Tropospheric correction'}))%change Tooltip any time user select Model type    
end

%*************DISABLE SOME TOOLs
%1.HYDROSTATIC & WET TROPO MODELs
set(handles.rb_Tropo_separate,'Value',0,'TooltipString','Select for Separate Tropospheric Correction');%Tropo_separate Radio Button
set(handles.pop_ITropo_hydrostatic,'enable','off','TooltipString','')%Tropo_hydrostatic popup menu
set(handles.pop_ITropo_wet,'enable','off','TooltipString','')%Tropo_wet popup menu
set(handles.text_Tropo_hydrostatic,'enable','off')%Tropo_hydrostatic static text
set(handles.text_Tropo_wet,'enable','off')%Tropo_wet static text
%2.MAPPING FUNCTION
set(handles.text_mapping_function,'enable','off');%Mapping function static text
set(handles.rb_model_MF,'enable','off','Value',0,'TooltipString','');%Model mapping function(MF) Radio Button
set(handles.rb_different_MF,'enable','off','Value',0,'TooltipString','');%Different MF Model Radio Button
set(handles.pop_IHydrostic_MF,'enable','off','TooltipString','')%Hydrostatic MF popup menu
set(handles.pop_IWet_MF,'enable','off','TooltipString','')%Wet MF popup menu
set(handles.text__Hydrostatic_MF,'enable','off');%Hydrostatic MF static text
set(handles.text__Wet_MF,'enable','off');%Hydrostatic MF static text

%SAVE SELECTED COMBINE TROPO MODEL
setappdata(0,'Tropo_Cmodel',tropoModel)

%SAVE OPTIONs FOR TROPO MODELLING
setappdata(0,'option_Cmodels',get(hObject,'Value'))%save option button for combine tropo models
setappdata(0,'option_Smodels',get(handles.rb_Tropo_separate,'Value'))%save option button for separate tropo models

%******SET MET PARAMETER SOURCE UICONTROLS VISIBLE/ENABLE ON

%GET SELECTED METEOROLOGICAL PARAMETER SOURCE
getconts_MET = cellstr(get(handles.pop_source_metpara,'String'));%get Tropo_combine popup menu contents as cell array
METpara = getconts_MET{get(handles.pop_source_metpara,'Value')};%get selected item from pop_lTropo_combine
metVAL=get(handles.pop_source_metpara,'Value');%GET SELECTED VALUE

%IF ANY OF THESE MODEL IS SELECTED,SET SOME UICONTROLS VISIBLE/ENABLE ON
if any([any([strncmpi(tropoModel,'no model',8),strncmpi(tropoModel,'GPT2 (5° x 5°)',14),strncmpi(tropoModel,'GPT2w (1° x 1°)',14),strncmpi(tropoModel,'GPT2w (5° x 5°)',14),...
             strncmpi(tropoModel,'GPT3 (1° x 1°)',14),strncmpi(tropoModel,'GPT3 (5° x 5°)',14),strncmpi(tropoModel,'UNB3m',5),...
             strncmpi(tropoModel,'EGNOS',5),strncmpi(tropoModel,'MOPS',4),strncmpi(tropoModel,'VMF gridded ZTD',15),strncmpi(tropoModel,'GTrop [Sun et al 2019]',22)]),...
             any(tropoModelVAL==[1,9,10,11,12,13,14,15,16,17,18])])
                 
   if strcmpi(get(handles.pop_source_metpara,'enable'),'On')
      set(handles.pop_source_metpara,'enable','off')
      set(handles.text_source_metpara,'enable','off')
 
   end
   
   if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
           any(metVAL==[7,8])])
       
      if strcmpi(get(handles.pop_met_manual,'visible'),'On')
          set(handles.pop_met_manual,'visible','off')
          set(handles.text_metmanual_source,'visible','Off')
      end
      
   elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
           any(metVAL==[4,5])])
          
          if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'on'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'on')])
             set(handles.rb_grid_resolution_1,'Visible','off')
             set(handles.rb_grid_resolution_5,'Visible','off')
             set(handles.text_grid_resolution,'Visible','off')  
          end 
   end
   
else
    
    %SET UICONTROLS VISIBLE/ENABLE ON
    set(handles.pop_source_metpara,'enable','on')
    set(handles.text_source_metpara,'enable','on')
    
    if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
            any(metVAL==[7,8])])
        
       if strcmpi(get(handles.pop_met_manual,'visible'),'Off') 
          set(handles.pop_met_manual,'visible','On')
          set(handles.text_metmanual_source,'visible','On')
       end
     
    elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
           any(metVAL==[4,5])])
          
          if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'off'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'off')])
             set(handles.rb_grid_resolution_1,'Visible','on')
             set(handles.rb_grid_resolution_5,'Visible','on')
             set(handles.text_grid_resolution,'Visible','on')  
          end
          
    end
       
end

%IF 'VMF gridded ZTD' EITHER SELECTED OR NOT, DO THE FF:
if ~any([strncmpi(tropoModel,'VMF gridded ZTD',15),tropoModelVAL==17])
   set(handles.cb_extract_VMF_ZTDs,'enable','on')
else
   set(handles.cb_extract_VMF_ZTDs,'value',0,'enable','off')
end

%================================END OF TROPOcorrection_combine_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%                       ----------------------------- 
%(2.2)******************SEPARATE MODELs(RADIO BUTTON)
%                       -----------------------------
% --- EXECUTES ON BUTTON PRESS IN rb_Tropo_separate.
function TROPOcorrection_separate_rb_Callback(hObject,handles)
%--------------------------------------------------------------------------
global dryModel wetModel MFh_model MFw_model 
%**************ENABLE SOME TOOLS
%1.HYDROSTATIC & WET TROPO MODELs
set(hObject,'Value',1,'TooltipString','Solution with Separate Tropospheric Corrections [ZHD,ZWD]');%Tropo_separate Radio Button
set(handles.pop_ITropo_hydrostatic,'enable','on')%Tropo_hydrostatic popup menu
set(handles.pop_ITropo_wet,'enable','on')%Tropo_wet popup menu
set(handles.text_Tropo_hydrostatic,'enable','on')%Tropo_hydrostatic static text
set(handles.text_Tropo_wet,'enable','on')%Tropo_wet static text

%SET TOOLTIP STRING
%1.HYDROSTATIC
getcontsH = cellstr(get(handles.pop_ITropo_hydrostatic,'String'));%get Tropo_hydrostatic popup menu contents as cell array
dryModel=getcontsH{get(handles.pop_ITropo_hydrostatic,'Value')};%get selected item from pop_Tropo_hydrostatic
set(handles.pop_ITropo_hydrostatic,'TooltipString',strjoin({'Solution using',dryModel,'model'}))%change Tooltip any time user select Model type

%2.WET
getcontsW = cellstr(get(handles.pop_ITropo_wet,'String'));%get Tropo_wet popup menu contents as cell array
wetModel=getcontsW{get(handles.pop_ITropo_wet,'Value')};%get selected item from pop_Tropo_wet
set(handles.pop_ITropo_wet,'TooltipString',strjoin({'Solution using',wetModel,'model'}))%change Tooltip any time user select Model type

%2.MAPPING FUNCTION
set(handles.text_mapping_function,'enable','on');%Mapping function static text
set(handles.rb_model_MF,'enable','on','Value',0,'TooltipString','Select to apply Model Mapping Function');%Model mapping function(MF) Radio Button
set(handles.rb_different_MF,'enable','on','Value',1,'TooltipString','Solution with different Mapping Function Models');%Different MF Model Radio Button
set(handles.pop_IHydrostic_MF,'enable','on')%Hydrostatic MF popup menu
set(handles.pop_IWet_MF,'enable','on')%Wet MF popup menu
set(handles.text__Hydrostatic_MF,'enable','on');%Hydrostatic MF static text
set(handles.text__Wet_MF,'enable','on');%Hydrostatic MF static text

%1.HYDROSTATIC MF
getcontsHMF = cellstr(get(handles.pop_IHydrostic_MF,'String'));%get Hydrostatic MF popup menu contents as cell array
MFh_model=getcontsHMF{get(handles.pop_IHydrostic_MF,'Value')};%get selected item from pop_IHydrostic_MF
set(handles.pop_IHydrostic_MF,'TooltipString',strjoin({'Solution using',MFh_model,'model'}))%change Tooltip any time user select Model type
%2.WET MF
getcontsWMF = cellstr(get(handles.pop_IWet_MF,'String'));%get Wet MF popup menu contents as cell array
MFw_model=getcontsWMF{get(handles.pop_IWet_MF,'Value')};%get selected item from pop_IWet_MF
set(handles.pop_IWet_MF,'TooltipString',strjoin({'Solution using',MFw_model,'model'}))%change Tooltip any time user select Model type

%*************DISABLE SOME TOOLs
%1.COMBINE TROPO MODELs
set(handles.rb_Tropo_combine,'Value',0,'TooltipString','Select for Combined Tropospheric Correction')%Tropo_combine Radio Button
set(handles.pop_lTropo_combine,'enable','off','TooltipString','')%Tropo combine popup menu
set(handles.text_Tropo_combine,'enable','off')%Tropo combine static text

% %SAVE SELECTED SEPARATE TROPO MODELs
 setappdata(0,'dryModel',dryModel)
 setappdata(0,'wetModel',wetModel)
 
%SAVE SELECTED MAPPING FUNCTION MODELs
setappdata(0,'MFh_model',MFh_model)
setappdata(0,'MFw_model',MFw_model)

%SAVE OPTIONs FOR TROPO MODELLING
setappdata(0,'option_Smodels',get(hObject,'Value'))%save option button for separate tropo models
setappdata(0,'option_Cmodels',get(handles.rb_Tropo_combine,'Value'))%save option button for combine tropo models

%SAVE OPTIONs FOR USE OF MODELs / OTHER MAPPING FUNCTION
setappdata(0,'option_model_MF',get(handles.rb_model_MF,'Value'))%save option button for models MF
setappdata(0,'option_different_MF',get(handles.rb_different_MF,'Value'))%save option button for different MF models

%******SET MET PARAMETER SOURCE UICONTROLS VISIBLE/ENABLE ON

%GET SELECTED METEOROLOGICAL PARAMETER SOURCE
getconts_MET = cellstr(get(handles.pop_source_metpara,'String'));%get MET popup menu contents as cell array
METpara = getconts_MET{get(handles.pop_source_metpara,'Value')};%get selected item from MET POPUP MENU
metVAL=get(handles.pop_source_metpara,'Value');%GET SELECTED VALUE

%GET SELECTED WET TROPO MODEL
getconts_wet = cellstr(get(handles.pop_ITropo_wet,'String'));%get Tropo_combine popup menu contents as cell array
wetModel = getconts_wet{get(handles.pop_ITropo_wet,'Value')};%get selected item from pop_lTropo_combine

%IF ANY OF THESE MODEL IS SELECTED,SET SOME UICONTROLS VISIBLE/ENABLE ON
if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
        strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),...
        strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),strncmpi(dryModel,'VMF gridded ZHD',15),strncmpi(dryModel,'GTrop [Sun et al 2019]',22)]) & ...   
    any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
        strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),...
        strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4),strncmpi(wetModel,'VMF gridded ZWD',15),strncmpi(wetModel,'GTrop [Sun et al 2019]',22)])
    
   if strcmpi(get(handles.pop_source_metpara,'enable'),'On')
      set(handles.pop_source_metpara,'enable','off')
      set(handles.text_source_metpara,'enable','off')
 
   end
   
   if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
           any(metVAL==[7,8])])
       
      if strcmpi(get(handles.pop_met_manual,'visible'),'On')
          set(handles.pop_met_manual,'visible','off')
          set(handles.text_metmanual_source,'visible','Off')
      end
      
   elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
           any(metVAL==[4,5])])
          
          if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'on'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'on')])
             set(handles.rb_grid_resolution_1,'Visible','off')
             set(handles.rb_grid_resolution_5,'Visible','off')
             set(handles.text_grid_resolution,'Visible','off')  
          end 
   end
   
else
    
    %SET UICONTROLS VISIBLE/ENABLE ON
    set(handles.pop_source_metpara,'enable','on')
    set(handles.text_source_metpara,'enable','on')
    
    if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
            any(metVAL==[7,8])])
        
       if strcmpi(get(handles.pop_met_manual,'visible'),'Off') 
          set(handles.pop_met_manual,'visible','On')
          set(handles.text_metmanual_source,'visible','On')
       end
     
    elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
           any(metVAL==[4,5])])
          
          if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'off'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'off')])
             set(handles.rb_grid_resolution_1,'Visible','on')
             set(handles.rb_grid_resolution_5,'Visible','on')
             set(handles.text_grid_resolution,'Visible','on')  
          end
          
    end
       
end

%IF ANY OF THESE MODEL IS SELECTED,SET cb_extract_VMF_ZTDs ENABLE OFF
if all([strncmpi(dryModel,'VMF gridded ZHD',15),strncmpi(wetModel,'VMF gridded ZWD',15)])
  
   set(handles.cb_extract_VMF_ZTDs,'value',0,'enable','off')

elseif any([all([~strncmpi(dryModel,'VMF gridded ZHD',15),strncmpi(wetModel,'VMF gridded ZWD',15)]),...
            all([strncmpi(dryModel,'VMF gridded ZHD',15),~strncmpi(wetModel,'VMF gridded ZWD',15)]),...
            all([~strncmpi(dryModel,'VMF gridded ZHD',15),~strncmpi(wetModel,'VMF gridded ZWD',15)])])

       set(handles.cb_extract_VMF_ZTDs,'enable','on')

end
   
%================================END OF TROPOcorrection_separate_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 

%                   --------------------------------
%(2.3)**************COMBINE MODELs(POPUP MENU(popm))
%                   --------------------------------
% --- EXECUTES ON SELECTION CHANGE in pop_lTropo_combine.
function TROPOmodels_combine_popm_Callback(hObject, handles)
%--------------------------------------------------------------------------
global tropoModel
getconts = cellstr(get(hObject,'String'));%get popup menu contents as cell array
tropoModel=getconts{get(hObject,'Value')};%get selected item from pop_lTropo_combine
Modelv=get(hObject,'Value');%get selected item value from pop_lTropo_combine

%GET SELECTED METEOROLOGICAL PARAMETER SOURCE
getconts_MET = cellstr(get(handles.pop_source_metpara,'String'));%get Tropo_combine popup menu contents as cell array
METpara = getconts_MET{get(handles.pop_source_metpara,'Value')};%get selected item from pop_lTropo_combine
metVAL=get(handles.pop_source_metpara,'Value');%GET SELECTED VALUE

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%*******FOR READING GLOBAL PRESSURE & TEMPERATURE(GPT) GRIDFILES
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%GET COMBINE TROPO MODEL & IT'S GRIDFILE(grid_c)
tropoModel1=getappdata(0,'Tropo_Cmodel');%combine topospheric model
grid_c1=getappdata(0,'grid_c');%Previous comnine gridfile
grid_res_c1=getappdata(0,'grid_res_c');%Previous comnine grid resolution

%GET STORED dryModel & IT'S GRIDFILE(grid_h) FROM THIS POPUP MENU 
dryModel1=getappdata(0,'dryModel');%previous hydrostatic model
grid_h1=getappdata(0,'grid_h');   %previous gridfile
grid_res_h1=getappdata(0,'grid_res_h');%previous grid resolution

%GET STORED wetModel & IT'S GRIDFILE(grid_w) FROM WET MODEL POPUP MENU 
wetModel1=getappdata(0,'wetModel');%previous wet model
grid_w1=getappdata(0,'grid_w');   %previous gridfile
grid_res_w1=getappdata(0,'grid_res_w');%previous grid resolution
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%--------------------------------------------------------------------------+

%IF ANY OF THESE MODEL IS SELECTED,SET SOME UICONTROLS VISIBLE/ENABLE ON
if any([strncmpi(tropoModel,'no model',8),strncmpi(tropoModel,'GPT2 (5° x 5°)',14),strncmpi(tropoModel,'GPT2w (1° x 1°)',14),strncmpi(tropoModel,'GPT2w (5° x 5°)',14),...
        strncmpi(tropoModel,'GPT3 (1° x 1°)',14),strncmpi(tropoModel,'GPT3 (5° x 5°)',14),strncmpi(tropoModel,'UNB3m',5),...
        strncmpi(tropoModel,'EGNOS',5),strncmpi(tropoModel,'MOPS',4),strncmpi(tropoModel,'VMF gridded ZTD',15),strncmpi(tropoModel,'GTrop [Sun et al 2019]',22),...
        Modelv==[1,9,10,11,12,13,14,15,16,17,18]])
    
    %******SET MET PARAMETER SOURCE UICONTROLS VISIBLE/ENABLE ON
                 
   if strcmpi(get(handles.pop_source_metpara,'enable'),'On')
      set(handles.pop_source_metpara,'enable','off')
      set(handles.text_source_metpara,'enable','off')
      
   else
       set(handles.pop_source_metpara,'enable','off')
       set(handles.text_source_metpara,'enable','off')
 
   end
   
   if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
           any(metVAL==[7,8])])
       
      if strcmpi(get(handles.pop_met_manual,'visible'),'On')
          set(handles.pop_met_manual,'visible','off')
          set(handles.text_metmanual_source,'visible','Off')
      end
      
   elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
           any(metVAL==[4,5])])
          
          if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'on'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'on')])
             set(handles.rb_grid_resolution_1,'Visible','off')
             set(handles.rb_grid_resolution_5,'Visible','off')
             set(handles.text_grid_resolution,'Visible','off')  
          end 
   end
   
else
    
    %SET UICONTROLS VISIBLE/ENABLE ON
    set(handles.pop_source_metpara,'enable','on')
    set(handles.text_source_metpara,'enable','on')
    
    if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
            any(metVAL==[7,8])])
        
       if strcmpi(get(handles.pop_met_manual,'visible'),'Off') 
          set(handles.pop_met_manual,'visible','On')
          set(handles.text_metmanual_source,'visible','On')
       end
     
    elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
           any(metVAL==[4,5])])
          
          if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'off'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'off')])
             set(handles.rb_grid_resolution_1,'Visible','on')
             set(handles.rb_grid_resolution_5,'Visible','on')
             set(handles.text_grid_resolution,'Visible','on')  
          end
          
    end
       
end


%IF 'VMF gridded ZTD' EITHER SELECTED OR NOT, DO THE FF:
if ~any([strncmpi(tropoModel,'VMF gridded ZTD',15),Modelv==17])
   set(handles.cb_extract_VMF_ZTDs,'enable','on')
else 
   set(handles.cb_extract_VMF_ZTDs,'value',0,'enable','off')
end


if Modelv == 1 %IF NO MODEL IS SELECTED
   set(hObject,'TooltipString',strjoin({'Solution without','Tropospheric correction'}))%change Tooltip any time user select Model type 

     %IF ANY OF THE CONVENTIONAL/TRADITIONAL MODEL IS SELECTED
elseif any(Modelv==[2,3,4,5,6,7,8])
        set(hObject,'TooltipString',strjoin({'Solution using',tropoModel,'model'}))%change Tooltip any time user select Model type
    
     %IF ANY OF THE BLIND(LATITUDE AND/OR LONGITUDE,TIME DEPENDENT) MODEL IS SELECTED    
elseif Modelv == 9 || strncmpi(tropoModel,'UNB3m',5)
       set(hObject,'TooltipString',strjoin({'Solution using','University of New Brunswick V3 modified(UNB3m)','model'}))%change Tooltip any time user select Model type 

elseif Modelv == 10 || strncmpi(tropoModel,'EGNOS',5)
       set(hObject,'TooltipString',strjoin({'Solution using','European Geo-stationary Navigation Overlay System(EGNOS)','model'}))%change Tooltip any time user select Model type     

elseif Modelv == 11 || strncmpi(tropoModel,'MOPS',4)
       set(hObject,'TooltipString',strjoin({'Solution using','Minimum Operational Performance Standards(MOPS)','model'}))%change Tooltip any time user select Model type 
       
elseif  Modelv == 17 || strncmpi(tropoModel,'VMF gridded ZTD',15)

        %******SET TOOL TIS STRING
        %STRINGS TO SET
        string1 = 'Solution using Vienna Mapping Function (VMF) gridded Zenith Total Delays.';
        string2 = 'Please make sure you have VMF grid files (VMF1 / VMF3) available in goGPS data folder.';
        string3 = 'Default folder from goGPS is : ''... / data / TropoGRIDS / VMF files / VMF grids''';

        %CONCATENATE STRINGS
        string = sprintf('%s\n%s\n%s',string1, string2,string3); 
        
        %SET TOOL TIP
        set(hObject,'TooltipString',string)
          
        %CHECK IF VMF GRID FILEs ARE AVAILABLE IN goGPS
        %Call for the 'SearchVMFgrids.m' fxn
        [VMF_grid_found,VMFgrids] = SearchVMFgrids(handles); 
                        
        %*********SAVE STRUCT GRID FILES
        setappdata(0,'VMFgrids',VMFgrids)
         
        %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
        setappdata(0,'VMF_grid_found',VMF_grid_found)
        
        if VMF_grid_found == 1 %IF VMF grid FILES FOUND,LET USER CHOOSE GRID TYPE & RESOLUTION
                               %I.E.['VMF1(2°x2.5°)','VMF3(1°x1°)','VMF3(5°x5°)']
                               
           [VMF_type,grid_res] = VMFgrid_type(); 
           
           %*********CHECK IF VMF_type,grid_res ARE BOTH EMPTY
           %NOTE:
           %     VMF_type,grid_res ARE ASSIGNED EMPTY WHEN USER CANCEL SELECTION
           %     OR CLOSE DIALOGUE BOX FOR THE SELECTION OF VMF GRID TYPE &
           %     RESOLUTION.IF IT HAPPENS SO, THE SELECTED TROPO DELAY MODEL,
           %     'VMF gridded ZTD' WOULD BE REPLACED WITH SAASTEMOIN MODEL
           %     OR THE POPUP MENU WILL BE SET TO SAASTAMOINEN AS DEFAULT
        
           if all([isempty(VMF_type),isempty(grid_res)])%IF USER CANCEL SELECTION/CLOSE DIALOGUE BOX
            
              %SET COMBINE TROPO MODEL POPUP MENU VALUE TO 2(I.E.SAASTAMOINEN)
              set(handles.pop_lTropo_combine,'Value',2)
              
              %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
              set(handles.cb_extract_VMF_ZTDs,'enable','on','value',0)
                            
              %******SET MET PARAMETER SOURCE UICONTROLS VISIBLE/ENABLE ON
              set(handles.pop_source_metpara,'enable','on')
              set(handles.text_source_metpara,'enable','on')
    
              if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
                              any(metVAL==[7,8])])
        
                 if strcmpi(get(handles.pop_met_manual,'visible'),'Off') 
                    set(handles.pop_met_manual,'visible','On')
                    set(handles.text_metmanual_source,'visible','On')
                 end    
     
              elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
                                  any(metVAL==[4,5])])
          
                     if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'off'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'off')])
                        set(handles.rb_grid_resolution_1,'Visible','on')
                        set(handles.rb_grid_resolution_5,'Visible','on')
                        set(handles.text_grid_resolution,'Visible','on')  
                     end    
          
              end    %//if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
                        %          any(metVAL==[7,8])])  
            
           end   %//if all([isempty(VMF_type),isempty(grid_res)])
        
           %SAVE GRID FILE VERSION & TYPE
           setappdata(0,'VMFgrid_type',VMF_type)
           setappdata(0,'VMFgrid_res',grid_res)
           
        end %//if VMF_grid_found == 1
        
       
elseif   any([strncmpi(tropoModel,'GTrop [Sun et al 2019]',22),Modelv == 18])
    
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'Solution using','Global Tropospheric model(GTrop)'}))%change Tooltip any time user select Model type 
       
       %GET & SAVE GTrop GRID VALUES
        if ~isempty(getappdata(0,'gridV_GTrop'))
            
           setappdata(0,'gridV_GTrop',getappdata(0,'gridV_GTrop')) 
           
        else
            %SEARCH FOR GTrop model GRID FILE(GTropCoefficient.mat) & DIRECTORY & LOAD VALUES               
            %Call the "SearchGTropgrid.m" fxn
            [~,gridV_GTrop] = SearchGTropgrid();%****LOADED GTrop MODEL COEFFICIENTS
       
            %SAVE GRID VALUES
            setappdata(0,'gridV_GTrop',gridV_GTrop)
        end       
        
        
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%***********GLOBAL PRESSURE  AND TEMPERATURE(GPT) GRID MODELS
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
elseif Modelv == 12 || strncmpi(tropoModel,'GPT2 (5° x 5°)',14)
       set(hObject,'TooltipString',strjoin({'Solution using','Global Pressure and Temperature v2 model', 'based on a 5° external grid file'}))%change Tooltip any time user select Model type 
          
       %GET GPT MODELs GRID FILE & RESOLUTION
       grid_res_c=5; %SET GRID RESOLUTION TO 5
       
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 
       %CHECK IF MODEL HAS BEEN RUN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       if ~isempty(grid_c1) & grid_res_c1 == 5
           
           if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(tropoModel,tropoModel1,14)
              grid_c = grid_c1;
              
           else 
                %****TRY FROM OTHER MODELS
                %1.WET MODELs
                if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(tropoModel,wetModel1,14)
                   if ~isempty(grid_w1) & grid_res_w1 == 5
                     grid_c = grid_w1;
                   else
                       grid_c=readGPTgrid('gpt2_5.mat','GPT2',5);%Combine model GRIDFILE
                   end
                %2.HYDROSTATIC MODELs
                elseif any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(tropoModel,dryModel1,14)
                       if ~isempty(grid_h1) & grid_res_h1 == 5
                          grid_c = grid_h1;
                          
                       else
                           grid_c=readGPTgrid('gpt2_5.mat','GPT2',5);%Combine model GRIDFILE
                       end
                       
                else %IF ALL FAIL,READ GRID
                     grid_c=readGPTgrid('gpt2_5.mat','GPT2',5);%Combine model GRIDFILE
                     
                end
           end
           
       else %IF COMBINE(DRY+WET)MODEL GRIDFILE IS EMPTY([])
           if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(tropoModel,wetModel1,14)
              if ~isempty(grid_w1) & grid_res_w1 == 5
                 grid_c = grid_w1;
                 
              else
                  grid_c=readGPTgrid('gpt2_5.mat','GPT2',5);%Combine model GRIDFILE
              end 
           
           elseif any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(tropoModel,dryModel1,14)
              if ~isempty(grid_h1) & grid_res_h1 == 5
                 grid_c = grid_h1;
              else
                  grid_c=readGPTgrid('gpt2_5.mat','GPT2',5);%Combine model GRIDFILE
              end 
           else 
               grid_c=readGPTgrid('gpt2_5.mat','GPT2',5);%Combine model GRIDFILE
           end  
       end
       
       %IF ALL OPTIONS FAIL
       if ~exist('grid_c','var')
          grid_c=readGPTgrid('gpt2_5.mat','GPT2',5);%Combine model GRIDFILE
       end
       
       if length(grid_c) > 8 
          grid_c=readGPTgrid('gpt2_5.mat','GPT2',5);%Combine model GRIDFILE
       end     
  
%SAVE GRIDFILE & RESOLUTION       
setappdata(0,'grid_c',grid_c)
setappdata(0,'grid_res_c',grid_res_c)      
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 

elseif Modelv == 13 || strncmpi(tropoModel,'GPT2w (1° x 1°)',14)
       set(hObject,'TooltipString',strjoin({'Solution using','Global Pressure and Temperature to wet model','based on a 1° external grid file'}))%change Tooltip any time user select Model type 
       
       %GET GPT MODELs GRID FILE & RESOLUTION
       grid_res_c=1; %SET GRID RESOLUTION TO 1
       
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       if ~isempty(grid_c1) & grid_res_c1 == 1
           
           if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(tropoModel,tropoModel1,14)
              grid_c = grid_c1;
              
           else 
                %****TRY FROM OTHER MODELS
                %1.WET MODELs
                if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(tropoModel,wetModel1,14)
                   if ~isempty(grid_w1) & grid_res_w1 == 1
                     grid_c = grid_w1;
                   else
                        grid_c=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Combine model GRIDFILE
                   end
                %2.HYDROSTATIC MODELs
                elseif any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(tropoModel,dryModel1,14)
                       if ~isempty(grid_h1) & grid_res_h1 == 1
                          grid_c = grid_h1;
                       else
                            grid_c=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Combine model GRIDFILE
                       end
                       
                else %IF ALL FAIL,READ GRID
                     grid_c=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Combine model GRIDFILE
                     
                end
           end
           
       else %IF COMBINE(DRY+WET)MODEL GRIDFILE IS EMPTY([])
           if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(tropoModel,wetModel1,14)
              if ~isempty(grid_w1) & grid_res_w1 == 1
                 grid_c = grid_w1;
              else
                   grid_c=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Combine model GRIDFILE
              end 
           
           elseif any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(tropoModel,dryModel1,14)
              if ~isempty(grid_h1) & grid_res_h1 == 1
                 grid_c = grid_h1;
              else
                   grid_c=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Combine model GRIDFILE
              end 
           else 
               grid_c=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Combine model GRIDFILE
           end  
       end
       
       %IF ALL OPTIONS FAIL
       if ~exist('grid_c','var')
          grid_c=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Combine model GRIDFILE
       end
           
       if any([length(grid_c) > 10,length(grid_c)< 10])
          grid_c=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Combine model GRIDFILE
       end    
       
%SAVE GRIDFILE & RESOLUTION       
setappdata(0,'grid_c',grid_c)
setappdata(0,'grid_res_c',grid_res_c)        
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

elseif Modelv == 14 || strncmpi(tropoModel,'GTP2w (5° x 5° )',14)
       set(hObject,'TooltipString',strjoin({'Solution using','Global Pressure and Temperature to wet model', 'based on a 5° external grid file'}))%change Tooltip any time user select Model type 
       
       %GET GPT MODELs GRID FILE & RESOLUTION
       grid_res_c = 5; %SET GRID RESOLUTION TO 5
       
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       if ~isempty(grid_c1) & grid_res_c1 == 5
           
           if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(tropoModel,tropoModel1,14)
              grid_c = grid_c1;
              
           else 
                %****TRY FROM OTHER MODELS
                %1.WET MODELs
                if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(tropoModel,wetModel1,14)
                   if ~isempty(grid_w1) & grid_res_w1 == 5
                      grid_c = grid_w1;
                   else
                       grid_c=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Combine model GRIDFILE
                   end
                %2.HYDROSTATIC MODELs
                elseif any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(tropoModel,dryModel1,14)
                       if ~isempty(grid_h1) & grid_res_h1 == 5
                          grid_c = grid_h1;
                       else
                           grid_c=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Combine model GRIDFILE
                       end
                       
                else %IF ALL FAIL,READ GRID
                     grid_c=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Combine model GRIDFILE
                     
                end
           end
           
       else %IF COMBINE(DRY+WET)MODEL GRIDFILE IS EMPTY([])
           if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(tropoModel,wetModel1,14)
              if ~isempty(grid_w1) & grid_res_w1 == 5
                 grid_c = grid_w1;
              else
                  grid_c=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Combine model GRIDFILE
              end 
           
           elseif any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(tropoModel,dryModel1,14)
              if ~isempty(grid_h1) & grid_res_h1 == 5
                 grid_c = grid_h1;
              else
                  grid_c=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Combine model GRIDFILE
              end
              
           else 
               grid_c=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Combine model GRIDFILE
               
           end  
       end
       
      %IF ALL OPTIONS FAIL
       if ~exist('grid_c','var')
          grid_c=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Combine model GRIDFILE
       end
           
       if any([length(grid_c) > 10,length(grid_c)< 10])
          grid_c=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Combine model GRIDFILE
       end    
       
%SAVE GRIDFILE & RESOLUTION       
setappdata(0,'grid_c',grid_c)
setappdata(0,'grid_res_c',grid_res_c)         
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+  

elseif  Modelv == 15 || strncmpi(tropoModel,'GPT3 (1° x 1°)',14)
       set(hObject,'TooltipString',strjoin({'Solution using','Global Pressure and Temperature v3 model', 'based on a 1° external grid file'}))%change Tooltip any time user select Model type 
       
       %GET GPT MODELs GRID FILE & RESOLUTION
       grid_res_c=1; %SET GRID RESOLUTION TO 1
       
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       if ~isempty(grid_c1) & grid_res_c1 == 1
           
           if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(tropoModel,tropoModel1,14)
              grid_c = grid_c1;
              
           else 
                %****TRY FROM OTHER MODELS
                %1.WET MODELs
                if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(tropoModel,wetModel1,14)
                   if ~isempty(grid_w1) & grid_res_w1 == 1
                      grid_c = grid_w1;
                   else
                       grid_c=readGPTgrid('gpt3_1.mat','GPT3',1);%Combine model GRIDFILE
                   end
                %2.HYDROSTATIC MODELs
                elseif any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(tropoModel,dryModel1,14)
                       if ~isempty(grid_h1) & grid_res_h1 == 1
                          grid_c = grid_h1;
                       else
                           grid_c=readGPTgrid('gpt3_1.mat','GPT3',1);%Combine model GRIDFILE
                       end
                       
                else %IF ALL FAIL,READ GRID
                     grid_c=readGPTgrid('gpt3_1.mat','GPT3',1);%Combine model GRIDFILE
                     
                end
           end
           
       else %IF COMBINE(DRY+WET)MODEL GRIDFILE IS EMPTY([])
           if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(tropoModel,wetModel1,14)
              if ~isempty(grid_w1) & grid_res_w1 == 1
                 grid_c = grid_w1;
              else
                  grid_c=readGPTgrid('gpt3_1.mat','GPT3',1);%Combine model GRIDFILE
              end 
           
           elseif any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(tropoModel,dryModel1,14)
              if ~isempty(grid_h1) & grid_res_h1 == 1
                 grid_c = grid_h1;
              else
                  grid_c=readGPTgrid('gpt3_1.mat','GPT3',1);%Combine model GRIDFILE
              end 
              
           else 
               grid_c=readGPTgrid('gpt3_1.mat','GPT3',1);%Combine model GRIDFILE
               
           end  
       end
       
       %IF ALL OPTIONS FAIL
       if ~exist('grid_c','var')
          grid_c=readGPTgrid('gpt3_1.mat','GPT3',1);%Combine model GRIDFILE
       end
  
       if length(grid_c)<14
          grid_c=readGPTgrid('gpt3_1.mat','GPT3',1);%Combine model GRIDFILE
       end       
       
%SAVE GRIDFILE & RESOLUTION       
setappdata(0,'grid_c',grid_c)
setappdata(0,'grid_res_c',grid_res_c)         
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 

elseif Modelv == 16 || strncmpi(tropoModel,'GPT3 (5° x 5°)',14)
       set(hObject,'TooltipString',strjoin({'Solution using','Global Pressure and Temperature v3 model', 'based on a 5° external grid file'}))%change Tooltip any time user select Model type 
       
       %GET GPT MODELs GRID FILE & RESOLUTION
       grid_res_c = 5; %SET GRID RESOLUTION TO 5
       
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       if ~isempty(grid_c1) & grid_res_c1 == 5
           
           if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(tropoModel,tropoModel1,14)
              grid_c = grid_c1;
              
           else 
                %****TRY FROM OTHER MODELS
                %1.WET MODELs
                if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(tropoModel,wetModel1,14)
                   if ~isempty(grid_w1) & grid_res_w1 == 5
                     grid_c = grid_w1;
                   else
                       grid_c=readGPTgrid('gpt3_5.mat','GPT3',5);%Combine model GRIDFILE
                   end
                %2.HYDROSTATIC MODELs
                elseif any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(tropoModel,dryModel1,14)
                       if ~isempty(grid_h1) & grid_res_h1 == 5
                          grid_c = grid_h1;
                       else
                           grid_c=readGPTgrid('gpt3_5.mat','GPT3',5);%Combine model GRIDFILE
                       end
                       
                else %IF ALL FAIL,READ GRID
                     grid_c=readGPTgrid('gpt3_5.mat','GPT3',5);%Combine model GRIDFILE
                     
                end
           end
           
       else %IF COMBINE(DRY+WET)MODEL GRIDFILE IS EMPTY([])
           if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(tropoModel,wetModel1,14)
              if ~isempty(grid_w1) & grid_res_w1 == 5
                 grid_c = grid_w1;
              else
                  grid_c=readGPTgrid('gpt3_5.mat','GPT3',5);%Combine model GRIDFILE
              end 
           
           elseif any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(tropoModel,dryModel1,14)
              if ~isempty(grid_h1) & grid_res_h1 == 5
                 grid_c = grid_h1;
              else
                  grid_c=readGPTgrid('gpt3_5.mat','GPT3',5);%Combine model GRIDFILE
              end 
              
           else 
               grid_c=readGPTgrid('gpt3_5.mat','GPT3',5);%Combine model GRIDFILE
               
           end  
       end
       
       %IF ALL OPTIONS FAIL
       if ~exist('grid_c','var')
          grid_c=readGPTgrid('gpt3_5.mat','GPT3',5);%Combine model GRIDFILE
       end
       
       if length(grid_c)<14
          grid_c=readGPTgrid('gpt3_5.mat','GPT3',5);%Combine model GRIDFILE
       end
      
%SAVE GRIDFILE & RESOLUTION       
setappdata(0,'grid_c',grid_c)
setappdata(0,'grid_res_c',grid_res_c)  
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

end

%SAVE SELECTED COMBINE TROPO MODEL
setappdata(0,'Tropo_Cmodel',tropoModel)

%CHECK IF ANY OF THE GPT MODEL(GPT2,GPT2w & GPT3) IS SELECTED,AND ALLOW
%USER TO INDICATE WHETHER TO APPLY TIME VARIATION OF NOT
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%NOTE:
%    GLOBAL PRESSURE & TEMPERATURE(GPT) MODELS SUCH AS GPT2,GPT2w & GPT3 
%    COMPUTES METEOROLOGICAL PARAMETERS EITHER IN STATIC MODE OR IN TIME
%    VARYING MODE:
%    case 1: no time variation but static quantities
%    case 0: with time variation (annual and semiannual terms)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
if any([strncmpi(tropoModel,'GPT2 (5° x 5°)',14),strncmpi(tropoModel,'GPT2w (1° x 1°)',14),strncmpi(tropoModel,'GPT2w (5° x 5°)',14),...
        strncmpi(tropoModel,'GPT3 (1° x 1°)',14),strncmpi(tropoModel,'GPT3 (5° x 5°)',14), Modelv==[12,13,14,15,16]]) 
    
if ~isempty(grid_c)
    
   TIMEvariation()
   uiwait(gcf)
   setappdata(0,'Timevar_c',getappdata(0,'Timevar'))
   
end
  
end


%--------------------------------------------------------------------------   
%NOTE:
%SAASTAMOINEN,MARINI,Askne & Nordius,UNB3m,EGNOS & MOPS MODELS REQUIRE +VE
%ORTHOMETRIC HEIGHT--------------------------------------------------------
%********GET goGPS DEFAULT GEOID MODEL (EGM2008)  

%CHECK SELECTED TROPO MODEL 

if any([strncmpi(tropoModel,'Saastamoinen',12),strncmpi(tropoModel,'Saastamoinen(Refined)',21),...
        strncmpi(tropoModel,'Askne & Nordius',15),strncmpi(tropoModel,'Marini',6),...
        strncmpi(tropoModel,'UNB3m',5),strncmpi(tropoModel,'EGNOS',5),strncmpi(tropoModel,'MOPS',4)])        

   if ~isempty(getappdata(0,'geoid'))
               
      %SAVE FOR FUTURE USE
      setappdata(0,'geoid',getappdata(0,'geoid')) 
              
   else  
       try %goGPS v0.5.2 beta1
          gs = Go_State.getInstance;
          geoid = gs.getRefGeoid(); 
      
       catch  %goGPS v0.4.3 
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
         
       end  %try 

    %SAVE FOR FUTURE USE
    setappdata(0,'geoid',geoid)
   
   end %//if ~isempty(getappdata(0,'geoid'))
   
end %//if any([strncmpi(tropoModel,'Saastamoinen',12),strncmpi(tropoModel,'Saastamoinen(Refined)',21),...
    %          strncmpi(tropoModel,'Askne & Nordius',15),strncmpi(tropoModel,'Marini',6),...
    %          strncmpi(tropoModel,'UNB3m',5),strncmpi(tropoModel,'EGNOS',5),strncmpi(tropoModel,'MOPS',4)])
     
%===============================END OF TROPOmodels_combine__popm_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 

%                     -----------------------------------
%(2.4)****************HYDROSTATIC MODELs(POPUP MENU(popm) 
%                     -----------------------------------
% --- EXECUTES ON SELECTION CHANGE in pop_ITropo_hydrostatic.
function TROPOmodels_hydrostatic_popm_Callback(hObject, handles)
%--------------------------------------------------------------------------
global dryModel 

%GET SELECTED HYDROSTATIC MODEL
getconts = cellstr(get(hObject,'String'));%get ITropo_hydrostatic popup menu contents as cell array
dryModel=getconts{get(hObject,'Value')};%get selected item from pop_ITropo_hydrostatic
Modelv=get(hObject,'Value');%get selected item value from pop_lTropo_hydrostatic

%GET SELECTED WET TROPO MODEL
getconts_wet = cellstr(get(handles.pop_ITropo_wet,'String'));%get Tropo_combine popup menu contents as cell array
wetModel = getconts_wet{get(handles.pop_ITropo_wet,'Value')};%get selected item from pop_lTropo_combine
Modelv_W=get(handles.pop_ITropo_wet,'Value');%get selected item value from pop_lTropo_wet

%GET SELECTED METEOROLOGICAL PARAMETER SOURCE
getconts_MET = cellstr(get(handles.pop_source_metpara,'String'));%get MET popup menu contents as cell array
METpara = getconts_MET{get(handles.pop_source_metpara,'Value')};%get selected item from MET POPUP MENU
metVAL=get(handles.pop_source_metpara,'Value');%GET SELECTED VALUE

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%*******FOR READING GLOBAL PRESSURE & TEMPERATURE(GPT) GRIDFILES
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%GET COMBINE TROPO MODEL & IT'S GRIDFILE(grid_c)
tropoModel1=getappdata(0,'Tropo_Cmodel');%combine topospheric model
grid_c=getappdata(0,'grid_c');%Previous comnine gridfile
grid_res_c=getappdata(0,'grid_res_c');%Previous comnine grid resolution

%GET STORED dryModel & IT'S GRIDFILE(grid_h) FROM THIS POPUP MENU 
dryModel1=getappdata(0,'dryModel');%previous hydrostatic model
grid_h1=getappdata(0,'grid_h');   %previous gridfile
grid_res_h1=getappdata(0,'grid_res_h');%previous grid resolution

%GET STORED wetModel & IT'S GRIDFILE(grid_w) FROM WET MODEL POPUP MENU 
wetModel1=getappdata(0,'wetModel');%previous wet model
grid_w1=getappdata(0,'grid_w');   %previous gridfile
grid_res_w1=getappdata(0,'grid_res_w');%previous grid resolution
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%--------------------------------------------------------------------------+

%IF ANY OF THESE MODEL IS SELECTED,SET SOME UICONTROLS VISIBLE/ENABLE ON
if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
        strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),...
        strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),strncmpi(dryModel,'VMF gridded ZHD',15),strncmpi(dryModel,'GTrop [Sun et al 2019]',22)]) & ...   
    any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
        strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),...
        strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4),strncmpi(wetModel,'VMF gridded ZWD',15),strncmpi(wetModel,'GTrop [Sun et al 2019]',22)])
    
   %******SET MET PARAMETER SOURCE UICONTROLS VISIBLE/ENABLE ON
   
   if strcmpi(get(handles.pop_source_metpara,'enable'),'On')
      set(handles.pop_source_metpara,'enable','off')
      set(handles.text_source_metpara,'enable','off')
 
   end
   
   if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
           any(metVAL==[7,8])])
       
      if strcmpi(get(handles.pop_met_manual,'visible'),'On')
          set(handles.pop_met_manual,'visible','off')
          set(handles.text_metmanual_source,'visible','Off')
      end
      
   elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
           any(metVAL==[4,5])])
          
          if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'on'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'on')])
             set(handles.rb_grid_resolution_1,'Visible','off')
             set(handles.rb_grid_resolution_5,'Visible','off')
             set(handles.text_grid_resolution,'Visible','off')  
          end 
   end
   
else
    
    %SET UICONTROLS VISIBLE/ENABLE ON
    set(handles.pop_source_metpara,'enable','on')
    set(handles.text_source_metpara,'enable','on')
    
    if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
            any(metVAL==[7,8])])
        
       if strcmpi(get(handles.pop_met_manual,'visible'),'Off') 
          set(handles.pop_met_manual,'visible','On')
          set(handles.text_metmanual_source,'visible','On')
       end
     
    elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
           any(metVAL==[4,5])])
          
          if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'off'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'off')])
             set(handles.rb_grid_resolution_1,'Visible','on')
             set(handles.rb_grid_resolution_5,'Visible','on')
             set(handles.text_grid_resolution,'Visible','on')  
          end
          
    end
       
end


%IF 'VMF gridded ZHD' EITHER SELECTED OR NOT, DO THE FF:
if ~any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv==15])
   set(handles.cb_extract_VMF_ZTDs,'enable','on')
else 
   set(handles.cb_extract_VMF_ZTDs,'value',0,'enable','off')
end


%IF ANY OF THE CONVENTIONAL/TRADITIONAL MODEL IS SELECTED
if any(Modelv==[1,2,3,4,5,6])  
   set(hObject,'TooltipString',strjoin({'Solution using',dryModel,'model'}))%change Tooltip any time user select Model type

   %IF ANY OF THE BLIND(LATITUDE AND/OR LONGITUDE,TIME DEPENDENT) MODEL IS SELECTED  
elseif Modelv == 7 || strncmpi(dryModel,'UNB3m',5)
       set(hObject,'TooltipString',strjoin({'Solution using','University of New Brunswick V3 modified(UNB3m)','model'}))%change Tooltip any time user select Model type 

elseif Modelv == 8 || strncmpi(dryModel,'EGNOS',5)
       set(hObject,'TooltipString',strjoin({'Solution using','European Geo-stationary Navigation Overlay System(EGNOS)','model'}))%change Tooltip any time user select Model type     

elseif Modelv == 9 || strncmpi(dryModel,'MOPS',4)
       set(hObject,'TooltipString',strjoin({'Solution using','Minimum Operational Performance Standards(MOPS)','model'}))%change Tooltip any time user select Model type 
       
elseif  Modelv == 15 || strncmpi(dryModel,'VMF gridded ZHD',15)
        
        %******SET TOOL TIS STRING
        %STRINGS TO SET
        string1 = 'Solution using Vienna Mapping Functions(VMF) gridded Zenith Hydrostatic Delays.';
        string2 = 'Please make sure you have VMF grid files (VMF1 / VMF3) available in goGPS data folder.';
        string3 = 'Default folder from goGPS is : ''... / data / TropoGRIDS / VMF files / VMF grids''';

        %CONCATENATE STRINGS
        string = sprintf('%s\n%s\n%s',string1, string2,string3); 
        
        %SET TOOL TIP
        set(hObject,'TooltipString',string)%change Tooltip any time user select Model type 
         
        %CHECK IF VMF GRID FILEs ARE AVAILABLE IN goGPS
        %Call for the 'SearchVMFgrids.m' fxn
        [VMF_grid_found,VMFgrids] = SearchVMFgrids(handles) ;
                        
        %*********SAVE STRUCT GRID FILES
        setappdata(0,'VMFgrids',VMFgrids)
         
        %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
        setappdata(0,'VMF_grid_found',VMF_grid_found)
        
        if VMF_grid_found == 1 %IF VMF grid FILES FOUND,LET USER CHOOSE GRID TYPE & RESOLUTION
                               %I.E.['VMF1(2°x2.5°)','VMF3(1°x1°)','VMF3(5°x5°)']
                               
           [VMF_type,grid_res] = VMFgrid_type(); 
           
           %*********CHECK IF VMF_type,grid_res ARE BOTH EMPTY
           %NOTE:
           %     VMF_type,grid_res ARE ASSIGNED EMPTY WHEN USER CANCEL SELECTION
           %     OR CLOSE DIALOGUE BOX FOR THE SELECTION OF VMF GRID TYPE &
           %     RESOLUTION.IF IT HAPPENS SO, THE SELECTED TROPO DELAY MODEL,
           %     'VMF gridded ZTD' WOULD BE REPLACED WITH SAASTEMOIN MODEL
           %     OR THE POPUP MENU WILL BE SET TO SAASTAMOINEN AS DEFAULT
        
           if all([isempty(VMF_type),isempty(grid_res)])%IF USER CANCEL SELECTION/CLOSE DIALOGUE BOX
        
              %SET DRY TROPO MODEL POPUP MENU VALUE TO 2(I.E.SAASTAMOINEN)
              set(handles.pop_ITropo_hydrostatic,'Value',1)
              
              if any([any([any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv == 15]),...
                           any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv_W == 20])]),...
                      all([any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv == 15]),...
                           any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv_W == 20])])])
                
                  %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
                  set(handles.cb_extract_VMF_ZTDs,'enable','on','value',0)
                  
              end
              
              %******SET MET PARAMETER SOURCE UICONTROLS VISIBLE/ENABLE ON
              set(handles.pop_source_metpara,'enable','on')
              set(handles.text_source_metpara,'enable','on')
    
              if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
                              any(metVAL==[7,8])])
        
                 if strcmpi(get(handles.pop_met_manual,'visible'),'Off') 
                    set(handles.pop_met_manual,'visible','On')
                    set(handles.text_metmanual_source,'visible','On')
                 end     
     
              elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
                                  any(metVAL==[4,5])])
          
                     if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'off'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'off')])
                        set(handles.rb_grid_resolution_1,'Visible','on')
                        set(handles.rb_grid_resolution_5,'Visible','on')
                        set(handles.text_grid_resolution,'Visible','on')  
                     end      
          
              end  %//if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
                   %                  any(metVAL==[7,8])])
                   
           else 
               
               if any([all([~any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv == 15]),...
                            any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv_W == 20])]),...
                       all([any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv == 15]),...
                            ~any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv_W == 20])])])
                
                  %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
                  set(handles.cb_extract_VMF_ZTDs,'enable','on','value',0)
                  
               elseif all([any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv == 15]),...
                           any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv_W == 20])])
                       
                      %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
                      set(handles.cb_extract_VMF_ZTDs,'enable','off','value',0)
                  
              end       
                        
           end  %//if all([isempty(VMF_type),isempty(grid_res)])
        
           %SAVE GRID FILE VERSION & TYPE
           setappdata(0,'VMFgrid_type',VMF_type)
           setappdata(0,'VMFgrid_res',grid_res)
           
        end  %//if VMF_grid_found == 1
           
elseif  any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),Modelv == 16])
        
        %SET TOOLTIP STRING
        set(hObject,'TooltipString','Solution using Global Tropospheric model(GTrop)')%change Tooltip any time user select Model type 
       
        %GET & SAVE GTrop GRID VALUES
        if ~isempty(getappdata(0,'gridV_GTrop'))
            
           setappdata(0,'gridV_GTrop',getappdata(0,'gridV_GTrop')) 
           
        else
            %SEARCH FOR GTrop model GRID FILE(GTropCoefficient.mat) & DIRECTORY & LOAD VALUES               
            %Call the "SearchGTropgrid.m" fxn
            [~,gridV_GTrop] = SearchGTropgrid();%****LOADED GTrop MODEL COEFFICIENTS
       
            %SAVE GRID VALUES
            setappdata(0,'gridV_GTrop',gridV_GTrop)
        end
        
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%***********GLOBAL PRESSURE  AND TEMPERATURE(GPT) GRID MODELS
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
elseif Modelv == 10 || strncmpi(dryModel,'GPT2 (5° x 5°)',14)
       set(hObject,'TooltipString',strjoin({'Solution using','Global Pressure and Temperature v2 model', 'based on a (5° x 5°) external grid file'}))%change Tooltip any time user select Model type 
       
       %GET GPT MODELs GRID FILE & RESOLUTION
       grid_res_h=5; %SET GRID RESOLUTION TO 5
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       if ~isempty(grid_h1) & grid_res_h1 == 5
           
           if any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(dryModel,dryModel1,14)
              grid_h = grid_h1;
              
           else 
                %****TRY FROM OTHER MODELS
                %1.WET MODELs
                if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(dryModel,wetModel1,14)
                   if ~isempty(grid_w1) & grid_res_w1 == 5
                     grid_h = grid_w1;
                   else
                       grid_h=readGPTgrid('gpt2_5.mat','GPT2',5);%HYDROSTATIC GRIDFILE
                   end
                %2.COMBINE MODELs
                elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(dryModel,tropoModel1,14)
                       if ~isempty(grid_c) & grid_res_c == 5
                          grid_h = grid_c;
                       else
                           grid_h=readGPTgrid('gpt2_5.mat','GPT2',5);%HYDROSTATIC GRIDFILE
                       end
                       
                else %IF ALL FAIL,READ GRID
                     grid_h=readGPTgrid('gpt2_5.mat','GPT2',5);%HYDROSTATIC GRIDFILE
                     
                end
           end
           
       else %IF HYDROSTATIC GRIDFILE IS EMPTY([])
           if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(dryModel,wetModel1,14)
              if ~isempty(grid_w1) & grid_res_w1 == 5
                 grid_h = grid_w1;
                 
              else
                  grid_h=readGPTgrid('gpt2_5.mat','GPT2',5);%HYDROSTATIC GRIDFILE
                  
              end 
           
           elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(dryModel,tropoModel1,14)
              if ~isempty(grid_c) & grid_res_c == 5
                 grid_h = grid_c;
              else 
                  grid_h=readGPTgrid('gpt2_5.mat','GPT2',5);%HYDROSTATIC GRIDFILE
              end 
           else 
               grid_h=readGPTgrid('gpt2_5.mat','GPT2',5);%HYDROSTATIC GRIDFILE
           end  
       end
       
       %IF ALL OPTIONS FAIL
       if ~exist('grid_h','var')
          grid_h=readGPTgrid('gpt2_5.mat','GPT2',5);%HYDROSTATIC GRIDFILE
       end
       
      if length(grid_h) > 8 
         grid_h=readGPTgrid('gpt2_5.mat','GPT2',5);%HYDROSTATIC model GRIDFILE
      end    
%SAVE HYDROSTATIC MODEL GRID FILE(grid_h) & RESOLUTION(grid_res_h)
setappdata(0,'grid_h',grid_h)
setappdata(0,'grid_res_h',grid_res_h)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

elseif Modelv == 11 || strncmpi(dryModel,'GPT2w (1° x 1°)',14)
       set(hObject,'TooltipString',strjoin({'Solution using','Global Pressure and Temperature to wet model', 'based on a (1° x 1°) external grid file'}))%change Tooltip any time user select Model type 
    
       %GET GPT MODELs GRID FILE & RESOLUTION
       grid_res_h=1; %SET GRID RESOLUTION TO 1
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       if ~isempty(grid_h1) & grid_res_h1 == 1
           
           if any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(dryModel,dryModel1,14)
              grid_h = grid_h1;
              
           else 
                %****TRY FROM OTHER MODELS
                %1.WET MODELs
                if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(dryModel,wetModel1,14)
                   if ~isempty(grid_w1) & grid_res_w1 == 1
                     grid_h = grid_w1;
                   else
                       grid_h=readGPTgrid('gpt2_1w.mat','GPT2w',1);%HYDROSTATIC GRIDFILE
                   end
                %2.COMBINE MODELs
                elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(dryModel,tropoModel1,14)
                       if ~isempty(grid_c) & grid_res_c == 1
                          grid_h = grid_c;
                       else
                           grid_h=readGPTgrid('gpt2_1w.mat','GPT2w',1);%HYDROSTATIC GRIDFILE
                       end
                       
                else %IF ALL FAIL,READ GRID
                     grid_h=readGPTgrid('gpt2_1w.mat','GPT2w',1);%HYDROSTATIC GRIDFILE
                     
                end
           end
           
       else %IF HYDROSTATIC GRIDFILE IS EMPTY([])
           
           if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(dryModel,wetModel1,14)
              if ~isempty(grid_w1) & grid_res_w1 == 1
                 grid_h = grid_w1;
              else
                  grid_h=readGPTgrid('gpt2_1w.mat','GPT2w',1);%HYDROSTATIC GRIDFILE
              end 
           
           elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(dryModel,tropoModel1,14)
                  if ~isempty(grid_c) & grid_res_c == 1
                     grid_h = grid_c;
                  else
                      grid_h=readGPTgrid('gpt2_1w.mat','GPT2w',1);%HYDROSTATIC GRIDFILE
                  end  
           else  
               grid_h=readGPTgrid('gpt2_1w.mat','GPT2w',1);%HYDROSTATIC GRIDFILE
           end   
       end
      
      %IF ALL OPTIONS FAIL
      if ~exist('grid_h','var')
         grid_h = readGPTgrid('gpt2_1w.mat','GPT2w',1);%HYDROSTATIC model GRIDFILE
      end   
       
     if any([length(grid_h) > 10,length(grid_h)< 10])
        grid_h = readGPTgrid('gpt2_1w.mat','GPT2w',1);%HYDROSTATIC model GRIDFILE
     end 
     
%SAVE HYDROSTATIC MODEL GRID FILE(grid_h) & RESOLUTION(grid_res_h)
setappdata(0,'grid_h',grid_h)
setappdata(0,'grid_res_h',grid_res_h)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

elseif Modelv == 12 || strncmpi(dryModel,'GTP2w (5° x 5° )',14)
       set(hObject,'TooltipString',strjoin({'Solution using','Global Pressure and Temperature to wet model', 'based on a (5° x 5°) external grid file'}))%change Tooltip any time user select Model type 
    
       %GET GPT MODELs GRID FILE & RESOLUTION
       grid_res_h=5; %SET GRID RESOLUTION TO 5
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       if ~isempty(grid_h1) & grid_res_h1 == 5
           
           if any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(dryModel,dryModel1,14)
              grid_h = grid_h1;
              
           else 
                %****TRY FROM OTHER MODELS
                %1.WET MODELs
                if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(dryModel,wetModel1,14)
                   if ~isempty(grid_w1) & grid_res_w1 == 5
                      grid_h = grid_w1;
                   else
                        grid_h=readGPTgrid('gpt2_5w.mat','GPT2w',5);%HYDROSTATIC GRIDFILE
                   end
                %2.COMBINE MODELs
                elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(dryModel,tropoModel1,14)
                       if ~isempty(grid_c) & grid_res_c == 5
                          grid_h = grid_c;
                       else
                           grid_h=readGPTgrid('gpt2_5w.mat','GPT2w',5);%HYDROSTATIC GRIDFILE
                       end
                       
                else %IF ALL FAIL,READ GRID
                     grid_h=readGPTgrid('gpt2_5w.mat','GPT2w',5);%HYDROSTATIC GRIDFILE
                     
                end
           end
           
       else %IF HYDROSTATIC GRIDFILE IS EMPTY([])
           
           if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(dryModel,wetModel1,14)
              if ~isempty(grid_w1) & grid_res_w1 == 5
                 grid_h = grid_w1;
              else
                  grid_h=readGPTgrid('gpt2_5w.mat','GPT2w',5);%HYDROSTATIC GRIDFILE
              end 
           
           elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(dryModel,tropoModel1,14)
                  if ~isempty(grid_c) & grid_res_c == 5
                     grid_h = grid_c;
                  else 
                      grid_h=readGPTgrid('gpt2_5w.mat','GPT2w',5);%HYDROSTATIC GRIDFILE
                  end   
           else  
               grid_h=readGPTgrid('gpt2_5w.mat','GPT2w',5);%HYDROSTATIC GRIDFILE
           end   
       end 
       
       %IF ALL OPTIONS FAIL
      if ~exist('grid_h','var')
         grid_h = readGPTgrid('gpt2_5w.mat','GPT2w',5);%HYDROSTATIC model GRIDFILE
      end   
       
     if any([length(grid_h) > 10,length(grid_h)< 10])
        grid_h = readGPTgrid('gpt2_5w.mat','GPT2w',5);%HYDROSTATIC model GRIDFILE
     end 
              
%SAVE HYDROSTATIC MODEL GRID FILE(grid_h) & RESOLUTION(grid_res_h)
setappdata(0,'grid_h',grid_h)
setappdata(0,'grid_res_h',grid_res_h)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  

elseif Modelv == 13 || strncmpi(dryModel,'GPT3 (1° x 1°)',14)
       set(hObject,'TooltipString',strjoin({'Solution using','Global Pressure and Temperature v3 model', 'based on a (1° x 1°) external grid file'}))%change Tooltip any time user select Model type  

       %GET GPT MODELs GRID FILE & RESOLUTION
       grid_res_h=1; %SET GRID RESOLUTION TO 1
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       if ~isempty(grid_h1) & grid_res_h1 == 1
           
           if any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(dryModel,dryModel1,14)
              grid_h = grid_h1;
              
           else 
                %****TRY FROM OTHER MODELS
                %1.WET MODELs
                if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(dryModel,wetModel1,14)
                   if ~isempty(grid_w1) & grid_res_w1 == 1
                     grid_h = grid_w1;
                   else
                       grid_h=readGPTgrid('gpt3_1.mat','GPT3',1);%HYDROSTATIC GRIDFILE
                   end
                %2.COMBINE MODELs
                elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(dryModel,tropoModel1,14)
                       if ~isempty(grid_c) & grid_res_c == 1
                          grid_h = grid_c;
                       else
                           grid_h=readGPTgrid('gpt3_1.mat','GPT3',1);%HYDROSTATIC GRIDFILE
                       end
                       
                else %IF ALL FAIL,READ GRID
                     grid_h=readGPTgrid('gpt3_1.mat','GPT3',1);%HYDROSTATIC GRIDFILE
                     
                end
           end
           
       else %IF HYDROSTATIC GRIDFILE IS EMPTY([])
           
           if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(dryModel,wetModel1,14)
              if ~isempty(grid_w1) & grid_res_w1 == 1
                 grid_h = grid_w1;
              else
                  grid_h=readGPTgrid('gpt3_1.mat','GPT3',1);%HYDROSTATIC GRIDFILE
              end 
           
           elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(dryModel,tropoModel1,14)
                  if ~isempty(grid_c) & grid_res_c == 1
                     grid_h = grid_c;
                  else 
                      grid_h=readGPTgrid('gpt3_1.mat','GPT3',1);%HYDROSTATIC GRIDFILE
                  end   
           else  
               grid_h=readGPTgrid('gpt3_1.mat','GPT3',1);%HYDROSTATIC GRIDFILE
           end   
       end 
       
       if ~exist('grid_h','var')
          grid_h=readGPTgrid('gpt3_1.mat','GPT3',1);%HYDROSTATIC GRIDFILE
       end
       
      if length(grid_h)<14
         grid_h=readGPTgrid('gpt3_1.mat','GPT3',1);%Combine model GRIDFILE
      end    
 
%SAVE HYDROSTATIC MODEL GRID FILE(grid_h) & RESOLUTION(grid_res_h)
setappdata(0,'grid_h',grid_h)
setappdata(0,'grid_res_h',grid_res_h)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

elseif Modelv == 14 || strncmpi(dryModel,'GPT3 (5° x 5°)',14)
       set(hObject,'TooltipString',strjoin({'Solution using','Global Pressure and Temperature v3 model', 'based on a (5° x 5°) external grid file'}))%change Tooltip any time user select Model type 
       
       %GET GPT MODELs GRID FILE & RESOLUTION
       grid_res_h=5; %SET GRID RESOLUTION TO 5
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       if ~isempty(grid_h1) & grid_res_h1 == 5
           
           if any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(dryModel,dryModel1,14)
              grid_h = grid_h1;
              
           else 
                %****TRY FROM OTHER MODELS
                %1.WET MODELs
                if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(dryModel,wetModel1,14)
                   if ~isempty(grid_w1) & grid_res_w1 == 5
                     grid_h = grid_w1;
                   else
                       grid_h=readGPTgrid('gpt3_5.mat','GPT3',5);%HYDROSTATIC GRIDFILE
                   end
                %2.COMBINE MODELs
                elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(dryModel,tropoModel1,14)
                       if ~isempty(grid_c) & grid_res_c == 5
                          grid_h = grid_c;
                       else
                           grid_h=readGPTgrid('gpt3_5.mat','GPT3',5);%HYDROSTATIC GRIDFILE
                       end
                       
                else %IF ALL FAIL,READ GRID
                     grid_h=readGPTgrid('gpt3_5.mat','GPT3',5);%HYDROSTATIC GRIDFILE
                     
                end
           end
           
       else %IF HYDROSTATIC GRIDFILE IS EMPTY([])
           
           if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(dryModel,wetModel1,14)
              if ~isempty(grid_w1) & grid_res_w1 == 5
                 grid_h = grid_w1;
              else
                  grid_h=readGPTgrid('gpt3_5.mat','GPT3',5);%HYDROSTATIC GRIDFILE
              end 
           
           elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(dryModel,tropoModel1,14)
                  if ~isempty(grid_c) & grid_res_c == 5
                     grid_h = grid_c;
                  else
                      grid_h=readGPTgrid('gpt3_5.mat','GPT3',5);%HYDROSTATIC GRIDFILE
                  end  
           else  
               grid_h=readGPTgrid('gpt3_5.mat','GPT3',5);%HYDROSTATIC GRIDFILE
           end   
       end 
       
       if ~exist('grid_h','var')
          grid_h=readGPTgrid('gpt3_5.mat','GPT3',5);%HYDROSTATIC GRIDFILE
       end
       
      if length(grid_h)<14
         grid_h=readGPTgrid('gpt3_5.mat','GPT3',5);%HYDROSTATIC model GRIDFILE
      end
       
%SAVE HYDROSTATIC MODEL GRID FILE(grid_h) & RESOLUTION(grid_res_h)
setappdata(0,'grid_h',grid_h)
setappdata(0,'grid_res_h',grid_res_h)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

end

%SAVE SELECTED HYDROSTATIC TROPO MODEL
setappdata(0,'dryModel',dryModel)


%CHECK IF ANY OF THE GPT MODEL(GPT2,GPT2w & GPT3) IS SELECTED,AND ALLOW
%USER TO INDICATE WHETHER TO APPLY TIME VARIATION OF NOT
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%NOTE:
%    GLOBAL PRESSURE & TEMPERATURE(GPT) MODELS SUCH AS GPT2,GPT2w & GPT3 
%    COMPUTES METEOROLOGICAL PARAMETERS EITHER IN STATIC MODE OR IN TIME
%    VARYING MODE:
%    case 1: no time variation but static quantities
%    case 0: with time variation (annual and semiannual terms)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
        strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]) & ...   
    any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
         strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])

if ~isempty(grid_h)
    
   TIMEvariation()
   uiwait(gcf)
  setappdata(0,'Timevar_dry',getappdata(0,'Timevar'))
  
end

elseif any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
            strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]) & ...   
       ~any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
             strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])

     if ~isempty(grid_h)
         
        TIMEvariation()
        uiwait(gcf)
        setappdata(0,'Timevar_dry',getappdata(0,'Timevar'))  
        
     end
end

%-------------------------------------------------------------------  
%NOTE"
%SAASTAMOINEN,DAVIS ET AL,Askne & Nordius,UNB3m,EGNOS & MOPS MODELS 
%REQUIRE +VE ORTHOMETRIC HEIGHT-------------------------------------
%********GET goGPS DEFAULT GEOID MODEL (EGM2008) 

%CHECK SELECTE  MODEL
if any([strncmpi(dryModel,'Saastamoinen',12),strncmpi(dryModel,'Davis et al)',12),strncmpi(dryModel,'Askne & Nordius',15),...
        strncmpi(dryModel,'UNB3m',5),strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4)])          
          
    %LOAD GEOID MODEL[%USING goGPS DEFAULT GEOID MODEL (EGM2008)] 
    if ~isempty(getappdata(0,'geoid'))
               
       %SAVE FOR FUTURE USE
       setappdata(0,'geoid',getappdata(0,'geoid')) 
              
    else           
         try %goGPS v0.5.2 beta1
             gs = Go_State.getInstance;
          geoid = gs.getRefGeoid(); 
      
         catch  %goGPS v0.4.3 
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

         %SAVE FOR FUTURE USE
         setappdata(0,'geoid',geoid)
             
    end %/if ~isempty(getappdata(0,'geoid'))    
       
end  %//if any([strncmpi(dryModel,'Saastamoinen',12),strncmpi(dryModel,'Davis et al)',12),strncmpi(dryModel,'Askne & Nordius',15),...
     %          strncmpi(dryModel,'UNB3m',5),strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4)])
          
%===========================END OF TROPOmodels_hydrostatic__popm_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 

%                    ----------------------------
%(2.5)************** WET MODELs(POPUP MENU(popm)) 
%                    ----------------------------
% --- EXECUTES ON SELECTION CHANGE in pop_ITropo_wet.
function TROPOmodels_wet_popm_Callback(hObject, handles)
%--------------------------------------------------------------------------
global wetModel

%GET SELECTED WET TROPO MODEL
getconts = cellstr(get(hObject,'String'));%get ITropo_wet popup menu contents as cell array
wetModel=getconts{get(hObject,'Value')};%get selected item from pop_ITropo_wet
Modelv=get(hObject,'Value');%get selected item value from pop_lTropo_hydrostatic

%GET SELECTED HYDROSTATIC MODEL
getconts_dry = cellstr(get(handles.pop_ITropo_hydrostatic,'String'));%get pop_ITropo_hydrostatic popup menu contents as cell array
dryModel = getconts_dry{get(handles.pop_ITropo_hydrostatic,'Value')};%get selected item from pop_ITropo_hydrostatic
Modelv_H=get(handles.pop_ITropo_hydrostatic,'Value');%get selected item value from pop_lTropo_hydrostatic

%GET SELECTED METEOROLOGICAL PARAMETER SOURCE
getconts_MET = cellstr(get(handles.pop_source_metpara,'String'));%get MET popup menu contents as cell array
METpara = getconts_MET{get(handles.pop_source_metpara,'Value')};%get selected item from MET POPUP MENU
metVAL=get(handles.pop_source_metpara,'Value');%GET SELECTED VALUE

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%*******FOR READING GLOBAL PRESSURE & TEMPERATURE(GPT) GRIDFILES
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%GET COMBINE TROPO MODEL & IT'S GRIDFILE(grid_c)
tropoModel1=getappdata(0,'Tropo_Cmodel');%combine topospheric model
grid_c=getappdata(0,'grid_c');%Previous comnine gridfile
grid_res_c=getappdata(0,'grid_res_c');%Previous comnine grid resolution

%GET STORED dryModel & IT'S GRIDFILE(grid_w) FROM THIS POPUP MENU 
wetModel1=getappdata(0,'wetModel');%previous hydrostatic model
grid_w1=getappdata(0,'grid_w');   %previous gridfile
grid_res_w1=getappdata(0,'grid_res_w');%previous grid resolution

%GET STORED dryModel & IT'S GRIDFILE(grid_h) FROM DRY MODEL POPUP MENU 
dryModel1=getappdata(0,'dryModel');%previous hydrostatic model
grid_h1=getappdata(0,'grid_h');   %previous gridfile
grid_res_h1=getappdata(0,'grid_res_h');%previous grid resolution
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%--------------------------------------------------------------------------+

%IF ANY OF THESE MODEL IS SELECTED,SET SOME UICONTROLS VISIBLE/ENABLE ON
if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
        strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),...
        strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),strncmpi(dryModel,'VMF gridded ZHD',7),strncmpi(dryModel,'GTrop [Sun et al 2019]',22)]) & ...   
    any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
         strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),...
         strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4),strncmpi(wetModel,'VMF gridded ZWD',7),strncmpi(wetModel,'GTrop [Sun et al 2019]',22)])
   
   %******SET MET PARAMETER SOURCE UICONTROLS VISIBLE/ENABLE ON 
   if strcmpi(get(handles.pop_source_metpara,'enable'),'On')
      set(handles.pop_source_metpara,'enable','off')
      set(handles.text_source_metpara,'enable','off')
 
   end
   
   if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
           any(metVAL==[7,8])])
       
      if strcmpi(get(handles.pop_met_manual,'visible'),'On')
          set(handles.pop_met_manual,'visible','off')
          set(handles.text_metmanual_source,'visible','Off')
      end
      
   elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
           any(metVAL==[4,5])])
          
          if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'on'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'on')])
             set(handles.rb_grid_resolution_1,'Visible','off')
             set(handles.rb_grid_resolution_5,'Visible','off')
             set(handles.text_grid_resolution,'Visible','off')  
          end 
   end
   
else
    
    %SET UICONTROLS VISIBLE/ENABLE ON
    set(handles.pop_source_metpara,'enable','on')
    set(handles.text_source_metpara,'enable','on')
    
    if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
            any(metVAL==[7,8])])
        
       if strcmpi(get(handles.pop_met_manual,'visible'),'Off') 
          set(handles.pop_met_manual,'visible','On')
          set(handles.text_metmanual_source,'visible','On')
       end
     
    elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
           any(metVAL==[4,5])])
          
          if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'off'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'off')])
             set(handles.rb_grid_resolution_1,'Visible','on')
             set(handles.rb_grid_resolution_5,'Visible','on')
             set(handles.text_grid_resolution,'Visible','on')  
          end
          
    end
       
end

%IF 'VMF gridded ZHD' EITHER SELECTED OR NOT, DO THE FF:
if ~any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv==20])
   set(handles.cb_extract_VMF_ZTDs,'enable','on')
else 
   set(handles.cb_extract_VMF_ZTDs,'value',0,'enable','off')
end


%IF ANY OF THE CONVENTIONAL/TRADITIONAL MODEL IS SELECTED
if any(Modelv==[1,2,3,4,5,6,7,8,9,10,11])  
   set(hObject,'TooltipString',strjoin({'Solution using',wetModel,'model'}))%change Tooltip any time user select Model type

    %IF ANY OF THE BLIND(LATITUDE AND/OR LONGITUDE,TIME DEPENDENT) MODEL IS SELECTED
elseif Modelv == 12 || strncmpi(wetModel,'UNB3m',5)
       set(hObject,'TooltipString',strjoin({'Solution using','University of New Brunswick V3 modified(UNB3m)','model'}))%change Tooltip any time user select Model type 

elseif Modelv == 13 || strncmpi(wetModel,'EGNOS',5)
       set(hObject,'TooltipString',strjoin({'Solution using','European Geo-stationary Navigation Overlay System(EGNOS)','model'}))%change Tooltip any time user select Model type     

elseif Modelv == 14 || strncmpi(wetModel,'MOPS',4)
       set(hObject,'TooltipString',strjoin({'Solution using','Minimum Operational Performance Standards(MOPS)','model'}))%change Tooltip any time user select Model type 
       
elseif Modelv == 20 || strncmpi(wetModel,'VMF gridded ZWD',15)
       %******SET TOOL TIS STRING
        %STRINGS TO SET
        string1 = 'Solution using Vienna Mapping Functions(VMF) gridded Zenith Wet Delays.';
        string2 = 'Please make sure you have VMF grid files (VMF1 / VMF3) available in goGPS data folder.';
        string3 = 'Default folder from goGPS is : ''... / data / TropoGRIDS / VMF files / VMF grids''';

        %CONCATENATE STRINGS
        string = sprintf('%s\n%s\n%s',string1, string2,string3); 
        
        %SET TOOL TIP
        set(hObject,'TooltipString',string)%change Tooltip any time user select Model type 
        
        %CHECK IF VMF GRID FILEs ARE AVAILABLE IN goGPS
        %Call for the 'SearchVMFgrids.m' fxn
        [VMF_grid_found,VMFgrids] = SearchVMFgrids(handles) ;
                        
        %*********SAVE STRUCT GRID FILES
        setappdata(0,'VMFgrids',VMFgrids)
         
        %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
        setappdata(0,'VMF_grid_found',VMF_grid_found)
        
        if VMF_grid_found == 1 %IF VMF grid FILES FOUND,LET USER CHOOSE GRID TYPE & RESOLUTION
                               %I.E.['VMF1(2°x2.5°)','VMF3(1°x1°)','VMF3(5°x5°)']
                               
           [VMF_type,grid_res] = VMFgrid_type(); 
           
           %*********CHECK IF VMF_type,grid_res ARE BOTH EMPTY
           %NOTE:
           %     VMF_type,grid_res ARE ASSIGNED EMPTY WHEN USER CANCEL SELECTION
           %     OR CLOSE DIALOGUE BOX FOR THE SELECTION OF VMF GRID TYPE &
           %     RESOLUTION.IF IT HAPPENS SO, THE SELECTED TROPO DELAY MODEL,
           %     'VMF gridded ZTD' WOULD BE REPLACED WITH SAASTEMOIN MODEL
           %     OR THE POPUP MENU WILL BE SET TO SAASTAMOINEN AS DEFAULT
        
           if all([isempty(VMF_type),isempty(grid_res)])%IF USER CANCEL SELECTION/CLOSE DIALOGUE BOX
            
              %SET WET TROPO MODEL POPUP MENU VALUE TO 1(I.E.SAASTAMOINEN)
              set(handles.pop_ITropo_wet,'Value',1)
              
              if any([any([any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv_H == 15]),...
                           any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv == 20])]),...
                      all([any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv_H == 15]),...
                           any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv == 20])])])
                
                  %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
                  set(handles.cb_extract_VMF_ZTDs,'enable','on','value',0)
                  
              end
              
              %******SET MET PARAMETER SOURCE UICONTROLS VISIBLE/ENABLE ON
              set(handles.pop_source_metpara,'enable','on')
              set(handles.text_source_metpara,'enable','on')
    
              if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
                              any(metVAL==[7,8])])
        
                 if strcmpi(get(handles.pop_met_manual,'visible'),'Off') 
                    set(handles.pop_met_manual,'visible','On')
                    set(handles.text_metmanual_source,'visible','On')
                 end     
     
              elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
                                  any(metVAL==[4,5])])
          
                     if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'off'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'off')])
                        set(handles.rb_grid_resolution_1,'Visible','on')
                        set(handles.rb_grid_resolution_5,'Visible','on')
                        set(handles.text_grid_resolution,'Visible','on')  
                     end      
          
              end  %//if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
                   %                  any(metVAL==[7,8])])
                      
           else
               
               if any([all([~any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv_H == 15]),...
                            any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv == 20])]),...
                       all([any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv_H == 15]),...
                            ~any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv == 20])])])
                
                  %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
                  set(handles.cb_extract_VMF_ZTDs,'enable','on','value',0)
                  
               elseif all([any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv_H == 15]),...
                           any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv == 20])])
                       
                      %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
                      set(handles.cb_extract_VMF_ZTDs,'enable','off','value',0)
                  
              end
                   
                   
           end  %//if all([isempty(VMF_type),isempty(grid_res)])
        
           %SAVE GRID FILE VERSION & TYPE
           setappdata(0,'VMFgrid_type',VMF_type)
           setappdata(0,'VMFgrid_res',grid_res)
           
        end  %//if VMF_grid_found == 1
                         
elseif  any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),Modelv == 21])
        
        %SET TOOLTIP STRING
        set(hObject,'TooltipString','Solution using Global Tropospheric model(GTrop)')%change Tooltip any time user select Model type 
       
        %GET & SAVE GTrop GRID VALUES
        if ~isempty(getappdata(0,'gridV_GTrop'))
            
           setappdata(0,'gridV_GTrop',getappdata(0,'gridV_GTrop')) 
           
        else
            %SEARCH FOR GTrop model GRID FILE(GTropCoefficient.mat) & DIRECTORY & LOAD VALUES               
            %Call the "SearchGTropgrid.m" fxn
            [~,gridV_GTrop] = SearchGTropgrid();%****LOADED GTrop MODEL COEFFICIENTS
       
            %SAVE GRID VALUES
            setappdata(0,'gridV_GTrop',gridV_GTrop) 
            
        end
                    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%***********GLOBAL PRESSURE  AND TEMPERATURE(GPT) GRID MODELS
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+    
elseif Modelv == 15 || strncmpi(wetModel,'GPT2 (5° x 5°)',14)
       set(hObject,'TooltipString',strjoin({'Solution using','Global Pressure and Temperature v2 model', 'based on a (5° x 5°) external grid file'}))%change Tooltip any time user select Model type 
       
       %GET GPT MODELs GRID FILE & RESOLUTION
       grid_res_w=5; %SET GRID RESOLUTION TO 5
       
      %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       if ~isempty(grid_w1) & grid_res_w1 == 5
           
           if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(wetModel,wetModel1,14)
              grid_w = grid_w1;
              
           else 
                %****TRY FROM OTHER MODELS
                %1.HYDROSTATIC MODELs
                if any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(wetModel,dryModel1,14)
                   if ~isempty(grid_h1) & grid_res_h1 == 5
                     grid_w = grid_h1;
                   else
                       grid_w=readGPTgrid('gpt2_5.mat','GPT2',5);%Wet model GRIDFILE
                   end
                %2.COMBINE MODELs
                elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(tropoModel1,'GPT2 (5° x 5°)',14)
                       if ~isempty(grid_c) & grid_res_c == 5
                          grid_w = grid_c;
                       else
                           grid_w=readGPTgrid('gpt2_5.mat','GPT2',5);%Wet model GRIDFILE
                       end
                       
                else %IF ALL FAIL,READ GRID
                     grid_w=readGPTgrid('gpt2_5.mat','GPT2',5);%WET GRIDFILE
                     
                end
           end
           
       else %IF WET GRIDFILE IS EMPTY([])
           if any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(wetModel,dryModel1,14)
              if ~isempty(grid_h1) & grid_res_h1 == 5
                 grid_w = grid_h1;
              else
                  grid_w=readGPTgrid('gpt2_5.mat','GPT2',5);%Wet model GRIDFILE
              end 
           
           elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(tropoModel1,'GPT2 (5° x 5°)',14)
              if ~isempty(grid_c) & grid_res_c == 5
                 grid_w = grid_c;
              else
                  grid_w=readGPTgrid('gpt2_5.mat','GPT2',5);%Wet model GRIDFILE
              end
              
           else 
               grid_w=readGPTgrid('gpt2_5.mat','GPT2',5);%WET GRIDFILE
               
           end  
       end
       
       %IF ALL OPTIONS FAIL
      if ~exist('grid_w','var')
         grid_w=readGPTgrid('gpt2_5.mat','GPT2',5);%Wet model GRIDFILE
      end   
      
       if length(grid_w) > 8 
          grid_w=readGPTgrid('gpt2_5.mat','GPT2',5);%Wet model GRIDFILE
       end 
       
%SAVE WET MODEL GRID FILE(grid_w) & RESOLUTION(grid_res_w)
setappdata(0,'grid_w',grid_w)
setappdata(0,'grid_res_w',grid_res_w)         
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+  

elseif Modelv == 16 || strncmpi(wetModel,'GPT2w (1° x 1°)',14)
       set(hObject,'TooltipString',strjoin({'Solution using','Global Pressure and Temperature to wet model', 'based on a (1° x 1°) external grid file'}))%change Tooltip any time user select Model type 
    
       %GET GPT MODELs GRID FILE & RESOLUTION
       grid_res_w=1; %SET GRID RESOLUTION TO 1
       
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       if ~isempty(grid_w1) & grid_res_w1 == 1
           
           if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(wetModel,wetModel1,14)
              grid_w = grid_w1;
              
           else 
                %****TRY FROM OTHER MODELS
                %1.HYDROSTATIC MODELs
                if any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(wetModel,dryModel1,14)
                   if ~isempty(grid_h1) & grid_res_h1 == 1
                      grid_w = grid_h1;
                   else
                       grid_w=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Wet model GRIDFILE
                   end
                %2.COMBINE MODELs
                elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(tropoModel1,'GPT2w (1° x 1°)',14)
                       if ~isempty(grid_c) & grid_res_c == 1
                          grid_w = grid_c;
                       else
                           grid_w=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Wet model GRIDFILE
                       end
                       
                else %IF ALL FAIL,READ GRID
                     grid_w=readGPTgrid('gpt2_1w.mat','GPT2w',1);%WET GRIDFILE
                     
                end
           end
           
       else %IF WET GRIDFILE IS EMPTY([])
           if any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(wetModel,dryModel1,14)
              if ~isempty(grid_h1) & grid_res_h1 == 1
                 grid_w = grid_h1;
              else
                  grid_w=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Wet model GRIDFILE
              end 
           
           elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(tropoModel1,'GPT2w (1° x 1°)',14)
              if ~isempty(grid_c) & grid_res_c == 1
                 grid_w = grid_c;
              else
                  grid_w=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Wet model GRIDFILE
              end 
              
           else 
               grid_w=readGPTgrid('gpt2_1w.mat','GPT2w',1);%WET GRIDFILE
               
           end  
       end
       
      %IF ALL OPTIONS FAIL
      if ~exist('grid_w','var')
         grid_w = readGPTgrid('gpt2_1w.mat','GPT2w',1);%Wet model GRIDFILE
      end   
       
      if any([length(grid_w) > 10,length(grid_w)< 10])
         grid_w=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Wet model GRIDFILE
      end      
       
%SAVE WET MODEL GRID FILE(grid_w) & RESOLUTION(grid_res_w)
setappdata(0,'grid_w',grid_w)
setappdata(0,'grid_res_w',grid_res_w)         
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+  

elseif Modelv == 17 || strncmpi(wetModel,'GTP2w (5° x 5° )',14)
       set(hObject,'TooltipString',strjoin({'Solution using','Global Pressure and Temperature to wet model', 'based on a (5° x 5°) external grid file'}))%change Tooltip any time user select Model type 
    
       %GET GPT MODELs GRID FILE & RESOLUTION
       grid_res_w=5; %SET GRID RESOLUTION TO 5
       
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       if ~isempty(grid_w1) & grid_res_w1 == 5
           
           if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(wetModel,wetModel1,14)
              grid_w = grid_w1;
              
           else 
                %****TRY FROM OTHER MODELS
                %1.HYDROSTATIC MODELs
                if any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(wetModel,dryModel1,14)
                   if ~isempty(grid_h1) & grid_res_h1 == 5
                      grid_w = grid_h1;
                   else
                       grid_w=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Wet model GRIDFILE
                   end
                %2.COMBINE MODELs
                elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(tropoModel1,'GPT2w (5° x 5°)',14)
                       if ~isempty(grid_c) & grid_res_c == 5
                          grid_w = grid_c;
                       else
                           grid_w=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Wet model GRIDFILE
                       end
                       
                else %IF ALL FAIL,READ GRID
                     grid_w=readGPTgrid('gpt2_5w.mat','GPT2w',5);%WET GRIDFILE
                     
                end
           end
           
       else %IF WET GRIDFILE IS EMPTY([])
           if any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(wetModel,dryModel1,14)
              if ~isempty(grid_h1) & grid_res_h1 == 5
                 grid_w = grid_h1;
              else
                  grid_w=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Wet model GRIDFILE
              end 
           
           elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(tropoModel1,'GPT2w (5° x 5°)',14)
              if ~isempty(grid_c) & grid_res_c == 5
                 grid_w = grid_c;
              else
                  grid_w=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Wet model GRIDFILE
              end
              
           else 
               grid_w=readGPTgrid('gpt2_5w.mat','GPT2w',5);%WET GRIDFILE
               
           end  
       end
       
      %IF ALL OPTIONS FAIL
      if ~exist('grid_w','var')
         grid_w = readGPTgrid('gpt2_5w.mat','GPT2w',5);%Wet model GRIDFILE
      end   
       
      if any([length(grid_w) > 10,length(grid_w)< 10])
         grid_w=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Wet model GRIDFILE
      end      
 
%SAVE WET MODEL GRID FILE(grid_w) & RESOLUTION(grid_res_w)
setappdata(0,'grid_w',grid_w)
setappdata(0,'grid_res_w',grid_res_w)         
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+  

elseif Modelv == 18 || strncmpi(wetModel,'GPT3 (1° x 1°)',14)
       set(hObject,'TooltipString',strjoin({'Solution using','Global Pressure and Temperature v3 model', 'based on a (1° x 1°) external grid file'}))%change Tooltip any time user select Model type  

       %GET GPT MODELs GRID FILE & RESOLUTION
       grid_res_w=1; %SET GRID RESOLUTION TO 1
       
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       if ~isempty(grid_w1) & grid_res_w1 == 1
           
           if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(wetModel,wetModel1,14)
              grid_w = grid_w1;
              
           else 
                %****TRY FROM OTHER MODELS
                %1.HYDROSTATIC MODELs
                if any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(wetModel,dryModel1,14)
                   if ~isempty(grid_h1) & grid_res_h1 == 1
                     grid_w = grid_h1;
                   else
                       grid_w=readGPTgrid('gpt3_1.mat','GPT3',1);%Wet model GRIDFILE
                   end
                %2.COMBINE MODELs
                elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(tropoModel1,'GPT3 (1° x 1°)',14)
                       if ~isempty(grid_c) & grid_res_c == 1
                          grid_w = grid_c;
                       else
                           grid_w=readGPTgrid('gpt3_1.mat','GPT3',1);%Wet model GRIDFILE
                       end
                       
                else %IF ALL FAIL,READ GRID
                     grid_w=readGPTgrid('gpt3_1.mat','GPT3',1);%WET GRIDFILE
                     
                end
           end
           
       else %IF WET GRIDFILE IS EMPTY([])
           if any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(wetModel,dryModel1,14)
              if ~isempty(grid_h1) & grid_res_h1 == 1
                 grid_w = grid_h1;
              else
                  grid_w=readGPTgrid('gpt3_1.mat','GPT3',1);%Wet model GRIDFILE
              end 
           
           elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(tropoModel1,'GPT3 (1° x 1°)',14)
              if ~isempty(grid_c) & grid_res_c == 1
                 grid_w = grid_c;
              else
                  grid_w=readGPTgrid('gpt3_1.mat','GPT3',1);%Wet model GRIDFILE
              end
              
           else 
               grid_w=readGPTgrid('gpt3_1.mat','GPT3',1);%WET GRIDFILE
               
           end  
       end
       
       %IF ALL OPTIONS FAIL
      if ~exist('grid_w','var')
         grid_w=readGPTgrid('gpt3_1.mat','GPT3',1);%Wet model GRIDFILE
      end   
   
      if length(grid_w)<14
         grid_w=readGPTgrid('gpt3_1.mat','GPT3',1);%Wet model GRIDFILE
      end   
 
%SAVE WET MODEL GRID FILE(grid_w) & RESOLUTION(grid_res_w)
setappdata(0,'grid_w',grid_w)
setappdata(0,'grid_res_w',grid_res_w)      
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+  

elseif Modelv == 19 || strncmpi(wetModel,'GPT3 (5° x 5°)',14)
       set(hObject,'TooltipString',strjoin({'Solution using','Global Pressure and Temperature v3 model', 'based on a (5° x 5°) external grid file'}))%change Tooltip any time user select Model type
       
       %GET GPT MODELs GRID FILE & RESOLUTION
       grid_res_w=5; %SET GRID RESOLUTION TO 5
       
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       if ~isempty(grid_w1) & grid_res_w1 == 5
           
           if any([~isempty(wetModel1) ~isnan(wetModel1)]) & strncmpi(wetModel,wetModel1,14)
              grid_w = grid_w1;
              
           else 
                %****TRY FROM OTHER MODELS
                %1.HYDROSTATIC MODELs
                if any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(wetModel,dryModel1,14)
                   if ~isempty(grid_h1) & grid_res_h1 == 5
                     grid_w = grid_h1;
                   end
                %2.COMBINE MODELs
                elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(tropoModel1,'GPT3 (5° x 5°)',14)
                       if ~isempty(grid_c) & grid_res_c == 5
                          grid_w = grid_c;
                       end
                       
                else %IF ALL FAIL,READ GRID
                     grid_w=readGPTgrid('gpt3_5.mat','GPT3',5);%WET GRIDFILE
                     
                end
           end
           
       else %IF WET GRIDFILE IS EMPTY([])
           if any([~isempty(dryModel1) ~isnan(dryModel1)]) & strncmpi(wetModel,dryModel1,14)
              if ~isempty(grid_h1) & grid_res_h1 == 5
                 grid_w = grid_h1;
              end 
           
           elseif any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & strncmpi(tropoModel1,'GPT3 (5° x 5°)',14)
              if ~isempty(grid_c) & grid_res_c == 5
                 grid_w = grid_c;
              end 
           else 
               grid_w=readGPTgrid('gpt3_5.mat','GPT3',5);%HYDROSTATIC GRIDFILE
           end  
       end
       
       if length(grid_w)<14
          grid_w=readGPTgrid('gpt3_5.mat','GPT3',5);%Combine model GRIDFILE
       end    
 
%SAVE WET MODEL GRID FILE(grid_w) & RESOLUTION(grid_res_w)
setappdata(0,'grid_w',grid_w)
setappdata(0,'grid_res_w',grid_res_w)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

end

%SAVE SELECTED WET TROPO MODEL
setappdata(0,'wetModel',wetModel)


%CHECK IF ANY OF THE GPT MODEL(GPT2,GPT2w & GPT3) IS SELECTED,AND ALLOW
%USER TO INDICATE WHETHER TO APPLY TIME VARIATION OF NOT
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%NOTE:
%    GLOBAL PRESSURE & TEMPERATURE(GPT) MODELS SUCH AS GPT2,GPT2w & GPT3 
%    COMPUTES METEOROLOGICAL PARAMETERS EITHER IN STATIC MODE OR IN TIME
%    VARYING MODE:
%    case 1: no time variation but static quantities
%    case 0: with time variation (annual and semiannual terms)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
        strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]) & ...   
    any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
         strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])


if ~isempty(grid_w)
    
   TIMEvariation()
   uiwait(gcf)
   setappdata(0,'Timevar_wet',getappdata(0,'Timevar'))
   
end

elseif ~any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
             strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14)]) & ...   
       any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
             strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14)])

     if ~isempty(grid_w) 
         
        TIMEvariation()
        uiwait(gcf)
        setappdata(0,'Timevar_wet',getappdata(0,'Timevar'))  
        
     end
end

%-------------------------------------------------------------------  
%NOTE"
%SAASTAMOINEN,DAVIS ET AL,Askne & Nordius,UNB3m,EGNOS & MOPS MODELS 
%REQUIRE +VE ORTHOMETRIC HEIGHT-------------------------------------
%********GET goGPS DEFAULT GEOID MODEL (EGM2008) 

%CHECK SELECTE  MODEL
if any([strncmpi(wetModel,'Saastamoinen',12),strncmpi(wetModel,'Askne & Nordius',15),...
        strncmpi(wetModel,'UNB3m',5),strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4)])          
          
    %LOAD GEOID MODEL[%USING goGPS DEFAULT GEOID MODEL (EGM2008)] 
    if ~isempty(getappdata(0,'geoid'))
               
       %SAVE FOR FUTURE USE
       setappdata(0,'geoid',getappdata(0,'geoid')) 
              
    else           
         try %goGPS v0.5.2 beta1
             gs = Go_State.getInstance;
          geoid = gs.getRefGeoid(); 
      
         catch  %goGPS v0.4.3 
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

         %SAVE FOR FUTURE USE
         setappdata(0,'geoid',geoid)
             
    end %/if ~isempty(getappdata(0,'geoid'))    
       
end  %//if any([strncmpi(wetModel,'Saastamoinen',12),strncmpi(wetModel,'Askne & Nordius',15),...
     %          strncmpi(wetModel,'UNB3m',5),strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4)])
          
%===================================END OF TROPOmodels_wet__popm_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 

%                              -----------------
%(2.6) ************************MAPPING FUNCTIONS
%                              -----------------

%                    ------------------------------
%(2.6.1)=============MAPPING FUNCTION RADIO BUTTONS========================
%                    ------------------------------
%            --------------------------
%(A)*********USE MODEL MAPPING FUNCTION
%            --------------------------
% --- EXECUTES ON BUTTON PRESS IN rb_model_MF.
function TROPOmodels_mappingfunction_rb_Callback(hObject,handles)
%--------------------------------------------------------------------------
set(hObject,'Value',1,'TooltipString','Solution with Delay Models Mapping Function')
set(handles.rb_different_MF,'Value',0,'TooltipString','Select to apply different Mapping Function Models');%Model mapping function(MF) Radio Button
%*********TURN OFF SOME TOOLs
%1.MAPPING FUNCTION
set(handles.pop_IHydrostic_MF,'enable','off')%Hydrostatic MF popup menu
set(handles.pop_IWet_MF,'enable','off')%Wet MF popup menu
set(handles.text__Hydrostatic_MF,'enable','off');%Hydrostatic MF static text
set(handles.text__Wet_MF,'enable','off');%Hydrostatic MF static text

%==========================END OF TROPOmodels_mappingfunction_rb_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 

%                ------------------------
%(B) ************USE OTHER MAPPING MODELS
%                -------------------------
% --- EXECUTES ON BUTTON PRESS IN rb_different_MF.
function OTHERmodels_mappingfunction_rb_Callback(hObject,handles)
%--------------------------------------------------------------------------
set(hObject,'Value',1,'TooltipString','Solution with different Mapping Function Models')
set(handles.rb_model_MF,'Value',0,'TooltipString','Select to apply Delay Model Mapping Function');%Model mapping function(MF) Radio Button
%1.MAPPING FUNCTION
set(handles.pop_IHydrostic_MF,'enable','on','TooltipString','')%Hydrostatic MF popup menu
set(handles.pop_IWet_MF,'enable','on','TooltipString','')%Wet MF popup menu
set(handles.text__Hydrostatic_MF,'enable','on');%Hydrostatic MF static text
set(handles.text__Wet_MF,'enable','on');%Hydrostatic MF static text

%==========================END OF OTHERmodels_mappingfunction_rb_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%                     ----------------------------
%(2.6.2)===============MAPPING FUNCTION POPUP MENUs========================
%                      -----------------------------

%                ----------------------------
%(A)*************HYDROSTATIC MAPPING FUNCTION
%                ----------------------------
% --- EXECUTES ON SELECTION CHANGE in pop_IHydrostic_MF.
function HYDROSTATIC_mappingfunction_popm_Callback(hObject,handles)
%--------------------------------------------------------------------------
global MFh_model

%*****GET VARIOUS UI CONTROL

%1.Hydrostatic MF popup menu contents
getcontsHMF = cellstr(get(hObject,'String'));%get Hydrostatic MF popup menu contents as cell array
MFh_model   = getcontsHMF{get(hObject,'Value')};%get selected item from pop_IHydrostic_MF
MFh_val     = get(hObject,'Value'); %get selected item(value) from pop_IHydrostic_MF

%2.MET PARAMETER SOURCE POPUP MENU CONTENTs
getMetconts = cellstr(get(handles.pop_source_metpara,'String'));%get Met Para popup menu contents as cell array
metPARA_str = getMetconts{get(handles.pop_source_metpara,'Value')};%get selected item from Met Para popup menu
metPARA_val = get(handles.pop_source_metpara,'Value');%get selected item from Met Para popup menu

%GET GRID RESOLUTIONs FORM RADIO BUTTOMS
grid_res_1=get(handles.rb_grid_resolution_1,'Value');%1° grid resolution 
grid_res_5=get(handles.rb_grid_resolution_5,'Value');%5° grid resolution 
gridModel    = getappdata(0,'gridModel');%grid model

%GET PREVIOUS FILEs
grid_MFh1    = getappdata(0,'grid_MFh');%previous gridfile
gridRES_MFh1 = getappdata(0,'gridRES_MFh');%previous grid resolution
 

%SET TOOLTIP STRING
set(hObject,'TooltipString',strjoin({'Solution using',MFh_model,'model'}))%change Tooltip any time user select Model type

if strncmpi(MFh_model,'GMF',3)
   set(hObject,'TooltipString',strjoin({'Solution using','Global Mapping Function(GMF)','model'}))%change Tooltip any time user select Model type
   
elseif strncmpi(MFh_model,'NMF',3)
       set(hObject,'TooltipString',strjoin({'Solution using','Neill Mapping Function(NMF)','model'}))%change Tooltip any time user select Model type 

elseif strncmpi(MFh_model,'VMF1(1° x 1°)',12)
       set(hObject,'TooltipString',strjoin({'Solution using','Vienna Mapping Function1( VMF1(1° x 1°) )','model'}))%change Tooltip any time user select Model type
 
elseif strncmpi(MFh_model,'VMF1(5° x 5°)',12)
       set(hObject,'TooltipString',strjoin({'Solution using','Vienna Mapping Function1( VMF1(5° x 5°) )','model'}))%change Tooltip any time user select Model type
    
elseif strncmpi(MFh_model,'VMF3(1° x 1°)',12)
       set(hObject,'TooltipString',strjoin({'Solution using','Vienna Mapping Function1( VMF3(1° x 1°) )','model'}))%change Tooltip any time user select Model type
           
elseif strncmpi(MFh_model,'VMF3(5° x 5°)',12)
       set(hObject,'TooltipString',strjoin({'Solution using','Vienna Mapping Function1( VMF3(5° x 5°) )','model'}))%change Tooltip any time user select Model type
       
end

%**********SEARCH FOR THE EXISTENCE OF VMF GRID FILES
if any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),...
        strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12)])
    
   
    %CHECK IF VMF GRID FILEs ARE AVAILABLE IN goGPS
    %Call for the 'SearchVMFgrids.m' fxn
    [VMF_grid_found,VMFgrids] = SearchVMFgrids(handles);
                       
    %*********SAVE STRUCT GRID FILES
    setappdata(0,'VMFgrids',VMFgrids)
         
    %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
    setappdata(0,'VMF_grid_found',VMF_grid_found)
    
    if VMF_grid_found == 1 %IF VMF grid FILES FOUND
                           
       if any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12)])
          
          VMF_type_mf = 'VMF1';
          grid_res_mf = 2; 
          
       elseif strncmpi(MFh_model,'VMF3(1° x 1°)',12)
              
              VMF_type_mf = 'VMF3';
              grid_res_mf = 1;
              
       elseif strncmpi(MFh_model,'VMF3(5° x 5°)',12)
                      
              VMF_type_mf = 'VMF3';
              grid_res_mf = 5;
              
       end %//if any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12)])
       
       
       %SAVE GRID FILE VERSION & TYPE
       setappdata(0,'VMFgrid_type_mf',VMF_type_mf)
       setappdata(0,'VMFgrid_res_mf',grid_res_mf)
       
    end %//if VMF_grid_found == 1
    
                   
end %//if any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),...
    %          strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12)])

    %GET GPT MODELS GRID VALUES IF VMF GRID FILES ARE NOT FOUND   
if strcmpi(get(hObject,'Enable'),'on') & strncmpi(MFh_model,'VMF1(1° x 1°)',12)
    
   if ~isempty(grid_MFh1) & ~isempty(gridRES_MFh1) 
       
       if gridRES_MFh1 == 1
           
          grid_MFh    = grid_MFh1;
          gridRES_MFh = gridRES_MFh1;
          
       else
           %*****READ GRID FILE
           grid_MFh    = readGPTgrid('gpt2_1w.mat','GPT2w',1);%Hydrostatic GRIDFILE
           gridRES_MFh = 1;%grid resolution
       end
       
   else
           
       if any(metPARA_val==[4,5])
      
          if grid_res_1 == 1  & strcmpi(gridModel,'GPT2w')
             grid_MFh    = getappdata(0,'gridFile'); %gridfile
             gridRES_MFh = getappdata(0,'gridRES');%grid resolution
         
          else 
              %*****READ GRID FILE
              grid_MFh    = readGPTgrid('gpt2_1w.mat','GPT2w',1);%Hydrostatic GRIDFILE
              gridRES_MFh = 1;%grid resolution
          end 
      
       else 
           %*****READ GRID FILE
           grid_MFh    = readGPTgrid('gpt2_1w.mat','GPT2w',1);%Hydrostatic GRIDFILE
           gridRES_MFh = 1;%grid resolution
       
       end %//if any(metPARA_val==[4,5])
       
   end %//if ~isempty(grid_MFh1) & ~isempty(gridRES_MFh1)
   
   if ((length(grid_MFh) > 10) | (length(grid_MFh)< 10))
      grid_MFh=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Combine model GRIDFILE
   end  
   
elseif strcmpi(get(hObject,'Enable'),'on') & strncmpi(MFh_model,'VMF1(5° x 5°)',12)
       
      if ~isempty(grid_MFh1) & ~isempty(gridRES_MFh1) 
       
         if gridRES_MFh1 == 5
           
             grid_MFh    = grid_MFh1;
             gridRES_MFh = gridRES_MFh1;
          
         else 
             %*****READ GRID FILE
             grid_MFh    = readGPTgrid('gpt2_5w.mat','GPT2w',5);%Hydrostatic GRIDFILE
             gridRES_MFh = 5;%grid resolution
         end 
       
      else 
           if any(metPARA_val==[4,5])
      
              if grid_res_5 == 1  & strcmpi(gridModel,'GPT2w')
                 grid_MFh    = getappdata(0,'gridFile'); %gridfile
                 gridRES_MFh = getappdata(0,'gridRES');%grid resolution
         
              else  
                  %*****READ GRID FILE
                  grid_MFh    = readGPTgrid('gpt2_5w.mat','GPT2w',5);%Hydrostatic GRIDFILE
                  gridRES_MFh = 5;%grid resolution
              end  
      
           else   
                %*****READ GRID FILE
                grid_MFh    = readGPTgrid('gpt2_5w.mat','GPT2w',5);%Hydrostatic GRIDFILE
                gridRES_MFh = 5;%grid resolution
       
           end  %//if any(metPARA_val==[4,5])
           
      end %//if ~isempty(grid_MFh1) & ~isempty(gridRES_MFh1)
       
       if ((length(grid_MFh) > 10) | (length(grid_MFh)< 10))
          grid_MFh=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Combine model GRIDFILE
       end  
   
elseif strcmpi(get(hObject,'Enable'),'on') & strncmpi(MFh_model,'VMF3(1° x 1°)',12)
   
      if ~isempty(grid_MFh1) & ~isempty(gridRES_MFh1) 
       
         if gridRES_MFh1 == 1
           
             grid_MFh    = grid_MFh1;
             gridRES_MFh = gridRES_MFh1;
          
         else 
             %*****READ GRID FILE
             grid_MFh    = readGPTgrid('gpt3_1.mat','GPT3',1);%Hydrostatic GRIDFILE
             gridRES_MFh = 1;%grid resolution
         end 
       
      else  
          if any(metPARA_val==[4,5])
      
             if grid_res_1 == 1  & strcmpi(gridModel,'GPT3')
                grid_MFh    = getappdata(0,'gridFile'); %gridfile
                gridRES_MFh = getappdata(0,'gridRES');%grid resolution
             else   
                 %*****READ GRID FILE
                 grid_MFh    = readGPTgrid('gpt3_1.mat','GPT3',1);%Hydrostatic GRIDFILE
                 gridRES_MFh = 1;%grid resolution
             end   
      
          else  
              %*****READ GRID FILE
              grid_MFh    = readGPTgrid('gpt3_1.mat','GPT3',1);%Hydrostatic GRIDFILE
              gridRES_MFh = 1;%grid resolution
       
          end %//if any(metPARA_val==[4,5]) 
          
      end %//if ~isempty(grid_MFh1) & ~isempty(gridRES_MFh1) 
      
       if length(grid_MFh)<14
          grid_MFh=readGPTgrid('gpt3_1.mat','GPT3',1);%Combine model GRIDFILE
       end    
   
elseif strcmpi(get(hObject,'Enable'),'on') & strncmpi(MFh_model,'VMF3(5° x 5°)',12)
       
      if ~isempty(grid_MFh1) & ~isempty(gridRES_MFh1) 
       
         if gridRES_MFh1 == 5
           
             grid_MFh    = grid_MFh1;
             gridRES_MFh = gridRES_MFh1;
          
         else  
             %*****READ GRID FILE
             grid_MFh    = readGPTgrid('gpt3_5.mat','GPT3',5);%Hydrostatic GRIDFILE
             gridRES_MFh = 5;%grid resolution
         end 
       
      else 
          
           if any(metPARA_val==[4,5])
      
             if grid_res_5 == 1  & strcmpi(gridModel,'GPT3')
                grid_MFh    = getappdata(0,'gridFile'); %gridfile
                gridRES_MFh = getappdata(0,'gridRES');%grid resolution
         
             else  
                 %*****READ GRID FILE
                 grid_MFh    = readGPTgrid('gpt3_5.mat','GPT3',5);%Hydrostatic GRIDFILE
                 gridRES_MFh = 5;%grid resolution
             end  
      
           else    
               %*****READ GRID FILE
               grid_MFh    = readGPTgrid('gpt3_5.mat','GPT3',5);%Hydrostatic GRIDFILE
               gridRES_MFh = 5;%grid resolution
       
           end %//if any(metPARA_val==[4,5])  
           
      end %//if ~isempty(grid_MFh1) & ~isempty(gridRES_MFh1)
    
      if length(grid_MFh)<14
         grid_MFh=readGPTgrid('gpt3_5.mat','GPT3',5);%Combine model GRIDFILE
      end     
      
      
end %//if strcmpi(get(hObject,'Enable'),'on') & strncmpi(MFh_model,'VMF1(1° x 1°)',12)

%SAVE GRID FILE & RESOLUTION 
if any([exist('grid_MFh','var') exist('gridRES_MFh','var')])
   setappdata(0,'grid_MFh',grid_MFh)
   setappdata(0,'gridRES_MFh',gridRES_MFh) 
end

%SAVE SELECTED HYDROSTATIC MAPPING FUNCTION MODEL
setappdata(0,'MFh_model',MFh_model)

%==========================END OF HYDROSTATIC_mappingfunction_popm_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%                      --------------------
%(B)*******************WET MAPPING FUNCTION
%                      ---------------------
% --- EXECUTES ON SELECTION CHANGE in pop_pop_IWet_MF.
function WET_mappingfunction_popm_Callback(hObject,handles)
%--------------------------------------------------------------------------
global MFw_model

%*****GET VARIOUS UI CONTROL

%1.Hydrostatic MF popup menu contents
getcontsWMF = cellstr(get(hObject,'String'));%get Wet MF popup menu contents as cell array
MFw_model=getcontsWMF{get(hObject,'Value')};%get selected item from pop_IWet_MF
MFw_val     = get(hObject,'Value'); %get selected item(value) from  pop_IWet_MF

%2.MET PARAMETER SOURCE POPUP MENU CONTENTs
getMetconts = cellstr(get(handles.pop_source_metpara,'String'));%get Met Para popup menu contents as cell array
metPARA_str = getMetconts{get(handles.pop_source_metpara,'Value')};%get selected item from Met Para popup menu
metPARA_val = get(handles.pop_source_metpara,'Value');%get selected item from Met Para popup menu

%GET GRID RESOLUTIONs FORM RADIO BUTTOMS
grid_res_1=get(handles.rb_grid_resolution_1,'Value');%1° grid resolution 
grid_res_5=get(handles.rb_grid_resolution_5,'Value');%5° grid resolution 
gridModel    = getappdata(0,'gridModel');%grid model

%GET PREVIOUS FILEs
grid_MFw1    = getappdata(0,'grid_MFw');%previous gridfile
gridRES_MFw1 = getappdata(0,'gridRES_MFw');%previous grid resolution

%SET TOOLTIP STRING
set(hObject,'TooltipString',strjoin({'Solution using',MFw_model,'model'}))%change Tooltip any time user select Model type

if strncmpi(MFw_model,'GMF',3)
   set(hObject,'TooltipString',strjoin({'Solution using','Global Mapping Function(GMF)','model'}))%change Tooltip any time user select Model type
   
elseif strncmpi(MFw_model,'NMF',3)
       set(hObject,'TooltipString',strjoin({'Solution using','Neill Mapping Function(NMF)','model'}))%change Tooltip any time user select Model type

elseif strncmpi(MFw_model,'VMF1(1° x 1°)',12)
       set(hObject,'TooltipString',strjoin({'Solution using','Vienna Mapping Function1( VMF1(1° x 1°) )','model'}))%change Tooltip any time user select Model type
 
elseif strncmpi(MFw_model,'VMF1(5° x 5°)',12)
       set(hObject,'TooltipString',strjoin({'Solution using','Vienna Mapping Function1( VMF1(5° x 5°) )','model'}))%change Tooltip any time user select Model type
    
elseif strncmpi(MFw_model,'VMF3(1° x 1°)',12)
       set(hObject,'TooltipString',strjoin({'Solution using','Vienna Mapping Function1( VMF3(1° x 1°) )','model'}))%change Tooltip any time user select Model type
           
elseif strncmpi(MFw_model,'VMF3(5° x 5°)',12)
       set(hObject,'TooltipString',strjoin({'Solution using','Vienna Mapping Function1( VMF3(5° x 5°) )','model'}))%change Tooltip any time user select Model type      
end

%*************SEARCH FOR THE EXISTENCE OF VMF GRID FILES
if any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),...
        strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12)])
   
  
    %CHECK IF VMF GRID FILEs ARE AVAILABLE IN goGPS
    %Call for the 'SearchVMFgrids.m' fxn
    [VMF_grid_found,VMFgrids] = SearchVMFgrids(handles);
                        
    %*********SAVE STRUCT GRID FILES
    setappdata(0,'VMFgrids',VMFgrids)
         
    %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
    setappdata(0,'VMF_grid_found',VMF_grid_found)
    
    if VMF_grid_found == 1 %IF VMF grid FILES FOUND
                           
       if any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12)])
          
          VMF_type_mf = 'VMF1';
          grid_res_mf = 2; 
          
       elseif strncmpi(MFw_model,'VMF3(1° x 1°)',12)
              
              VMF_type_mf = 'VMF3';
              grid_res_mf = 1;
              
       elseif strncmpi(MFw_model,'VMF3(5° x 5°)',12)
                      
              VMF_type_mf = 'VMF3';
              grid_res_mf = 5;
              
       end %//if any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12)])
       
       
       %SAVE GRID FILE VERSION & TYPE
       setappdata(0,'VMFgrid_type_mf',VMF_type_mf)
       setappdata(0,'VMFgrid_res_mf',grid_res_mf)
       
    end %//if VMF_grid_found == 1
                  
end %//if any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),...
    %          strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12)])
 
  %GET GPT MODELS GRID VALUES IF VMF GRID FILES ARE NOT FOUND   
if strcmpi(get(hObject,'Enable'),'on') & strncmpi(MFw_model,'VMF1(1° x 1°)',12)
    
   if ~isempty(grid_MFw1) & ~isempty(gridRES_MFw1) 
       
       if gridRES_MFw1 == 1
           
          grid_MFw    = grid_MFw1;
          gridRES_MFw = gridRES_MFw1;
          
       else
           %*****READ GRID FILE
           grid_MFw    = readGPTgrid('gpt2_1w.mat','GPT2w',1);%Hydrostatic GRIDFILE
           gridRES_MFw = 1;%grid resolution
       end
       
   else
           
       if any(metPARA_val==[4,5])
      
          if grid_res_1 == 1  & strcmpi(gridModel,'GPT2w')
             grid_MFw    = getappdata(0,'gridFile'); %gridfile
             gridRES_MFw = getappdata(0,'gridRES');%grid resolution
         
          else 
              %*****READ GRID FILE
              grid_MFw    = readGPTgrid('gpt2_1w.mat','GPT2w',1);%Hydrostatic GRIDFILE
              gridRES_MFw = 1;%grid resolution
          end 
      
       else 
           %*****READ GRID FILE
           grid_MFw    = readGPTgrid('gpt2_1w.mat','GPT2w',1);%Hydrostatic GRIDFILE
           gridRES_MFw = 1;%grid resolution
       
       end %//if any(metPARA_val==[4,5])
       
   end %//if ~isempty(grid_MFh1) & ~isempty(gridRES_MFh1)
   
   if ((length(grid_MFw) > 10) | (length(grid_MFw)< 10))
      grid_MFw=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Combine model GRIDFILE
   end  
   
elseif strcmpi(get(hObject,'Enable'),'on') & strncmpi(MFw_model,'VMF1(5° x 5°)',12)
       
      if ~isempty(grid_MFw1) & ~isempty(gridRES_MFw1) 
       
         if gridRES_MFw1 == 5
           
             grid_MFw    = grid_MFw1;
             gridRES_MFw = gridRES_MFw1;
          
         else 
             %*****READ GRID FILE
             grid_MFw    = readGPTgrid('gpt2_5w.mat','GPT2w',5);%Hydrostatic GRIDFILE
             gridRES_MFw = 5;%grid resolution
         end 
       
      else 
           if any(metPARA_val==[4,5])
      
              if grid_res_5 == 1  & strcmpi(gridModel,'GPT2w')
                 grid_MFw    = getappdata(0,'gridFile'); %gridfile
                 gridRES_MFw = getappdata(0,'gridRES');%grid resolution
         
              else  
                  %*****READ GRID FILE
                  grid_MFw    = readGPTgrid('gpt2_5w.mat','GPT2w',5);%Hydrostatic GRIDFILE
                  gridRES_MFw = 5;%grid resolution
              end  
      
           else   
                %*****READ GRID FILE
                grid_MFw    = readGPTgrid('gpt2_5w.mat','GPT2w',5);%Hydrostatic GRIDFILE
                gridRES_MFw = 5;%grid resolution
       
           end  %//if any(metPARA_val==[4,5])
           
      end %//if ~isempty(grid_MFh1) & ~isempty(gridRES_MFh1)
      
      if ((length(grid_MFw) > 10) | (length(grid_MFw)< 10))
         grid_MFw=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Combine model GRIDFILE
      end    
            
elseif strcmpi(get(hObject,'Enable'),'on') & strncmpi(MFw_model,'VMF3(1° x 1°)',12)
   
      if ~isempty(grid_MFw1) & ~isempty(gridRES_MFw1) 
       
         if gridRES_MFw1 == 1
           
             grid_MFw    = grid_MFw1;
             gridRES_MFw = gridRES_MFw1;
          
         else 
             %*****READ GRID FILE
             grid_MFw    = readGPTgrid('gpt3_1.mat','GPT3',1);%Hydrostatic GRIDFILE
             gridRES_MFw = 1;%grid resolution
         end 
       
      else  
          if any(metPARA_val==[4,5])
      
             if grid_res_1 == 1  & strcmpi(gridModel,'GPT3')
                grid_MFw    = getappdata(0,'gridFile'); %gridfile
                gridRES_MFw = getappdata(0,'gridRES');%grid resolution
             else   
                 %*****READ GRID FILE
                 grid_MFw    = readGPTgrid('gpt3_1.mat','GPT3',1);%Hydrostatic GRIDFILE
                 gridRES_MFw = 1;%grid resolution
             end   
      
          else  
              %*****READ GRID FILE
              grid_MFw    = readGPTgrid('gpt3_1.mat','GPT3',1);%Hydrostatic GRIDFILE
              gridRES_MFw = 1;%grid resolution
       
          end %//if any(metPARA_val==[4,5]) 
          
      end %//if ~isempty(grid_MFh1) & ~isempty(gridRES_MFh1) 
      
      if length(grid_MFw)<14
         grid_MFw=readGPTgrid('gpt3_1.mat','GPT3',1);%Combine model GRIDFILE
      end    
   
elseif strcmpi(get(hObject,'Enable'),'on') & strncmpi(MFw_model,'VMF3(5° x 5°)',12)
       
      if ~isempty(grid_MFw1) & ~isempty(gridRES_MFw1) 
       
         if gridRES_MFw1 == 5
           
             grid_MFw    = grid_MFw1;
             gridRES_MFw = gridRES_MFw1;
          
         else  
             %*****READ GRID FILE
             grid_MFw    = readGPTgrid('gpt3_5.mat','GPT3',5);%Hydrostatic GRIDFILE
             gridRES_MFw = 5;%grid resolution
         end 
       
      else 
          
           if any(metPARA_val==[4,5])
      
             if grid_res_5 == 1  & strcmpi(gridModel,'GPT3')
                grid_MFw    = getappdata(0,'gridFile'); %gridfile
                gridRES_MFw = getappdata(0,'gridRES');%grid resolution
         
             else  
                 %*****READ GRID FILE
                 grid_MFw    = readGPTgrid('gpt3_5.mat','GPT3',5);%Hydrostatic GRIDFILE
                 gridRES_MFw = 5;%grid resolution
             end  
      
           else    
               %*****READ GRID FILE
               grid_MFw    = readGPTgrid('gpt3_5.mat','GPT3',5);%Hydrostatic GRIDFILE
               gridRES_MFw = 5;%grid resolution
       
           end %//if any(metPARA_val==[4,5])  
           
      end %//if ~isempty(grid_MFh1) & ~isempty(gridRES_MFh1)
      
      if length(grid_MFw)<14
         grid_MFw=readGPTgrid('gpt3_5.mat','GPT3',5);%Combine model GRIDFILE
      end    
            
end %//if strcmpi(get(hObject,'Enable'),'on') & strncmpi(MFh_model,'VMF1(1° x 1°)',12)

%SAVE GRID FILE & RESOLUTION 
if any([exist('grid_MFw','var') exist('gridRES_MFw','var')])
   setappdata(0,'grid_MFw',grid_MFw)
   setappdata(0,'gridRES_MFw',gridRES_MFw) 
end

%SAVE SELECTED WET MAPPING FUNCTION MODEL 
setappdata(0,'MFw_model',MFw_model)

%==========================END OF WET_mappingfunction_popm_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%                            --------------------------
%(3)**************************METEOROLOGICAL PARAMETERS
%                            --------------------------
% --- EXECUTES ON SELECTION CHANGE in pop_source_metpara.
function METparameters_popm_Callback(hObject,handles)
%--------------------------------------------------------------------------
global METpara
getconts = cellstr(get(hObject,'String'));%get popup menu contents as cell array
metPARA=getconts{get(hObject,'Value')};%get selected item from pop_lTropo_combine
metVAL=get(hObject,'Value');%GET SELECTED VALUE

%HID TOGGLE STATE OF GRID RESOLUTION TYPE (1 OR 5)
if any(metVAL==[1,2,3,6,7,8])
   set(handles.rb_grid_resolution_1,'Visible','off')
   set(handles.rb_grid_resolution_5,'Visible','off')
   set(handles.text_grid_resolution,'Visible','off')
   set(handles.pop_met_manual,'visible','off')
   set(handles.text_metmanual_source,'visible','off')
   
   if any(metVAL==[7,8])
      set(handles.pop_met_manual,'visible','on')
      set(handles.text_metmanual_source,'visible','on')
   end
   
else %MAKE TOGGLE STATE OF GRID RESOLUTION TYPE (1 OR 5) VISIBLE
     set(handles.rb_grid_resolution_1,'Visible','on')
     if get(handles.rb_grid_resolution_1,'value')==0 & get(handles.rb_grid_resolution_5,'value')==0
        set(handles.rb_grid_resolution_5,'Visible','on','Value',1,'ForegroundColor',[0,0.75,0.75])%set as default
     else
          set(handles.rb_grid_resolution_5,'Visible','on')
     end
     set(handles.text_grid_resolution,'Visible','on') 
     set(handles.pop_met_manual,'visible','off')
     set(handles.text_metmanual_source,'visible','off')
end

%GET TOGGLE GRID RESOLUTION TYPE
grid_res1=get(handles.rb_grid_resolution_1,'value');%for 1° resolution
grid_res5=get(handles.rb_grid_resolution_5,'value');%for 5° resolution

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%*******FOR READING GLOBAL PRESSURE & TEMPERATURE(GPT) GRIDFILES
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%GET PREVIOUS SELECTED MET SOURCE MODEL
metPARA1=getappdata(0,'METpara');

%GET STORED GRIDFILE & RESOLUTION FOR MET PARAMETERS
grid_met1=getappdata(0,'grid_met');%Previous comnine gridfile
grid_res_met1=getappdata(0,'grid_res_met');%Previous comnine grid resolution

%GET COMBINE TROPO MODEL GRIDFILE(grid_c)
tropoModel1=getappdata(0,'Tropo_Cmodel');%combine topospheric model
grid_c1=getappdata(0,'grid_c');%Previous comnine gridfile
grid_res_c1=getappdata(0,'grid_res_c');%Previous comnine grid resolution

%GET STORED GRIDFILE(grid_h) FROM THIS POPUP MENU 
dryModel1=getappdata(0,'dryModel');%previous hydrostatic model
grid_h1=getappdata(0,'grid_h');   %previous gridfile
grid_res_h1=getappdata(0,'grid_res_h');%previous grid resolution

%GET STORED GRIDFILE(grid_w) FROM WET MODEL POPUP MENU 
wetModel1=getappdata(0,'wetModel');%previous wet model
grid_w1=getappdata(0,'grid_w');   %previous gridfile
grid_res_w1=getappdata(0,'grid_res_w');%previous grid resolution
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

if any([strcmpi(metPARA,'1.Standard Meteorological Parameters'),isequal(metVAL,1)])
   set(hObject,'TooltipString','Solution Using Standard Meteorological Parameters')%change Tooltip any time user select Model type 
   METpara='standard';
   
elseif any([strcmpi(metPARA,'2.Generate from GPT model'),isequal(metVAL,2)])
       set(hObject,'TooltipString','Solution Using Meteorological Parameters Generated from GPT model')%change Tooltip any time user select Model type 
       METpara='GPT';
       
       %****GET SUPPLEMENTARY FILES FOR SUPPLEMENTARY MET PARAMETERS
       %-------------------------------------------------------------------
       %GPT MODEL PROVIDES ONLY PRESSURE(P) & TEMPERATURE(T).WE NEED OTHER
       %MET PARAMETERS LIKE WATER VAPOUR PARTIAL PRESSURE(e) TO SUPPLEMENT
       %WHAT THE GPT MODEL PROVIDES----------------------------------------
       
       %FIRST CHECK IF ALREADY STORED GRID FILE IS EMPTY([]) OR NOT
       if ~isempty(getappdata(0,'grid_met_s'))
          
          %GET & STORE GRID VALUES,GRID RESOLUTION & GPT MODEL TYPE(GPT2w/GPT3) 
          setappdata(0,'grid_met_s',getappdata(0,'grid_met_s'))
          setappdata(0,'grid_res_met_s',getappdata(0,'grid_res_met_s'))
          setappdata(0,'grid_METmodel',getappdata(0,'grid_METmodel'))
          
       else   
           %READ FROM GRID FILE(using GPT2w model)
           grid_met_ss = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
           grid_res_met_ss = 5;%GRID RESOLUTION 
           grid_MET_model = 'GPT2w';
       
           if isempty(grid_met_ss) %if GPT2w model grid is unavailable
              %READ FROM GRID FILE(using GPT3 model)
              grid_met_ss=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE 
              grid_MET_model = 'GPT3';
           end   
   
           setappdata(0,'grid_met_s',grid_met_ss) %SAVE GRIDFILE FOR MET PARAMETER
           setappdata(0,'grid_res_met_s',grid_res_met_ss)%SAVE GRID RESOLUTION FOR MET PARAMETER
           setappdata(0,'grid_METmodel',grid_MET_model)%SAVE GRID MODEL
           
       end %if ~isempty(getappdata(0,'grid_met_s'))
       
elseif any([strcmpi(metPARA,'3.Generate from GPT2 model'),isequal(metVAL,3)])
       set(hObject,'TooltipString','Solution Using Meteorological Parameters Generated from GPT2 model')%change Tooltip any time user select Model type 
       METpara='GPT2';
       
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       %FIRST, CHECK FROM PREVIOUS STORED FILE FOR MET PARAMETERS
       if ~isempty(grid_met1) & grid_res_met1 == 5 & ~isempty(metPARA1) %#ok<*AND2>
           
         if strcmpi(metPARA1,'GPT2')% IF MET PARAMETER SOURCE IS GPT2
            grid_met     = grid_met1;
            grid_res_met = 5;
            
         else  
             %GET GRIDFILE FROM TROPO MODELs
             if ~isempty(grid_c1) & grid_res_c1 == 5
           
               if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT2 (5° x 5°)',14) strncmpi(tropoModel1,'GPT2 (5° x 5°)',14)])
                  grid_met     = grid_c1;
                  grid_res_met = 5;
               end 
           
             elseif ~isempty(grid_h1) & grid_res_h1 == 5
                    if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT2 (5° x 5°)',14) strncmpi(dryModel1,'GPT2 (5° x 5°)',14)])
                       grid_met     = grid_h1;
                       grid_res_met = 5;
                    end 
             
             elseif ~isempty(grid_w1) & grid_res_w1 == 5
                    if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT2 (5° x 5°)',14) strncmpi(wetModel1,'GPT2 (5° x 5°)',14)])
                       grid_met     = grid_w1;
                       grid_res_met = 5;
                    end 
              
             else %IF ALL FAIL,READ GRID
                  grid_met=readGPTgrid('gpt2_5.mat','GPT2',5);%GRIDFILE
                  grid_res_met = 5;
             end 
         end
         
         if any ([~exist('grid_met','var') ~exist('grid_res_met','var')])
            grid_met=readGPTgrid('gpt2_5.mat','GPT2',5);%GRIDFILE
            grid_res_met = 5; 
         end
         
       else
           %GET GRIDFILE FROM TROPO MODELs
           if ~isempty(grid_c1) & grid_res_c1 == 5
               
              if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT2 (5° x 5°)',14) strncmpi(tropoModel1,'GPT2 (5° x 5°)',14)])
                 grid_met     = grid_c1;
                 grid_res_met = 5;
              end  
           
           elseif ~isempty(grid_h1) & grid_res_h1 == 5
                  if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT2 (5° x 5°)',14) strncmpi(dryModel1,'GPT2 (5° x 5°)',14)])
                     grid_met     = grid_h1;
                     grid_res_met = 5;
                  end  
             
           elseif ~isempty(grid_w1) & grid_res_w1 == 5
                  if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT2 (5° x 5°)',14) strncmpi(wetModel1,'GPT2 (5° x 5°)',14)])
                     grid_met     = grid_w1;
                     grid_res_met = 5;
                  end  
              
           else %IF ALL FAIL,READ GRID
                grid_met=readGPTgrid('gpt2_5.mat','GPT2',5);%GRIDFILE
                grid_res_met = 5;
           end
           
           if any ([~exist('grid_met','var') ~exist('grid_res_met','var')])
            grid_met=readGPTgrid('gpt2_5.mat','GPT2',5);%GRIDFILE
            grid_res_met = 5 ;
           end
       end
       
       if length(grid_met) > 8 
          grid_met=readGPTgrid('gpt2_5.mat','GPT2',5);%Combine model GRIDFILE
       end   
%--------------------------------------------------------------------------       
%SAVE GRID FILE & RESOLUTION      
setappdata(0,'grid_met',grid_met); %gridfile
setappdata(0,'grid_res_met',grid_res_met);%grid resolution         
%--------------------------------------------------------------------------
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 
elseif any([strcmpi(metPARA,'4.Generate from GPT2w model'),isequal(metVAL,4)])
       set(hObject,'TooltipString','Solution Using Meteorological Parameters Generated from GPT2w model')%change Tooltip any time user select Model type 
       METpara='GPT2w';
       
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       %FIRST, CHECK FROM PREVIOUS STORED FILE FOR MET PARAMETERS
       if grid_res5 == 1 & grid_res1 == 0
           
          if ~isempty(grid_met1) & grid_res_met1 == 5 & ~isempty(metPARA1)
           
             if strcmpi(metPARA1,'GPT2w')% IF MET PARAMETER SOURCE IS GPT2
                grid_met     = grid_met1;
                grid_res_met = 5;
            
             else   
                 %GET GRIDFILE FROM TROPO MODELs
                 if ~isempty(grid_c1) & grid_res_c1 == 5
           
                    if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT2w (5° x 5°)',14) strncmpi(tropoModel1,'GPT2w (5° x 5°)',14)])
                       grid_met     = grid_c1;
                       grid_res_met = 5; 
                    end  
           
                 elseif ~isempty(grid_h1) & grid_res_h1 == 5
                        if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT2w (5° x 5°)',14) strncmpi(dryModel1,'GPT2w (5° x 5°)',14)])
                           grid_met     = grid_h1;
                           grid_res_met = 5;
                        end  
             
             elseif ~isempty(grid_w1) & grid_res_w1 == 5
                    if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT2w (5° x 5°)',14) strncmpi(wetModel1,'GPT2w (5° x 5°)',4)])
                       grid_met     = grid_w1;
                       grid_res_met = 5;
                    end 
              
                 else  %IF ALL FAIL,READ GRID
                      grid_met=readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
                      grid_res_met = 5;%GRID RESOLUTION
                 end  
             end 
         
            if any ([~exist('grid_met','var') ~exist('grid_res_met','var')])
               grid_met=readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
               grid_res_met = 5;%GRID RESOLUTION 
               
            elseif any ([exist('grid_met','var') exist('grid_res_met','var')]) & isequal(grid_res_met,1)
                   grid_met=readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
                   grid_res_met = 5;%GRID RESOLUTION
            end  
         
          else 
              %GET GRIDFILE FROM TROPO MODELs
              if ~isempty(grid_c1) & grid_res_c1 == 5
           
                if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT2w (5° x 5°)',14) strncmpi(tropoModel1,'GPT2w (5° x 5°)',14)])
                   grid_met     = grid_c1;
                   grid_res_met = 5;
                end   
           
              elseif  ~isempty(grid_h1) & grid_res_h1 == 5
                     if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT2w (5° x 5°)',14) strncmpi(dryModel1,'GPT2w (5° x 5°)',14)])
                        grid_met     = grid_h1;
                        grid_res_met = 5;
                     end   
             
              elseif  ~isempty(grid_w1) & grid_res_w1 == 5
                      if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT2w (5° x 5°)',14) strncmpi(wetModel1,'GPT2w (5° x 5°)',14)])
                         grid_met     = grid_w1;
                         grid_res_met = 5;
                      end   
              
              else %IF ALL FAIL,READ GRID
                   grid_met=readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
                   grid_res_met = 5;
              end 
           
              if any ([~exist('grid_met','var') ~exist('grid_res_met','var')])
                 grid_met=readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
                 grid_res_met = 5;%GRID RESOLUTION
                 
              elseif any ([exist('grid_met','var') exist('grid_res_met','var')]) & isequal(grid_res_met,1)
                     grid_met=readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
                     grid_res_met = 5;%GRID RESOLUTION
              end   
          end
          
          if any([length(grid_met) > 10,length(grid_met)< 10])
             grid_met=readGPTgrid('gpt2_5w.mat','GPT2w',5);%Combine model GRIDFILE
          end     
          
       elseif grid_res5 == 0 & grid_res1 == 1
               
              if ~isempty(grid_met1) & grid_res_met1 == 1 & ~isempty(metPARA1)
           
                 if strcmpi(metPARA1,'GPT2w')% IF MET PARAMETER SOURCE IS GPT2
                    grid_met     = grid_met1;
                    grid_res_met = 1;
            
                 else     
                     %GET GRIDFILE FROM TROPO MODELs
                     if ~isempty(grid_c1) & grid_res_c1 == 1
           
                        if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT2w (1° x 1°)',14) strncmpi(tropoModel1,'GPT2w (1° x 1°)',4)])
                           grid_met     = grid_c1;
                           grid_res_met = 1; 
                        end   
           
                 elseif ~isempty(grid_h1) & grid_res_h1 == 1
                        if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strcmpi(dryModel1,'GPT2w (1° x 1°)') strncmpi(dryModel1,'GPT2w (1° x 1°)',14)])
                           grid_met     = grid_h1;
                           grid_res_met = 1;
                        end  
             
             elseif ~isempty(grid_w1) & grid_res_w1 == 1
                    if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT2w (1° x 1°)',14) strncmpi(wetModel1,'GPT2w (1° x 1°)',4)])
                       grid_met     = grid_w1;
                       grid_res_met = 1;
                    end 
              
                 else  %IF ALL FAIL,READ GRID
                      grid_met=readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
                      grid_res_met = 1;%GRID RESOLUTION
                     end   
                 end  
         
            if any ([~exist('grid_met','var') ~exist('grid_res_met','var')])
               grid_met=readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
               grid_res_met = 1;%GRID RESOLUTION 
               
            elseif  any ([exist('grid_met','var') exist('grid_res_met','var')]) & isequal(grid_res_met,5)
                     grid_met=readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
                     grid_res_met = 1;%GRID RESOLUTION
            end  
         
          else 
              %GET GRIDFILE FROM TROPO MODELs
              if ~isempty(grid_c1) & grid_res_c1 == 1
           
                if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT2w (1° x 1°)',14) strncmpi(tropoModel1,'GPT2w (1° x 1°)',14)])
                   grid_met     = grid_c1;
                   grid_res_met = 1;
                end   
           
              elseif  ~isempty(grid_h1) & grid_res_h1 == 1
                     if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT2w (1° x 1°)',14) strncmpi(dryModel1,'GPT2w (1° x 1°)',14)])
                        grid_met     = grid_h1;
                        grid_res_met = 1;
                     end   
             
              elseif  ~isempty(grid_w1) & grid_res_w1 == 1
                      if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT2w (1° x 1°)',14) strncmpi(wetModel1,'GPT2w (1° x 1°)',14)])
                         grid_met     = grid_w1;
                         grid_res_met = 1;
                      end   
              
              else %IF ALL FAIL,READ GRID
                   grid_met=readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
                   grid_res_met = 1;
              end 
           
              if any ([~exist('grid_met','var') ~exist('grid_res_met','var')])
                 grid_met=readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
                 grid_res_met = 1;%GRID RESOLUTION 
                 
              elseif any ([exist('grid_met','var') exist('grid_res_met','var')]) & isequal(grid_res_met,5)
                     grid_met=readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
                     grid_res_met = 1;%GRID RESOLUTION
              end
              
              end
              
              if any([length(grid_met) > 10,length(grid_met)< 10])
                 grid_met=readGPTgrid('gpt2_1w.mat','GPT2w',1);%Combine model GRIDFILE
              end     

       end %if grid_res5 == 1 & grid_res1 == 0
%--------------------------------------------------------------------------       
%SAVE GRID FILE & RESOLUTION      
setappdata(0,'grid_met',grid_met); %gridfile
setappdata(0,'grid_res_met',grid_res_met);%grid resolution        
%--------------------------------------------------------------------------
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       
elseif any([strcmpi(metPARA,'5.Generate from GPT3 model'),isequal(metVAL,5)])
       set(hObject,'TooltipString','Solution Using Meteorological Parameters Generated from GPT3 model')%change Tooltip any time user select Model type 
       METpara='GPT3';
       
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 
       %CHECK IF MODEL HAS BEEN RAN ALREADY, HERE OR ELSEWHERE AND RETRIEVE
       %FILE. THIS IS TO HELP AVOID READING SAME FILE OVER & OVER AGAIN
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       %FIRST, CHECK FROM PREVIOUS STORED FILE FOR MET PARAMETERS
      if grid_res5 == 1 & grid_res1 == 0
       
       if ~isempty(grid_met1) & grid_res_met1 == 5 & ~isempty(metPARA1)
           
         if strcmpi(metPARA1,'GPT3')% IF MET PARAMETER SOURCE IS GPT2
            grid_met     = grid_met1;
            grid_res_met = 5;
            
         else  
             %GET GRIDFILE FROM TROPO MODELs
             if ~isempty(grid_c1) & grid_res_c1 == 5
           
               if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT3 (5° x 5°)',14) strncmpi(tropoModel1,'GPT3 (5° x 5°)',14)])
                  grid_met     = grid_c1;
                  grid_res_met = 5;
               end 
           
             elseif ~isempty(grid_h1) & grid_res_h1 == 5
                    if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT3 (5° x 5°)',14) strncmpi(dryModel1,'GPT3 (5° x 5°)',14)])
                       grid_met     = grid_h1;
                       grid_res_met = 5;
                    end 
             
             elseif ~isempty(grid_w1) & grid_res_w1 == 5
                    if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT3 (5° x 5°)',14) strncmpi(wetModel1,'GPT3 (5° x 5°)',14)])
                       grid_met     = grid_w1;
                       grid_res_met = 5;
                    end 
              
             else %IF ALL FAIL,READ GRID
                  grid_met=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE
                  grid_res_met = 5;
             end 
         end
         
         if any ([~exist('grid_met','var') ~exist('grid_res_met','var')])
            grid_met=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE
            grid_res_met = 5;%GRID RESOLUTION 
            
         elseif any ([exist('grid_met','var') exist('grid_res_met','var')]) & isequal(grid_res_met,1)
                grid_met=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE
                grid_res_met = 5;%GRID RESOLUTION 
         end 
         
       else
           %GET GRIDFILE FROM TROPO MODELs
           if ~isempty(grid_c1) & grid_res_c1 == 5
           
              if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT3 (5° x 5°)',14) strncmpi(tropoModel1,'GPT3 (5° x 5°)',14)])
                 grid_met     = grid_c1;
                 grid_res_met = 5;
              end  
           
           elseif ~isempty(grid_h1) & grid_res_h1 == 5
                  if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT3 (5° x 5°)',14) strncmpi(dryModel1,'GPT3 (5° x 5°)',14)])
                     grid_met     = grid_h1;
                     grid_res_met = 5;
                  end  
             
           elseif ~isempty(grid_w1) & grid_res_w1 == 5
                  if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT3 (5° x 5°)',14) strncmpi(wetModel1,'GPT3 (5° x 5°)',14)])
                     grid_met     = grid_w1;
                     grid_res_met = 5;
                  end  
              
           else %IF ALL FAIL,READ GRID
                grid_met=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE
                grid_res_met = 5;
           end
           
           if any ([~exist('grid_met','var') ~exist('grid_res_met','var')])
              grid_met=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE
              grid_res_met = 5;%GRID RESOLUTION 
              
           elseif any ([exist('grid_met','var') exist('grid_res_met','var')]) & isequal(grid_res_met,1)
                  grid_met=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE
                  grid_res_met = 5;%GRID RESOLUTION 
           end 
         
       end %//if ~isempty(grid_met1) & grid_res_met1 == 5 & ~isempty(metPARA1)
       
       if length(grid_met)<14
          grid_met=readGPTgrid('gpt3_5.mat','GPT3',5);%Combine model GRIDFILE
       end    
       
      elseif grid_res5 == 0 & grid_res1 == 1
          
             if ~isempty(grid_met1) & grid_res_met1 == 1 & ~isempty(metPARA1)
           
                 if strcmpi(metPARA1,'GPT3')% IF MET PARAMETER SOURCE IS GPT2
                    grid_met     = grid_met1;
                    grid_res_met = 1;
            
                 else     
                     %GET GRIDFILE FROM TROPO MODELs
                     if ~isempty(grid_c1) & grid_res_c1 == 1
           
                        if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT3 (1° x 1°)',14) strncmpi(tropoModel1,'GPT3 (1° x 1°)',14)])
                           grid_met     = grid_c1;
                           grid_res_met = 1; 
                        end   
           
                 elseif ~isempty(grid_h1) & grid_res_h1 == 1
                        if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT3 (1° x 1°)',14) strncmpi(dryModel1,'GPT3 (1° x 1°)',14)])
                           grid_met     = grid_h1;
                           grid_res_met = 1;
                        end  
             
             elseif ~isempty(grid_w1) & grid_res_w1 == 1
                    if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT3 (1° x 1°)',14) strncmpi(wetModel1,'GPT3 (1° x 1°)',14)])
                       grid_met     = grid_w1;
                       grid_res_met = 1;
                    end 
              
                 else  %IF ALL FAIL,READ GRID
                      grid_met=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE
                      grid_res_met = 1;%GRID RESOLUTION
                     end   
                 end  
         
            if any ([~exist('grid_met','var') ~exist('grid_res_met','var')])
               grid_met=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE
               grid_res_met = 1;%GRID RESOLUTION 
               
            elseif any ([exist('grid_met','var') exist('grid_res_met','var')]) & isequal(grid_res_met,5)
                   grid_met=readGPTgrid('gpt3_5.mat','GPT3',1);%GRIDFILE
                   grid_res_met = 1;%GRID RESOLUTION 
            end  
         
          else 
              %GET GRIDFILE FROM TROPO MODELs
              if ~isempty(grid_c1) & grid_res_c1 == 1
           
                if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT3 (1° x 1°)',14) strncmpi(tropoModel1,'GPT3 (1° x 1°)',14)])
                   grid_met     = grid_c1;
                   grid_res_met = 1;
                end   
           
              elseif  ~isempty(grid_h1) & grid_res_h1 == 1
                     if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT3 (1° x 1°)',14) strncmpi(dryModel1,'GPT3 (1° x 1°)',14)])
                        grid_met     = grid_h1;
                        grid_res_met = 1;
                     end   
             
              elseif  ~isempty(grid_w1) & grid_res_w1 == 1
                      if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT3 (1° x 1°)',14) strncmpi(wetModel1,'GPT3 (1° x 1°)',14)])
                         grid_met     = grid_w1;
                         grid_res_met = 1;
                      end   
              
              else %IF ALL FAIL,READ GRID
                   grid_met=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE
                   grid_res_met = 1;
              end 
           
              if any ([~exist('grid_met','var') ~exist('grid_res_met','var')])
                 grid_met=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE
                 grid_res_met = 1;%GRID RESOLUTION 
                 
              elseif any ([exist('grid_met','var') exist('grid_res_met','var')]) & isequal(grid_res_met,5)
                     grid_met=readGPTgrid('gpt3_5.mat','GPT3',1);%GRIDFILE
                     grid_res_met = 1;%GRID RESOLUTION 
              end 
              
             end %//if ~isempty(grid_met1) & grid_res_met1 == 1 & ~isempty(metPARA1)
             
           if length(grid_met)<14
              grid_met=readGPTgrid('gpt3_1.mat','GPT3',1);%Combine model GRIDFILE
           end   
          
      end %if grid_res5 == 1 & grid_res1 == 0
           
%--------------------------------------------------------------------------          
%SAVE GRID FILE & RESOLUTION      
setappdata(0,'grid_met',grid_met); %gridfile
setappdata(0,'grid_res_met',grid_res_met);%grid resolution  
%--------------------------------------------------------------------------
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       
elseif any([strcmpi(metPARA,'6.Generate from UNB3m model'),isequal(metVAL,6)])
       set(hObject,'TooltipString','Solution Using Meteorological Parameters Generated from UNB3m model')%change Tooltip any time user select Model type 
       METpara='UNB3m';
       
elseif any([strcmpi(metPARA,'7.Enter Meteorological Parameters'),isequal(metVAL,7)])
       set(hObject,'TooltipString','Solution Using User input/defined Meteorological Parameters')%change Tooltip any time user select Model type 
       
       %CREATE A UITABLE FOR USER INPUT
       metpara_table()
       METpara=getappdata(0,'metPARA_table');
       
elseif any([strcmpi(metPARA,'8.Import from csv file [T P RH]'),isequal(metVAL,8)])
       set(hObject,'TooltipString','Solution Using User input/defined Meteorological Parameters')%change Tooltip any time user select Model type 
       
       %OPEN A DIALOGUE FOR USER FILE IMPORT
       importMETPARA_csv();
       METpara=getappdata(0,'metPARA_csv');
       
elseif any([strcmpi(metPARA,'9.MET file'),isequal(metVAL,9)])
       set(hObject,'TooltipString','Solution Using Meteorological Parameters Observed at GNSS site')%change Tooltip any time user select Model type 
       METpara='MET file';
       
       %USE GPT2W 1° x 1° MET DATA AS ALTERNATIVE IF 'MET file' isempty([]) OR NAN
       
       %GET PREVIOUS SAVED GRID VALUES & GRID RESOLUTION 
       if ~isempty(getappdata(0,'grid_METfile'))%IF STORED GRID VALUES IS NOT EMPTY
          grid_METfile = getappdata(0,'grid_METfile');
          grid_res_METfile = getappdata(0,'grid_res_METfile');
          
       else %*****READ GRID FILE
            grid_METfile =readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
            grid_res_METfile = 1;%GRID RESOLUTION
            
       end
       
            %SAVE GRID VALUES & GRID RESOLUTION
            setappdata(0,'grid_METfile',grid_METfile)
            setappdata(0,'grid_res_METfile',grid_res_METfile)                          
end
    
%SAVE SELECTED COMBINE TROPO MODEL
setappdata(0,'METpara',METpara)

%CHECK IF ANY OF THE GPT MODEL(GPT2,GPT2w & GPT3) IS SELECTED,AND ALLOW
%USER TO INDICATE WHETHER TO APPLY TIME VARIATION OF NOT
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%NOTE:
%    GLOBAL PRESSURE & TEMPERATURE(GPT) MODELS SUCH AS GPT2,GPT2w & GPT3 
%    COMPUTES METEOROLOGICAL PARAMETERS EITHER IN STATIC MODE OR IN TIME
%    VARYING MODE:
%    case 1: no time variation but static quantities
%    case 0: with time variation (annual and semiannual terms)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
 if any(strcmpi(METpara,{'GPT2','GPT2w','GPT3'}))
    
    if ~isempty(grid_met)
       TIMEvariation()
       uiwait(gcf)
       setappdata(0,'Timevar_met',getappdata(0,'Timevar'))
    end
 
 end
%--------------------------------------------------------------------------
 
%                    -------------------------------------------
%(3.1).**************SELECTION OF GRID RESOLUTION(RADIO BUTTON)
%                    -------------------------------------------

%                ------------------       
%(3.1.1)******** 1° GRID RESOLUTION
%                ------------------
% --- EXECUTES ON BUTTON PRESS IN rb_grid_resolution_1.
function GridResolution_1_rb_Callback(hObject,handles)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%*******FOR READING GLOBAL PRESSURE & TEMPERATURE(GPT) GRIDFILES
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%GET SOURCE OF MET PARAMETERS popup menu contents
getconts = cellstr(get(handles.pop_source_metpara,'String'));%get popup menu contents as cell array
metPARA=getconts{get(handles.pop_source_metpara,'Value')};%get selected item from pop_lTropo_combine
metVAL=get(handles.pop_source_metpara,'Value');%GET SELECTED VALUE

%GET PREVIOUS SELECTED MET SOURCE MODEL
metPARA1=getappdata(0,'METpara');%FROM SOURCE OF MET PARAMETER POPUP MENU

%GET STORED GRIDFILE & RESOLUTION FOR MET PARAMETERS
grid_met1=getappdata(0,'grid_met');%Previous comnine gridfile
grid_res_met1=getappdata(0,'grid_res_met');%Previous comnine grid resolution

grid_met_rb1_GPT2W=getappdata(0,'grid_met_rb1_gpt2w');%Previous gridfile from 1° grid resolution radio button using gpt2w model
grid_res_met_rb1_GPT2W=getappdata(0,'grid_res_met_rb1_gpt2w');%Previous grid resolution from 1° grid resolution radio button using gpt2w model

grid_met_rb1_GPT3=getappdata(0,'grid_met_rb1_gpt3');%Previous gridfile from 1° grid resolution radio button using gpt3 model
grid_res_met_rb1_GPT3=getappdata(0,'grid_res_met_rb1_gpt3');%Previous grid resolution from 1° grid resolution radio button using gpt3 model

%GET COMBINE TROPO MODEL GRIDFILE(grid_c)
tropoModel1=getappdata(0,'Tropo_Cmodel');%combine topospheric model
grid_c1=getappdata(0,'grid_c');%Previous comnine gridfile
grid_res_c1=getappdata(0,'grid_res_c');%Previous comnine grid resolution

%GET STORED GRIDFILE(grid_h) FROM THIS POPUP MENU 
dryModel1=getappdata(0,'dryModel');%previous hydrostatic model
grid_h1=getappdata(0,'grid_h');   %previous gridfile
grid_res_h1=getappdata(0,'grid_res_h');%previous grid resolution

%GET STORED GRIDFILE(grid_w) FROM WET MODEL POPUP MENU 
wetModel1=getappdata(0,'wetModel');%previous wet model
grid_w1=getappdata(0,'grid_w');   %previous gridfile
grid_res_w1=getappdata(0,'grid_res_w');%previous grid resolution
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%SET 1° GRID RESOLUTION RADIO BUTTON TO 1
set(hObject,'Value',1,'ForegroundColor',[0,0.75,0.75])

%SET 5° GRID RESOLUTION RADIO BUTTON TO 0
set(handles.rb_grid_resolution_5,'Value',0,'ForegroundColor','black')

if get(hObject,'Value') ==1 %IF USER SELECTS 1° resolution radio button
    
   %CHECK IF SELECTED SOURCE OF MET PARAMETERS IS GPT2w MODEL    
   if any([strcmpi(metPARA,'4.Generate from GPT2w model') isequal(metVAL,4)] ) 

      if ~isempty(grid_met1) & grid_res_met1 == 1 & ~isempty(metPARA1)
          
         if strcmpi(metPARA1,'GPT2w')% IF MET PARAMETER SOURCE IS GPT2w
            grid_met_rb1_gpt2w     = grid_met1; %Gridfile for 1° resolution radio button
            grid_res_met_rb1_gpt2w = 1;  %Grid resolution for 1° resolution radio button 
         else        
             %GET GRIDFILE FROM COMBINE(DRY+WET)TROPO MODELs
             if ~isempty(grid_c1) & grid_res_c1 == 1
           
                if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT2w (1° x 1°)',14) strncmpi(tropoModel1,'GPT2w (1° x 1°)',14)])
                   grid_met_rb1_gpt2w     = grid_c1;%Gridfile for 1° resolution radio button
                   grid_res_met_rb1_gpt2w = 1; %Grid resolution for 1° resolution radio button 
                end     
         
         %GET GRIDFILE FROM HYDROSTATIC TROPO MODELs           
         elseif ~isempty(grid_h1) & grid_res_h1 == 1
                if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT2w (1° x 1°)',14) strncmpi(dryModel1,'GPT2w (1° x 1°)',14)])
                   grid_met_rb1_gpt2w     = grid_h1;
                   grid_res_met_rb1_gpt2w = 1;%Grid resolution for 1° resolution radio button 
                end   
         
         %GET GRIDFILE FROM WET TROPO MODELs        
         elseif ~isempty(grid_w1) & grid_res_w1 == 1
                if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT2w (1° x 1°)',14) strncmpi(wetModel1,'GPT2w (1° x 1°)',14)])
                  grid_met_rb1_gpt2w     = grid_w1;
                  grid_res_met_rb1_gpt2w = 1;
                end  
              
         else  %IF ALL FAIL,READ GRID
               grid_met_rb1_gpt2w=readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
               grid_res_met_rb1_gpt2w = 1;%GRID RESOLUTION
             end     
         end   
      
     %IF ALL OPTIONS FAIL
     if any ([~exist('grid_met_rb1_gpt2w','var') ~exist('grid_res_met_rb1_gpt2w','var')])
        
        %READ FROM GRID FILE 
        grid_met_rb1_gpt2w=readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
        grid_res_met_rb1_gpt2w = 1;%GRID RESOLUTION 
        
     elseif  any ([exist('grid_met_rb1_gpt2w','var') exist('grid_res_met_rb1_gpt2w','var')]) & any ([isempty(grid_met_rb1_gpt2w) isempty(grid_res_met_rb1_gpt2w)])
             
             %READ FROM GRID FILE 
             grid_met_rb1_gpt2w=readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
             grid_res_met_rb1_gpt2w = 1;%GRID RESOLUTION
     end   
    
      elseif  ~isempty(grid_met_rb1_GPT2W) & ~isempty(grid_res_met_rb1_GPT2W) 
             grid_met_rb1_gpt2w     = grid_met_rb1_GPT2W;
             grid_res_met_rb1_gpt2w = 1;
          
      else    
          %GET GRIDFILE FROM COMBINE(DRY+WET)TROPO MODELs
          if ~isempty(grid_c1) & grid_res_c1 == 1
           
            if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT2w (1° x 1°)',14) strncmpi(tropoModel1,'GPT2w (1° x 1°)',14)])
             grid_met_rb1_gpt2w     = grid_c1;
             grid_res_met_rb1_gpt2w = 1;
            end     
        
       %GET GRIDFILE FROM HYDROSTATIC TROPO MODELs   
          elseif   ~isempty(grid_h1) & grid_res_h1 == 1
               if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT2w (1° x 1°)',14) strncmpi(dryModel1,'GPT2w (1° x 1°)',14)])
                  grid_met_rb1_gpt2w     = grid_h1;
                  grid_res_met_rb1_gpt2w = 1;
               end
               
       %GET GRIDFILE FROM WET TROPO MODELs      
       elseif  ~isempty(grid_w1) & grid_res_w1 == 1
               if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT2w (1° x 1°)',14) strncmpi(wetModel1,'GPT2w (1° x 1°)',14)])
                  grid_met_rb1_gpt2w     = grid_w1;
                  grid_res_met_rb1_gpt2w = 1;
               end    
              
       else %IF ALL FAIL,READ GRID
            grid_met_rb1_gpt2w=readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
            grid_res_met_rb1_gpt2w = 1;
          end  
       
      %IF ALL OPTIONS FAIL     
      if any ([~exist('grid_met_rb1_gpt2w','var') ~exist('grid_res_met_rb1_gpt2w','var')])
         
         %READ FROM GRID FILE  
         grid_met_rb1_gpt2w=readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
         grid_res_met_rb1_gpt2w = 1;%GRID RESOLUTION 
         
      elseif any ([exist('grid_met_rb1_gpt2w','var') exist('grid_res_met_rb1_gpt2w','var')]) & any ([isempty(grid_met_rb1_gpt2w) isempty(grid_res_met_rb1_gpt2w)])
             
             %READ FROM GRID FILE 
             grid_met_rb1_gpt2w=readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
             grid_res_met_rb1_gpt2w = 1;%GRID RESOLUTION     
      end 
      
      end %if ~isempty(grid_met1) & grid_res_met1 == 1 & ~isempty(metPARA1) 
      
      %CHECK AGAIN TO BE SURE
        if ~exist('grid_met_rb1_gpt2w','var')
           grid_met_rb1_gpt2w=readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
        
        elseif all([exist('grid_met_rb1_gpt2w','var'),any([length(grid_met_rb1_gpt2w) > 10,length(grid_met_rb1_gpt2w)< 10])])
               grid_met_rb1_gpt2w=readGPTgrid('gpt2_1w.mat','GPT2w',1);% GRIDFILE
        end     
        
        
      %SAVE GRID FILE & RESOLUTION      
      setappdata(0,'grid_met_rb1_gpt2w',grid_met_rb1_gpt2w); %gridfile
      setappdata(0,'grid_res_met_rb1_gpt2w',grid_res_met_rb1_gpt2w);%grid resolution 
      
      %CHECK IF SELECTED SOURCE OF MET PARAMETERS IS GPT3 MODEL   
   elseif any([strcmpi(metPARA,'5.Generate from GPT3 model') isequal(metVAL,5)] )
          
        if ~isempty(grid_met1) & grid_res_met1 == 1 & ~isempty(metPARA1)
           
         if strcmpi(metPARA1,'GPT3')% IF MET PARAMETER SOURCE IS GPT3
            grid_met_rb1_gpt3     = grid_met1;
            grid_res_met_rb1_gpt3 = 1;
            
         else  
             %GET GRIDFILE FROM COMBINE(DRY+WET)TROPO MODELs
             if ~isempty(grid_c1) & grid_res_c1 == 1
           
               if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT3 (1° x 1°)',14) strncmpi(tropoModel1,'GPT3 (1° x 1°)',14)])
                  grid_met_rb1_gpt3     = grid_c1;
                  grid_res_met_rb1_gpt3 = 1;
               end 
           
             %GET GRIDFILE FROM HYDROSTATIC TROPO MODELs  
             elseif ~isempty(grid_h1) & grid_res_h1 == 1
                    if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT3 (1° x 1°)',14) strncmpi(dryModel1,'GPT3 (1° x 1°)',14)])
                       grid_met_rb1_gpt3     = grid_h1;
                       grid_res_met_rb1_gpt3 = 1;
                    end 
             
             %GET GRIDFILE FROM WET TROPO MODELs       
             elseif ~isempty(grid_w1) & grid_res_w1 == 1
                    if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT3 (1° x 1°)',14) strncmpi(wetModel1,'GPT3 (1° x 1°)',14)])
                       grid_met_rb1_gpt3     = grid_w1;
                       grid_res_met_rb1_gpt3 = 1;
                    end 
              
             else %IF ALL FAIL,READ GRID
                  grid_met_rb1_gpt3=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE
                  grid_res_met_rb1_gpt3 = 1;%GRID RESOLUTION
             end 
         end
         
         %IF ALL OPTIONS FAIL
         if any ([~exist('grid_met_rb1_gpt3','var') ~exist('grid_res_met_rb1_gpt3','var')])
             
            %READ FROM GRID FILE
            grid_met_rb1_gpt3=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE
            grid_res_met_rb1_gpt3 = 1;%GRID RESOLUTION 
            
         elseif any ([exist('grid_met_rb1_gpt3','var') exist('grid_res_met_rb1_gpt3','var')]) & any ([isempty(grid_met_rb1_gpt3) isempty(grid_res_met_rb1_gpt3)])
                
                %READ FROM GRID FILE
                grid_met_rb1_gpt3=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE
                grid_res_met_rb1_gpt3 = 1;%GRID RESOLUTION
         end 
        
        elseif  ~isempty(grid_met_rb1_GPT3) & ~isempty(grid_res_met_rb1_GPT3) 
                grid_met_rb1_gpt3     = grid_met_rb1_GPT3;
                grid_res_met_rb1_gpt3 = 1;      
        else  
           %GET GRIDFILE FROM COMBINE(DRY+WET)TROPO MODELs
           if ~isempty(grid_c1) & grid_res_c1 == 1
           
              if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT3 (1° x 1°)',14) strncmpi(tropoModel1,'GPT3 (1° x 1°)',14)])
                 grid_met_rb1_gpt3     = grid_c1;
                 grid_res_met_rb1_gpt3 = 1;
              end  
           
           %GET GRIDFILE FROM HYDROSTATIC TROPO MODELs   
           elseif ~isempty(grid_h1) & grid_res_h1 == 1
                  if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT3 (1° x 1°)',14) strncmpi(dryModel1,'GPT3 (1° x 1°)',14)])
                     grid_met_rb1_gpt3     = grid_h1;
                     grid_res_met_rb1_gpt3 = 1;
                  end  
           
           %GET GRIDFILE FROM WET TROPO MODELs       
           elseif ~isempty(grid_w1) & grid_res_w1 == 1
                  if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT3 (1° x 1°)',14) strncmpi(wetModel1,'GPT3 (1° x 1°)',14)])
                     grid_met_rb1_gpt3     = grid_w1;
                     grid_res_met_rb1_gpt3 = 1;
                  end  
              
           else %IF ALL FAIL,READ GRID
                grid_met_rb1_gpt3=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE
                grid_res_met_rb1_gpt3 = 1;
           end
           
           %IF ALL OPTIONS FAIL
           if any ([~exist('grid_met_rb1_gpt3','var') ~exist('grid_res_met_rb1_gpt3','var')])
               
              %READ FROM GRID FILE 
              grid_met_rb1_gpt3=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE
              grid_res_met_rb1_gpt3 = 1;%GRID RESOLUTION 
              
           elseif any ([exist('grid_met_rb1_gpt3','var') exist('grid_res_met_rb1_gpt3','var')]) & any ([isempty(grid_met_rb1_gpt3) isempty(grid_res_met_rb1_gpt3)])
                  
                  %READ FROM GRID FILE
                  grid_met_rb1_gpt3=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE
                  grid_res_met_rb1_gpt3 = 1;%GRID RESOLUTION
           end  
         
        end %if ~isempty(grid_met1) & grid_res_met1 == 1 & ~isempty(metPARA1) 
        
        %CHECK AGAIN TO BE SURE
        if ~exist('grid_met_rb1_gpt3','var')
           grid_met_rb1_gpt3=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE
        
        elseif all([exist('grid_met_rb1_gpt3','var'),length(grid_met_rb1_gpt3)<14])
               grid_met_rb1_gpt3=readGPTgrid('gpt3_1.mat','GPT3',1);%Combine model GRIDFILE
        end     
           
        %SAVE GRID FILE & RESOLUTION      
        setappdata(0,'grid_met_rb1_gpt3',grid_met_rb1_gpt3); %gridfile
        setappdata(0,'grid_res_met_rb1_gpt3',grid_res_met_rb1_gpt3);%grid resolution 
      
   end %if any([strcmpi(metPARA,'4.Generate from GPT2w model') isequal(metVAL,4)] )      
%-------------------------------------------------------------------------- 
%******SAVE GRID FILE & RESOLUTION
%IF MODEL IS GPT2w
if any([strcmpi(metPARA,'4.Generate from GPT2w model') isequal(metVAL,4)] )
   setappdata(0,'grid_met_rb1',grid_met_rb1_gpt2w); %gridfile
   setappdata(0,'grid_res_met_rb1',grid_res_met_rb1_gpt2w);%grid resolution
   setappdata(0,'grid_model','GPT2w');%grid model
   
   %*****ASSIGNMENT
   gridFile  = grid_met_rb1_gpt2w;
   gridRES   = grid_res_met_rb1_gpt2w;
   gridModel = 'GPT2w'; 
   
   %IF MODEL IS GPT3
elseif any([strcmpi(metPARA,'5.Generate from GPT3 model') isequal(metVAL,5)] )
       setappdata(0,'grid_met_rb1',grid_met_rb1_gpt3); %gridfile
       setappdata(0,'grid_res_met_rb1',grid_res_met_rb1_gpt3);%grid resolution
       setappdata(0,'grid_model','GPT3');%grid model
       
       %*****ASSIGNMENT
       gridFile  = grid_met_rb1_gpt3;
       gridRES   = grid_res_met_rb1_gpt3;
       gridModel = 'GPT3'; 
       
else
     gridFile  = [];
     gridRES   = [];
     gridModel = []; 
    
end

if any([exist('gridFile','var') exist('gridRES','var') exist('gridModel','var')])
   %SAVE TO BE CALLED LATER FOR MF POPUP MENUs
   setappdata(0,'gridFile',gridFile); %gridfile
   setappdata(0,'gridRES',gridRES);%grid resolution
   setappdata(0,'gridModel',gridModel);%grid model 
end

end  %if get(hObject,'Value') ==1

%==========================END OF GridResolution_1_rb_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%                   ------------------
%(3.1.2)*********** 5° GRID RESOLUTION 
%                   ------------------
% --- EXECUTES ON BUTTON PRESS IN rb_grid_resolution_1.
function GridResolution_5_rb_Callback(hObject,handles)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%*******FOR READING GLOBAL PRESSURE & TEMPERATURE(GPT) GRIDFILES
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%GET SOURCE OF MET PARAMETERS popup menu contents
getconts = cellstr(get(handles.pop_source_metpara,'String'));%get popup menu contents as cell array
metPARA=getconts{get(handles.pop_source_metpara,'Value')};%get selected item from pop_lTropo_combine
metVAL=get(handles.pop_source_metpara,'Value');%GET SELECTED VALUE

%GET PREVIOUS SELECTED MET SOURCE MODEL
metPARA1=getappdata(0,'METpara');

%GET STORED GRIDFILE & RESOLUTION FOR MET PARAMETERS
grid_met1=getappdata(0,'grid_met');%Previous comnine gridfile
grid_res_met1=getappdata(0,'grid_res_met');%Previous comnine grid resolution

grid_met_rb5_GPT2W=getappdata(0,'grid_met_rb5_gpt2w');%Previous gridfile from 5° grid resolution radio button using gpt2w model
grid_res_met_rb5_GPT2W=getappdata(0,'grid_res_met_rb5_gpt2w');%Previous grid resolution from 1° grid resolution radio button using gpt2w model

grid_met_rb5_GPT3=getappdata(0,'grid_met_rb5_gpt3');%Previous gridfile from 5° grid resolution radio button using gpt3 model
grid_res_met_rb5_GPT3=getappdata(0,'grid_res_met_rb5_gpt3');%Previous grid resolution from 1° grid resolution radio button using gpt3 model

%GET COMBINE TROPO MODEL GRIDFILE(grid_c)
tropoModel1=getappdata(0,'Tropo_Cmodel');%combine topospheric model
grid_c1=getappdata(0,'grid_c');%Previous comnine gridfile
grid_res_c1=getappdata(0,'grid_res_c');%Previous comnine grid resolution

%GET STORED GRIDFILE(grid_h) FROM THIS POPUP MENU 
dryModel1=getappdata(0,'dryModel');%previous hydrostatic model
grid_h1=getappdata(0,'grid_h');   %previous gridfile
grid_res_h1=getappdata(0,'grid_res_h');%previous grid resolution

%GET STORED GRIDFILE(grid_w) FROM WET MODEL POPUP MENU 
wetModel1=getappdata(0,'wetModel');%previous wet model
grid_w1=getappdata(0,'grid_w');   %previous gridfile
grid_res_w1=getappdata(0,'grid_res_w');%previous grid resolution
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%SET 5° GRID RESOLUTION RADIO BUTTON TO 1
set(hObject,'Value',1,'ForegroundColor',[0,0.75,0.75])

%SET 1° GRID RESOLUTION RADIO BUTTON TO 0
set(handles.rb_grid_resolution_1,'Value',0,'ForegroundColor','black')

if get(hObject,'Value') ==1
    
   %CHECK IF SELECTED SOURCE OF MET PARAMETERS IS GPT2w MODEL    
   if any([strcmpi(metPARA,'4.Generate from GPT2w model') isequal(metVAL,4)] )
    
      if ~isempty(grid_met1) & grid_res_met1 == 5 & ~isempty(metPARA1)
           
         if strcmpi(metPARA1,'GPT2w')% IF MET PARAMETER SOURCE IS GPT2w
            grid_met_rb5_gpt2w     = grid_met1; %Gridfile for 5° resolution radio button
            grid_res_met_rb5_gpt2w = 5;  %Grid resolution for 5° resolution radio button 
     else      
         %GET GRIDFILE FROM COMBINE(DRY+WET) TROPO MODELs
         if ~isempty(grid_c1) & grid_res_c1 == 5
           
            if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT2w (5° x 5°)',14) strncmpi(tropoModel1,'GPT2w (5° x 5°)',14)])
               grid_met_rb5_gpt2w     = grid_c1;%Gridfile for 5° resolution radio button
               grid_res_met_rb5_gpt2w = 5; %Grid resolution for 5° resolution radio button
            end    
           
         %GET GRIDFILE FROM HYDROSTATIC TROPO MODELs     
         elseif ~isempty(grid_h1) & grid_res_h1 == 5
                if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT2w (5° x 5°)',14) strncmpi(dryModel1,'GPT2w (5° x 5°)',14)])
                   grid_met_rb5_gpt2w     = grid_h1;
                   grid_res_met_rb5_gpt2w = 5;
                end   
         
         %GET GRIDFILE FROM WET TROPO MODELs       
         elseif ~isempty(grid_w1) & grid_res_w1 == 5
                if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT2w (5° x 5°)',14) strncmpi(wetModel1,'GPT2w (5° x 5°)',14)])
                  grid_met_rb5_gpt2w     = grid_w1;
                  grid_res_met_rb5_gpt2w = 5;
                end  
              
         else  %IF ALL FAIL,READ GRID
               grid_met_rb5_gpt2w=readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
               grid_res_met_rb5_gpt2w = 5;%GRID RESOLUTION
         end    
         end   
       
     %IF ALL OPTIONS FAIL
     if any ([~exist('grid_met_rb5_gpt2w','var') ~exist('grid_res_met_rb5_gpt2w','var')])
        
        %READ FROM GRID FILE 
        grid_met_rb5_gpt2w=readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
        grid_res_met_rb5_gpt2w = 5;%GRID RESOLUTION 
        
     elseif  any ([exist('grid_met_rb5_gpt2w','var') exist('grid_res_met_rb5_gpt2w','var')]) & any ([isempty(grid_met_rb5_gpt2w) isempty(grid_res_met_rb5_gpt2w)])
             
             %READ FROM GRID FILE 
             grid_met_rb5_gpt2w=readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
             grid_res_met_rb5_gpt2w = 5;%GRID RESOLUTION
     end 
     
      elseif ~isempty(grid_met_rb5_GPT2W) & ~isempty(grid_res_met_rb5_GPT2W) 
             grid_met_rb5_gpt2w     = grid_met_rb5_GPT2W;
             grid_res_met_rb5_gpt2w = 5;
         
  else  
       %GET GRIDFILE FROM COMBINE(DRY+WET)TROPO MODELs
       if ~isempty(grid_c1) & grid_res_c1 == 5
           
          if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT2w (5° x 5°)',14) strncmpi(tropoModel1,'GPT2w (5° x 5°)',14)])
             grid_met_rb5_gpt2w     = grid_c1;
             grid_res_met_rb5_gpt2w = 5;
          end    
       
       %GET GRIDFILE FROM HYDROSTATIC TROPO MODELs    
       elseif  ~isempty(grid_h1) & grid_res_h1 == 5
               if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT2w (5° x 5°)',14) strncmpi(dryModel1,'GPT2w (5° x 5°)',14)])
                  grid_met_rb5_gpt2w     = grid_h1;
                  grid_res_met_rb5_gpt2w = 5;
               end    
       
       %GET GRIDFILE FROM WET TROPO MODELs
       elseif  ~isempty(grid_w1) & grid_res_w1 == 5
               if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT2w (5° x 5°)',14) strncmpi(wetModel1,'GPT2w (5° x 5°)',14)])
                  grid_met_rb5_gpt2w     = grid_w1;
                  grid_res_met_rb5_gpt2w = 5;
               end    
              
       else %IF ALL FAIL,READ GRID
            grid_met_rb5_gpt2w=readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
            grid_res_met_rb5_gpt2w = 5;
       end   
       
      %IF ALL OPTIONS FAIL
      if any ([~exist('grid_met_rb5_gpt2w','var') ~exist('grid_res_met_rb5_gpt2w','var')])
        
        %READ FROM GRID FILE 
        grid_met_rb5_gpt2w=readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
        grid_res_met_rb5_gpt2w = 5;%GRID RESOLUTION 
        
     elseif  any ([exist('grid_met_rb5_gpt2w','var') exist('grid_res_met_rb5_gpt2w','var')]) & any ([isempty(grid_met_rb5_gpt2w) isempty(grid_res_met_rb5_gpt2w)])
             
             %READ FROM GRID FILE 
             grid_met_rb5_gpt2w=readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
             grid_res_met_rb5_gpt2w = 5;%GRID RESOLUTION 
      end 
      
      end %if ~isempty(grid_met1) & grid_res_met1 == 5 & ~isempty(metPARA1)
      
       %CHECK AGAIN TO BE SURE
       if ~exist('grid_met_rb5_gpt2w','var')
          grid_met_rb5_gpt2w=readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
        
       elseif all([exist('grid_met_rb5_gpt2w','var'),any([length(grid_met_rb5_gpt2w) > 10,length(grid_met_rb5_gpt2w)< 10])])
              grid_met_rb5_gpt2w=readGPTgrid('gpt2_5w.mat','GPT2w',5);% GRIDFILE
       end      
      
      %SAVE GRID FILE & RESOLUTION      
      setappdata(0,'grid_met_rb5_gpt2w',grid_met_rb5_gpt2w); %gridfile
      setappdata(0,'grid_res_met_rb5_gpt2w',grid_res_met_rb5_gpt2w);%grid resolution 
      
   %CHECK IF SELECTED SOURCE OF MET PARAMETERS IS GPT3 MODEL   
   elseif any([strcmpi(metPARA,'5.Generate from GPT3 model') isequal(metVAL,5)] )
       
        if ~isempty(grid_met1) & grid_res_met1 == 5 & ~isempty(metPARA1)
           
         if strcmpi(metPARA1,'GPT3')% IF MET PARAMETER SOURCE IS GPT3
            grid_met_rb5_gpt3     = grid_met1;
            grid_res_met_rb5_gpt3 = 5;
            
         else  
             %GET GRIDFILE FROM COMBINE(DRY+WET)TROPO MODELs
             if ~isempty(grid_c1) & grid_res_c1 == 5
           
               if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT3 (5° x 5°)',14) strncmpi(tropoModel1,'GPT3 (5° x 5°)',14)])
                  grid_met_rb5_gpt3     = grid_c1;
                  grid_res_met_rb5_gpt3 = 5;
               end 
           
             %GET GRIDFILE FROM HYDROSTATIC TROPO MODELs  
             elseif ~isempty(grid_h1) & grid_res_h1 == 5
                    if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT3 (5° x 5°)',14) strncmpi(dryModel1,'GPT3 (5° x 5°)',14)])
                       grid_met_rb5_gpt3     = grid_h1;
                       grid_res_met_rb5_gpt3 = 5;
                    end 
             
             %GET GRIDFILE FROM WET TROPO MODELs       
             elseif ~isempty(grid_w1) & grid_res_w1 == 5
                    if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT3 (5° x 5°)',14) strncmpi(wetModel1,'GPT3 (5° x 5°)',14)])
                       grid_met_rb5_gpt3     = grid_w1;
                       grid_res_met_rb5_gpt3 = 5;
                    end 
              
             else %IF ALL FAIL,READ GRID
                  grid_met_rb5_gpt3=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE
                  grid_res_met_rb5_gpt3 = 5;
             end 
         end
         
         %IF ALL OPTIONS FAIL
         if any ([~exist('grid_met_rb5_gpt3','var') ~exist('grid_res_met_rb5_gpt3','var')])
             
            %READ FROM GRID FILE
            grid_met_rb5_gpt3=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE
            grid_res_met_rb5_gpt3 = 5;%GRID RESOLUTION 
            
         elseif any ([exist('grid_met_rb5_gpt3','var') exist('grid_res_met_rb5_gpt3','var')]) & any ([isempty(grid_met_rb5_gpt3) isempty(grid_res_met_rb5_gpt3)])
                
                %READ FROM GRID FILE
                grid_met_rb5_gpt3=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE
                grid_res_met_rb5_gpt3 = 5;%GRID RESOLUTION
         end
         
        elseif  ~isempty(grid_met_rb5_GPT3) & ~isempty(grid_res_met_rb5_GPT3) 
                grid_met_rb5_gpt3     = grid_met_rb5_GPT3;
                grid_res_met_rb5_gpt3 = 5;      
         
       else
           %GET GRIDFILE FROM COMBINE(DRY+WET)TROPO MODELs
           if ~isempty(grid_c1) & grid_res_c1 == 5
           
              if any([~isempty(tropoModel1) ~isnan(tropoModel1)]) & any([strncmpi(tropoModel1,'GPT3 (5° x 5°)',14) strncmpi(tropoModel1,'GPT3 (5° x 5°)',14)])
                 grid_met_rb5_gpt3     = grid_c1;
                 grid_res_met_rb5_gpt3 = 5;
              end  
           
           %GET GRIDFILE FROM HYDROSTATIC TROPO MODELs   
           elseif ~isempty(grid_h1) & grid_res_h1 == 5
                  if any([~isempty(dryModel1) ~isnan(dryModel1)]) & any([strncmpi(dryModel1,'GPT3 (5° x 5°)',14) strncmpi(dryModel1,'GPT3 (5° x 5°)',14)])
                     grid_met_rb5_gpt3     = grid_h1;
                     grid_res_met_rb5_gpt3 = 5;
                  end  
           
           %GET GRIDFILE FROM WET TROPO MODELs       
           elseif ~isempty(grid_w1) & grid_res_w1 == 5
                  if any([~isempty(wetModel1) ~isnan(wetModel1)]) & any([strncmpi(wetModel1,'GPT3 (5° x 5°)',14) strncmpi(wetModel1,'GPT3 (5° x 5°)',14)])
                     grid_met_rb5_gpt3     = grid_w1;
                     grid_res_met_rb5_gpt3 = 5;
                  end  
              
           else %IF ALL FAIL,READ GRID
                grid_met_rb5_gpt3=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE
                grid_res_met_rb5_gpt3 = 5;
           end
           
           %IF ALL OPTIONS FAIL
           if any ([~exist('grid_met_rb5_gpt3','var') ~exist('grid_res_met_rb5_gpt3','var')])
             
              %READ FROM GRID FILE
              grid_met_rb5_gpt3=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE
              grid_res_met_rb5_gpt3 = 5;%GRID RESOLUTION 
            
           elseif any ([exist('grid_met_rb5_gpt3','var') exist('grid_res_met_rb5_gpt3','var')]) & any ([isempty(grid_met_rb5_gpt3) isempty(grid_res_met_rb5_gpt3)])
                
                  %READ FROM GRID FILE
                  grid_met_rb5_gpt3=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE
                  grid_res_met_rb5_gpt3 = 5;%GRID RESOLUTION
           end  
         
        end
        
        %CHECK AGAIN TO BE SURE
        if ~exist('grid_met_rb5_gpt3','var')
           grid_met_rb5_gpt3=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE
        
        elseif all([exist('grid_met_rb5_gpt3','var'),length(grid_met_rb5_gpt3)<14])
               grid_met_rb5_gpt3=readGPTgrid('gpt3_5.mat','GPT3',5);%Combine model GRIDFILE
        end            
        
       %SAVE GRID FILE & RESOLUTION      
       setappdata(0,'grid_met_rb5_gpt3',grid_met_rb5_gpt3); %gridfile
       setappdata(0,'grid_res_met_rb5_gpt3',grid_res_met_rb5_gpt3);%grid resolution 
        
   end %if any([strcmpi(metPARA,'4.Generate from GPT2w model') isequal(metVAL,4)] )
%--------------------------------------------------------------------------       
%******SAVE GRID FILE & RESOLUTION
%IF MODEL IS GPT2w
if any([strcmpi(metPARA,'4.Generate from GPT2w model') isequal(metVAL,4)] )
   setappdata(0,'grid_met_rb5',grid_met_rb5_gpt2w); %gridfile
   setappdata(0,'grid_res_met_rb5',grid_res_met_rb5_gpt2w);%grid resolution
   setappdata(0,'grid_model','GPT2w');%grid model
   
   %*****ASSIGNMENT
   gridFile  = grid_met_rb5_gpt2w;
   gridRES   = grid_res_met_rb5_gpt2w;
   gridModel = 'GPT2w'; 
   
   %IF MODEL IS GPT3
elseif any([strcmpi(metPARA,'5.Generate from GPT3 model') isequal(metVAL,5)] )
       setappdata(0,'grid_met_rb5',grid_met_rb5_gpt3); %gridfile
       setappdata(0,'grid_res_met_rb5',grid_res_met_rb5_gpt3);%grid resolution
       setappdata(0,'grid_model','GPT3');%grid model
       
       %*****ASSIGNMENT
       gridFile  = grid_met_rb5_gpt3;
       gridRES   = grid_res_met_rb5_gpt3;
       gridModel = 'GPT3';
       
else 
     gridFile  = [];
     gridRES   = [];
     gridModel = []; 
       
end

if any([exist('gridFile','var') exist('gridRES','var') exist('gridModel','var')])
    
   %SAVE TO BE CALLED LATER FOR MF POPUP MENUs
   setappdata(0,'gridFile',gridFile); %gridfile
   setappdata(0,'gridRES',gridRES);%grid resolution
   setappdata(0,'gridModel',gridModel);%grid model 
end   

end %if get(hObject,'Value') ==1

%===============================END OF GridResolution_5_rb_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%                          =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================END OF METparameters_popm_Callback.m
%                          =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%               ------------------------------------------------
%4.***************ADD VMD GRIDDED ZENITH DELAYS TO THE OUPT FILE
%               ------------------------------------------------
% --- EXECUTES ON BUTTON PRESS IN cb_extract_VMF_ZTDs
function extract_VMF_ZTDs_cb_Callback(hObject,handles)

if get(hObject,'Value') == 1

   
    %CHECK IF VMF GRID FILEs ARE AVAILABLE IN goGPS
    %Call for the 'SearchVMFgrids.m' fxn
    [VMF_grid_found,VMFgrids] = SearchVMFgrids(handles); 
                        
    %*********SAVE STRUCT GRID FILES
    setappdata(0,'VMFgrids',VMFgrids)
         
    %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
    setappdata(0,'VMF_grid_found',VMF_grid_found)
    
    
    if VMF_grid_found == 0 %IF VMF grid FILES NOT FOUND
           
        beep %Give a beep sound
        errmsg302{1}=sprintf('No VMF grid file(s) found in goGPS directory.\n');
        errmsg302{2}=sprintf('Please Provide VMF grid file(s) & Try Again.\n');
        errmsg302{3}=sprintf('VMF grid(s) can as well be downloaded from : http://vmf.geo.tuwien.ac.at/trop_products/GRID/\n');
        warndlg(errmsg302,'VMF grid file Error','modal')

        %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
        set(handles.cb_extract_VMF_ZTDs,'value',0)
               
    else %IF VMF grid FILES FOUND,LET USER CHOOSE GRID TYPE & RESOLUTION
         %I.E.['VMF1(2°x2.5°)','VMF3(1°x1°)','VMF3(5°x5°)']
         
         [VMF_type,grid_res] = VMFgrid_type(); 
           
           %*********CHECK IF VMF_type,grid_res ARE BOTH EMPTY
           %NOTE:
           %     VMF_type,grid_res ARE ASSIGNED EMPTY WHEN USER CANCEL SELECTION
           %     OR CLOSE DIALOGUE BOX FOR THE SELECTION OF VMF GRID TYPE &
           %     RESOLUTION.IF IT HAPPENS SO, THE SELECTED TROPO DELAY MODEL,
           %     'VMF gridded ZTD' WOULD BE REPLACED WITH SAASTEMOIN MODEL
           %     OR THE POPUP MENU WILL BE SET TO SAASTAMOINEN AS DEFAULT
        
           if all([isempty(VMF_type),isempty(grid_res)])%IF USER CANCEL SELECTION/CLOSE DIALOGUE BOX
            
              %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
              set(handles.cb_extract_VMF_ZTDs,'value',0)
              
           end
           
           %SAVE GRID FILE VERSION & TYPE
           setappdata(0,'VMFgrid_type',VMF_type)
           setappdata(0,'VMFgrid_res',grid_res)
               
    end %//if VMF_grid_found == 0 
       
end  %//if get(hObject,'Value') == 1


%             ----------------------------------------------------
%5.***********RETRIEVE PRECIPITABLE WATER VAPOUR(PWV) CHECKBOX(cb)
%             ----------------------------------------------------
% --- EXECUTES ON BUTTON PRESS IN cb_retrievePWV.
function retrievePWV_cb_Callback(hObject,handles)
 
%GET Tm popup menu contents
getconts = cellstr(get(handles.pop_Tm,'String'));%get popup menu contents as cell array
TmMODEL  = getconts{get(handles.pop_Tm,'Value')};%get selected item from pop_Tm
TmVAL    = get(handles.pop_Tm,'Value');%GET SELECTED VALUE    
    
if get(hObject,'Value') == 1
    
   %CHANGE TOOL TIP STRING
   set(hObject,'TooltipString','Selected to retrieve Precipitable Water Vapur(PWV)')
 
   %SET WEIGHTED MEAN TEMPERATURE(Tm) MODELs POPUP MENU VISIBLE ON
   set(handles.pop_Tm,'Enable','On','visible','On')
   set(handles.text_Tm,'Enable','On','visible','On')
   set(handles.text_grid_resolution_pwv,'visible',get(handles.text_grid_resolution_pwv,'visible'))
   set(handles.text_grid_resolution_pwv,'visible',get(handles.rb_grid_resolution_1_pwv,'visible'),...
                                        'value',get(handles.rb_grid_resolution_1_pwv,'value'))
   set(handles.text_grid_resolution_pwv,'visible',get(handles.rb_grid_resolution_5_pwv,'visible'),...
                                        'value',get(handles.rb_grid_resolution_5_pwv,'value'))
 
                                    
   if any([strncmpi(TmMODEL,'GPT2w model [Böhm et al 2014]',28),strncmpi(TmMODEL,'GPT3 model [Landskron & Böhm, 2018]',35),...
            any(TmVAL==[11,12])])
      %GRID RESOLUTION RADIO BUTTON & TEXT ENABLE ON
      set(handles.rb_grid_resolution_1_pwv,'enable','on')
      set(handles.rb_grid_resolution_5_pwv,'enable','on')
      set(handles.text_grid_resolution_pwv,'enable','on')
        
    end                                 
                                    
else 
    %CHANGE TOOL TIP STRING
    set(hObject,'TooltipString','Select to retrieve Precipitable Water Vapur(PWV)')
 
    %SET WEIGHTED MEAN TEMPERATURE(Tm) MODELs POPUP MENU VISIBLE ON
    set(handles.pop_Tm,'Enable','Off')
    set(handles.text_Tm,'Enable','Off')
    
    if any([strncmpi(TmMODEL,'GPT2w model [Böhm et al 2014]',28),strncmpi(TmMODEL,'GPT3 model [Landskron & Böhm, 2018]',35),...
            any(TmVAL==[11,12])])
        %GRID RESOLUTION RADIO BUTTON & TEXT ENABLE OFF
        set(handles.rb_grid_resolution_1_pwv,'enable','off')
        set(handles.rb_grid_resolution_5_pwv,'enable','off')
        set(handles.text_grid_resolution_pwv,'enable','off')
        
    end
    
end
%===============================END OF retrievePWV_cb_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+


%             ----------------------------------------------------
%6.***********WEIGHTED MEAN TEMPERATURE(Tm) MODELS POPUP MENU
%             ----------------------------------------------------

% ---EXECUTES ON SELECTION CHANGE IN pop_Tm.
function Tm_model_pop_Callback(hObject,handles)

%GET Tm popup menu contents
getconts = cellstr(get(handles.pop_Tm,'String'));%get popup menu contents as cell array
TmMODEL  = getconts{get(handles.pop_Tm,'Value')};%get selected item from pop_Tm
TmVAL    = get(handles.pop_Tm,'Value');%GET SELECTED VALUE

%GET TOGGLE GRID RESOLUTION TYPE
grid_res1_pwv = get(handles.rb_grid_resolution_1_pwv,'value');%for 1° resolution
grid_res5_pwv = get(handles.rb_grid_resolution_5_pwv,'value');%for 5° resolution


%HID TOGGLE STATE OF GRID RESOLUTION TYPE (1 OR 5)
if (any(TmVAL==[1:10,13,14,15]))
    
   %GRID RESOLUTION RADIO BUTTON & TEXT VISIBLE OFF
   set(handles.rb_grid_resolution_1_pwv,'visible','off')
   set(handles.rb_grid_resolution_5_pwv,'visible','off')
   set(handles.text_grid_resolution_pwv,'visible','off')
   
else %IF SELECTED IS EITHER GPT2w or GPT3 
      
    %GET STORED GRID VALUES & RESOLUTION FROM Tm MODELS(I.E. GPT MODELS)
    grid_Tm1     = getappdata(0,'grid_Tm');%Previous gridvalues
    grid_res_Tm1 = getappdata(0,'grid_res_Tm');%Previous grid resolution
    Tm_model1    = getappdata(0,'Tm_model');% Tm MODEL IF MODEL IS GPT2w or GPT3

    %GET STORED GRID VALUES & RESOLUTION FROM GPT MODELS 1° GRIDRESOLUTION RADIO BUTTON
    grid_Tm_1     = getappdata(0,'grid_Tm_1');%Previous gridvalues from 1° rb
    grid_res_Tm_1 = getappdata(0,'grid_res_Tm_1');%Previous grid resolution 1° rb
    TmModel_1     = getappdata(0,'TmModel_1');% Tm MODEL IF MODEL IS GPT2w or GPT3
      
    %GET STORED GRID VALUES & RESOLUTION FROM GPT MODELS 5° GRIDRESOLUTION RADIO BUTTON
    grid_Tm_5     = getappdata(0,'grid_Tm_5');%Previous gridvalues from 5° rb
    grid_res_Tm_5 = getappdata(0,'grid_res_Tm_5');%Previous grid resolution 5° rb
    TmModel_5     = getappdata(0,'TmModel_5');% Tm MODEL IF MODEL IS GPT2w or GPT3
    
    %*******MAKE TOGGLE STATE OF GRID RESOLUTION TYPE (1 OR 5) VISIBLE
    %GRID RESOLUTION RADIO BUTTON & TEXT VISIBLE ON
    set(handles.rb_grid_resolution_1_pwv,'visible','on','enable','on')
     
    if get(handles.rb_grid_resolution_1_pwv,'value')==0 & get(handles.rb_grid_resolution_5_pwv,'value')==0
       set(handles.rb_grid_resolution_5_pwv,'Visible','on','Value',1,'enable','on','ForegroundColor',[0,0.75,0.75])%set as default
    else  
        set(handles.rb_grid_resolution_5_pwv,'visible','on','enable','on')
    end 
         
    set(handles.text_grid_resolution_pwv,'visible','on','enable','on')
     
    %TIME VARIATION INDICATION
    TIMEvariation() %Time variation sub-routine
    uiwait(gcf)
    setappdata(0,'Timevar_Tm',getappdata(0,'Timevar'))%save for use
     
    
    switch getappdata(0,'Timevar')
        case 0
               Timevar_gpt = 'time variation (annual and semiannual terms)';
                         
        case  1
               Timevar_gpt = 'no time variation but static quantities';
    end          
    
end  

%*********************SET TOOLTIP STRINGS 
if any([strncmpi(TmMODEL,'Bevis et al (1992) model',24),TmVAL == 1])
   %SET TOOLTIP STRING
   set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','Bevis et al (1992,1994) Tm model','i.e. Tm = 70.2 + 0.72.*T '}))%change Tooltip any time user select Tm Model type 
   
   %CREATE OUTPUT STRING FOR Tm MODEL TYPE
   TmModel_type = 'Bevis et al (1992) [i.e. Tm = 70.2 + 0.72.*T]' ;%Needed for goGPS OUTPUT file
   
elseif any([strncmpi(TmMODEL,'Bevis et al (1995) model',24),TmVAL == 2])
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','Bevis et al (1995) Tm model','i.e. Tm = 86.63 + 0.668.*T '}))%change Tooltip any time user select Tm Model type 
       
       %CREATE OUTPUT STRING FOR Tm MODEL TYPE
       TmModel_type = 'Bevis et al (1995) [i.e. Tm = 86.63 + 0.668.*T]' ;%Needed for goGPS OUTPUT file
       
elseif any([strncmpi(TmMODEL,'Mendes et al (2000) model',25),TmVAL == 3])
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','Mendes et al (2000) Tm model','i.e. Tm = 50.4 + 0.789.*T '}))%change Tooltip any time user select Tm Model type 
       
       %CREATE OUTPUT STRING FOR Tm MODEL TYPE
       TmModel_type = 'Mendes et al (2000) [i.e. Tm = 50.4 + 0.789.*T]' ;%Needed for goGPS OUTPUT file
       
elseif any([strncmpi(TmMODEL,'Schueler et al (2001) model',27),TmVAL == 4])
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','Schueler et al (2001) Tm model','i.e. Tm = 86.9 + 0.647.*T '}))%change Tooltip any time user select Tm Model type 
       
       %CREATE OUTPUT STRING FOR Tm MODEL TYPE
       TmModel_type = 'Schueler et al (2001) [i.e. Tm = 86.9 + 0.647.*T]' ;%Needed for goGPS OUTPUT file
       
elseif any([strncmpi(TmMODEL,'Yao et al (2014) model',22),TmVAL == 5])
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','Yao et al (2014) Tm model','i.e. Tm = 43.69 + 0.8116.*T '}))%change Tooltip any time user select Tm Model type 
       
       %CREATE OUTPUT STRING FOR Tm MODEL TYPE
       TmModel_type = 'Yao et al (2014) [i.e. Tm = 43.69 + 0.8116.*T]' ;%Needed for goGPS OUTPUT file
       
elseif any([strncmpi(TmMODEL,'GTm-I model [Yao et al 2012]',28),TmVAL == 6])
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','Global Weighted Mean Temperature (GWMT-I/GTm-I) model'}))%change Tooltip any time user select Tm Model type 
       
       %CREATE OUTPUT STRING FOR Tm MODEL TYPE
       TmModel_type = 'Global Weighted Mean Temperature I (GWMT-I/GTm-I) [Yao et al 2012]';%Needed for goGPS OUTPUT file
       
elseif any([strncmpi(TmMODEL,'GTm-II model [Yao et al 2013]',29),TmVAL == 7])
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','Global Weighted Mean Temperature (GWMT-II/GTm-II) model'}))%change Tooltip any time user select Tm Model type 
       
       %CREATE OUTPUT STRING FOR Tm MODEL TYPE
       TmModel_type = 'Global Weighted Mean Temperature II (GWMT-II/GTm-II) [Yao et al 2013]';%Needed for goGPS OUTPUT file
       
elseif any([strncmpi(TmMODEL,'GTm-III model [Yao et al 2014]',30),TmVAL == 8])
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','Global Weighted Mean Temperature (GWMT-III/GTm-III) model'}))%change Tooltip any time user select Tm Model type 
       
       %CREATE OUTPUT STRING FOR Tm MODEL TYPE
       TmModel_type = 'Global Weighted Mean Temperature III (GWMT-III/GTm-III) [Yao et al 2014]';%Needed for goGPS OUTPUT file
       
elseif any([strncmpi(TmMODEL,'GTVR-Tm model [Yao et al 2018]',30),TmVAL == 9])
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','Global Weighted Mean Temperature Vertical Correction (GTVR-Tm) model'}))%change Tooltip any time user select Tm Model type 
       
       %GET & SAVE GTVR-Tm model GRID VALUES
        if ~isempty(getappdata(0,'gridV_GTVR'))
            
           setappdata(0,'gridV_GTVR',getappdata(0,'gridV_GTVR'))
           
        else   
            %SEARCH FOR GTVR-Tm model GRID FILE(GTVR-Tm.txt) & DIRECTORY & IMPORT VALUES               
            %Call the "searchGTVRgrid.m" function
            [~,gridV_GTVR] = SearchGTVRgrid();%****IMPORTED TVGG-Tm MODEL COEFFICIENTS
       
            %SAVE GRID VALUES
            setappdata(0,'gridV_GTVR',gridV_GTVR)
        end
       
       %GET & SAVE GPT2w model GRID VALUES
       if ~isempty(getappdata(0,'GPT2w_gridV'))
           
          setappdata(0,'GPT2w_gridV',getappdata(0,'GPT2w_gridV'))
          
       else    
           %READ GPT2w GRID FILE
           GPT2w_gridV     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
       
          %SAVE GPT2w GRID VALUES
          setappdata(0,'GPT2w_gridV',GPT2w_gridV)
          
       end
       
       %CREATE OUTPUT STRING FOR Tm MODEL TYPE
       TmModel_type = 'Global Weighted Mean Temperature Vertical Correction (GTVR-Tm) model [Yao et al 2018]';%Needed for goGPS OUTPUT file

elseif any([strncmpi(TmMODEL,'TVGG-Tm model [Jiang et al 2018]',32),TmVAL == 10])
    
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','Time-Varying Global Gridded (TVGG) Tm model'}))%change Tooltip any time user select Tm Model type 
       
       %GET & SAVE TVGG-Tm model GRID VALUES
       if ~isempty(getappdata(0,'gridV_TVGG'))
       
          setappdata(0,'gridV_TVGG',getappdata(0,'gridV_TVGG'))
          
       else          
           %SEARCH FOR TVGG-Tm model GRID FILE(Coeff_TVGG_ERAI.mat) & DIRECTORY & LOAD VALUES               
           %Call the "SearchTVGGgrid.m" fxn
           [~,gridV_TVGG] = SearchTVGGgrid();%****LOADED TVGG-Tm MODEL COEFFICIENTS
       
           %SAVE GRID VALUES
           setappdata(0,'gridV_TVGG',gridV_TVGG)
           
       end
          
       %CREATE OUTPUT STRING FOR Tm MODEL TYPE
       TmModel_type = 'Time-Varying Global Gridded (TVGG) model [Jiang et al 2018]';%Needed for goGPS OUTPUT file
             
elseif any([strncmpi(TmMODEL,'GPT2w model [Böhm et al 2014]',28),TmVAL == 11])
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','Global Pressure & Temperature to wet (GPT2w) model'}))%change Tooltip any time user select Tm Model type 
        
        
       %*************EXTRACT GRID VALUEs 
       if grid_res1_pwv == 1 & grid_res5_pwv == 0
          
          %CREATE OUTPUT STRING FOR Tm MODEL TYPE
          TmModel_type =  strjoin({'GPT2w based on a 1° x 1° grid file with',Timevar_gpt,'[Böhm et al 2014]'});%Needed for goGPS OUTPUT file
           
          if ~isempty(grid_Tm1) & grid_res_Tm1 == 1 & ~isempty(Tm_model1)
          
             if strcmpi(Tm_model1,'GPT2w')% IF Tm MODEL IS GPT2w
                grid_Tm     = grid_Tm1; %Gridfile for 1° resolution radio button
                grid_res_Tm = 1;  %Grid resolution for 1° resolution radio button 
            
             else  
                 grid_Tm     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                 grid_res_Tm = 1;%GRID RESOLUTION
                 
             end
             
          elseif ~isempty(grid_Tm_1) & grid_res_Tm_1 == 1 & ~isempty(TmModel_1)
             
                 if strcmpi(TmModel_1,'GPT2w')% IF Tm MODEL IS GPT2w
                    grid_Tm     = grid_Tm_1; %Gridvalues for 1° resolution radio button
                    grid_res_Tm = 1;  %Grid resolution for 1° resolution radio button
                
                 else      
                     grid_Tm     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                     grid_res_Tm = 1;%GRID RESOLUTION
                 
                 end 
         
          else 
              grid_Tm     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
              grid_res_Tm = 1;%GRID RESOLUTION
              
          end %if ~isempty(grid_Tm1) & grid_res_Tm1 == 1 & ~isempty(Tm_model1)
             
                      
       elseif grid_res1_pwv == 0 & grid_res5_pwv == 1
              
              %CREATE OUTPUT STRING FOR Tm MODEL TYPE
              TmModel_type =  strjoin({'GPT2w based on a 5° x 5° grid file with',Timevar_gpt,'[Böhm et al 2014]'});%Needed for goGPS OUTPUT file
           
              if ~isempty(grid_Tm1) & grid_res_Tm1 == 5 & ~isempty(Tm_model1)
          
                 if strcmpi(Tm_model1,'GPT2w')% IF Tm MODEL IS GPT2w
                    grid_Tm     = grid_Tm1; %Gridfile for 5° resolution radio button
                    grid_res_Tm = 5;  %Grid resolution for 5° resolution radio button 
            
                 else   
                     grid_Tm     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                     grid_res_Tm = 5;%GRID RESOLUTION
                 
                 end 
             
          elseif ~isempty(grid_Tm_5) & grid_res_Tm_5 == 5 & ~isempty(TmModel_5)
             
                 if strcmpi(TmModel_5,'GPT2w')% IF Tm MODEL IS GPT2w
                    grid_Tm     = grid_Tm_5; %Gridvalues for 5° resolution radio button
                    grid_res_Tm = 5;  %Grid resolution for 5° resolution radio button
                
                 else      
                     grid_Tm     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                     grid_res_Tm = 5;%GRID RESOLUTION
                 
                 end 
         
              else  
                  grid_Tm     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                  grid_res_Tm = 5;%GRID RESOLUTION
              
              end  %if ~isempty(grid_Tm1) & grid_res_Tm1 == 1 & ~isempty(Tm_model1)
                      
       end
       
       if all([~exist('grid_Tm','var'),exist('grid_res_Tm','var')])
          
          if grid_res_Tm == 1
              
             grid_Tm = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
             
             if ~exist(' TmModel_type','var')
                 
                %CREATE OUTPUT STRING FOR Tm MODEL TYPE
                TmModel_type =  strjoin({'GPT2w based on a 1° x 1° grid file with',Timevar_gpt,'[Böhm et al 2014]'});%Needed for goGPS OUTPUT file
           
             end
             
          elseif grid_res_Tm == 5 
              
                 grid_Tm = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                 
                 if ~exist(' TmModel_type','var')
                 
                    %CREATE OUTPUT STRING FOR Tm MODEL TYPE
                    TmModel_type =  strjoin({'GPT2w based on a 5° x 5° grid file with',Timevar_gpt,'[Böhm et al 2014]'});%Needed for goGPS OUTPUT file
           
                 end 
                 
          end
          
       elseif all([~exist('grid_Tm','var'),~exist('grid_res_Tm','var')])
             
              %SET DEFAULT GRAD VALUES
              grid_Tm = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
               
              grid_res_Tm = 1;%GRID RESOLUTION
              
              if ~exist(' TmModel_type','var')
                 
                %CREATE OUTPUT STRING FOR Tm MODEL TYPE
                TmModel_type =  strjoin({'GPT2w based on a 1° x 1° grid file with',Timevar_gpt,'[Böhm et al 2014]'});%Needed for goGPS OUTPUT file
           
              end 
              
       end 
              
       %SAVE GRID FILE & RESOLUTION      
       setappdata(0,'grid_Tm',grid_Tm); %gridvalues
       setappdata(0,'grid_res_Tm',grid_res_Tm);%grid resolution
       setappdata(0,'Tm_model','GPT2w')%SAVE Tm MODEL IF MODEL IS GPT2w
       
elseif any([strncmpi(TmMODEL,'GPT3 model [Landskron & Böhm, 2018]',35),TmVAL == 12])
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','Global Pressure & Temperature 3 (GPT3) model'}))%change Tooltip any time user select Tm Model type 
        
       %EXTRACT GRID VALUEs
       if grid_res1_pwv == 1 & grid_res5_pwv == 0
          
          %CREATE OUTPUT STRING FOR Tm MODEL TYPE
          TmModel_type =  strjoin({'GPT3 based on a 1° x 1° grid file with',Timevar_gpt,'[Landskron & Böhm, 2018]'});%Needed for goGPS OUTPUT file
          
          if ~isempty(grid_Tm1) & grid_res_Tm1 == 1 & ~isempty(Tm_model1)
          
             if strcmpi(Tm_model1,'GPT3')% IF Tm MODEL IS GPT3
                grid_Tm     = grid_Tm1; %Gridfile for 1° resolution radio button
                grid_res_Tm = 1;  %Grid resolution for 1° resolution radio button 
            
             else  
                 grid_Tm     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                 grid_res_Tm = 1;%GRID RESOLUTION
                 
             end
             
          elseif ~isempty(grid_Tm_1) & grid_res_Tm_1 == 1 & ~isempty(TmModel_1)
             
                 if strcmpi(TmModel_1,'GPT3')% IF Tm MODEL IS GPT3
                    grid_Tm     = grid_Tm_1; %Gridvalues for 1° resolution radio button
                    grid_res_Tm = 1;  %Grid resolution for 1° resolution radio button
                
                 else      
                     grid_Tm     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                     grid_res_Tm = 1;%GRID RESOLUTION
                 
                 end 
         
          else  
              grid_Tm     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
              grid_res_Tm = 1;%GRID RESOLUTION
              
          end  %if ~isempty(grid_Tm1) & grid_res_Tm1 == 1 & ~isempty(Tm_model1)
                      
       elseif grid_res1_pwv == 0 & grid_res5_pwv == 1
              
              %CREATE OUTPUT STRING FOR Tm MODEL TYPE
              TmModel_type =  strjoin({'GPT3 based on a 5° x 5° grid file with',Timevar_gpt,'[Landskron & Böhm, 2018]'});%Needed for goGPS OUTPUT file
           
              
              if ~isempty(grid_Tm1) & grid_res_Tm1 == 5 & ~isempty(Tm_model1)
          
                 if strcmpi(Tm_model1,'GPT3')% IF Tm MODEL IS GPT3
                    grid_Tm     = grid_Tm1; %Gridfile for 5° resolution radio button
                    grid_res_Tm = 5;  %Grid resolution for 5° resolution radio button 
            
                 else   
                     grid_Tm     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                     grid_res_Tm = 5;%GRID RESOLUTION
                 
                 end 
             
          elseif ~isempty(grid_Tm_5) & grid_res_Tm_5 == 5 & ~isempty(TmModel_5)
             
                 if strcmpi(TmModel_5,'GPT3')% IF Tm MODEL IS GPT3
                    grid_Tm     = grid_Tm_5; %Gridvalues for 5° resolution radio button
                    grid_res_Tm = 5;  %Grid resolution for 5° resolution radio button
                
                 else      
                     grid_Tm     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                     grid_res_Tm = 5;%GRID RESOLUTION
                 
                 end 
         
              else   
                  grid_Tm     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                  grid_res_Tm = 5;%GRID RESOLUTION
              
              end %if ~isempty(grid_Tm1) & grid_res_Tm1 == 1 & ~isempty(Tm_model1)
              
       end %//if grid_res1_pwv == 1 & grid_res5_pwv == 0
       
       if all([~exist('grid_Tm','var'),exist('grid_res_Tm','var')])
          
          if grid_res_Tm == 1
              
             grid_Tm = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
               
             if ~exist(' TmModel_type','var')
                 
                %CREATE OUTPUT STRING FOR Tm MODEL TYPE
                TmModel_type =  strjoin({'GPT3 based on a 1° x 1° grid file with',Timevar_gpt,'[Landskron & Böhm, 2018]'});%Needed for goGPS OUTPUT file
           
             end  
             
          elseif grid_res_Tm == 5 
              
                 grid_Tm = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                 
                 if ~exist(' TmModel_type','var')
                 
                    %CREATE OUTPUT STRING FOR Tm MODEL TYPE
                    TmModel_type =  strjoin({'GPT3 based on a 5° x 5° grid file with',Timevar_gpt,'[Landskron & Böhm, 2018]'});%Needed for goGPS OUTPUT file
           
                 end   
                 
          end
          
       elseif all([~exist('grid_Tm','var'),~exist('grid_res_Tm','var')])
             
              grid_Tm = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
               
              grid_res_Tm = 1;%GRID RESOLUTION
              
              if ~exist(' TmModel_type','var')
                 
                 %CREATE OUTPUT STRING FOR Tm MODEL TYPE
                 TmModel_type =  strjoin({'GPT3 based on a 1° x 1° grid file with',Timevar_gpt,'[Landskron & Böhm, 2018]'});%Needed for goGPS OUTPUT file
           
              end   
              
       end
       
       %SAVE GRID FILE & RESOLUTION      
       setappdata(0,'grid_Tm',grid_Tm); %gridvalues
       setappdata(0,'grid_res_Tm',grid_res_Tm);%grid resolution
       setappdata(0,'Tm_model','GPT3')%SAVE Tm MODEL IF MODEL IS GPT3
       
elseif any([strncmpi(TmMODEL,'UNB3m model [Leandro et al 2006]',32),TmVAL == 13])
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','UNB3m [Leandro et al 2006] model'}))%change Tooltip any time user select Tm Model type 
       
       %CREATE OUTPUT STRING FOR Tm MODEL TYPE
       TmModel_type = 'University of New Brunswick v3 modified(UNB3m) [Leandro et al 2006]';%Needed for goGPS OUTPUT file
       
elseif any([strncmpi(TmMODEL,'Simplified PI  model [Manandhar et al 2017]',43),TmVAL == 14])
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','Simplified PI  model [Manandhar et al 2017]'}))%change Tooltip any time user select Tm Model type 
       
       %CREATE OUTPUT STRING FOR Tm MODEL TYPE
       TmModel_type = 'Simplified PI  model [Manandhar et al 2017]';%Needed for goGPS OUTPUT file
       
elseif  any([strncmpi(TmMODEL,'GTrop model [Sun et al 2019]',28),TmVAL == 15])
    
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','Global Tropospheric model(GTrop)'}))%change Tooltip any time user select Tm Model type 
       
       %GET & SAVE GTrop GRID VALUES
        if ~isempty(getappdata(0,'gridV_GTrop'))
            
           setappdata(0,'gridV_GTrop',getappdata(0,'gridV_GTrop')) 
           
        else
            %SEARCH FOR GTrop model GRID FILE(GTropCoefficient.mat) & DIRECTORY & LOAD VALUES               
            %Call the "SearchGTropgrid.m" fxn
            [~,gridV_GTrop] = SearchGTropgrid();%****LOADED GTrop MODEL COEFFICIENTS
       
            %SAVE GRID VALUES
            setappdata(0,'gridV_GTrop',gridV_GTrop)
        end
          
       %CREATE OUTPUT STRING FOR Tm MODEL TYPE
       TmModel_type = 'Global Tropospheric model(GTrop) [Sun et al 2019]';%Needed for goGPS OUTPUT file
end 

%SAVE TmModel_type
setappdata(0,'TmModel_type',TmModel_type)

%===============================END OF Tm_model_pop_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%                   ------------------
%(6.1)*********** 1° GRID RESOLUTION 
%                   ------------------
% --- EXECUTES ON BUTTON PRESS IN rb_grid_resolution_1.
function GridResolution_1_rb_pwv_Callback(hObject,handles)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%*******FOR READING GLOBAL PRESSURE & TEMPERATURE(GPT) GRIDFILES
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%GET Tm popup menu contents
getconts = cellstr(get(handles.pop_Tm,'String'));%get popup menu contents as cell array
TmMODEL  = getconts{get(handles.pop_Tm,'Value')};%get selected item from pop_Tm
TmVAL    = get(handles.pop_Tm,'Value');%GET SELECTED VALUE

%GET STORED GRIDVALUES & RESOLUTION FROM Tm MODELS(I.E. GPT MODELS)
grid_Tm1     = getappdata(0,'grid_Tm');%Previous gridvalues
grid_res_Tm1 = getappdata(0,'grid_res_Tm');%Previous grid resolution
Tm_model1    = getappdata(0,'Tm_model');% Tm MODEL IF MODEL IS GPT2w or GPT3

%GET STORED GRIDVALUES & RESOLUTION FROM GPT MODELS 1° GRIDRESOLUTION RADIO BUTTON
grid_Tm_1     = getappdata(0,'grid_Tm_1');%Previous gridvalues from 1° rb
grid_res_Tm_1 = getappdata(0,'grid_res_Tm_1');%Previous grid resolution 1° rb
TmModel_1     = getappdata(0,'TmModel_1');% Tm MODEL IF MODEL IS GPT2w or GPT3

%SET 1° GRID RESOLUTION RADIO BUTTON TO 1
set(hObject,'Value',1,'ForegroundColor',[0,0.75,0.75])

%SET 5° GRID RESOLUTION RADIO BUTTON TO 0
set(handles.rb_grid_resolution_5_pwv,'Value',0,'ForegroundColor','black')

if get(hObject,'Value') ==1 %IF USER SELECTS 1° resolution radio button
    
   %CHECK IF SELECTED Tm MODEL IS GPT2w MODEL  
   if any([strncmpi(TmMODEL,'GPT2w model [Böhm et al 2014]',28),TmVAL == 11])
      
      TmModel = 'GPT2w';%INDICATE Tm model type
       
      if ~isempty(grid_Tm1) & grid_res_Tm1 == 1 & ~isempty(Tm_model1)
          
         if strcmpi(Tm_model1,'GPT2w')% IF Tm MODEL IS GPT2w
            grid_Tm     = grid_Tm1; %Gridfile for 1° resolution radio button
            grid_res_Tm = 1;  %Grid resolution for 1° resolution radio button 
            
         else
              grid_Tm     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
              grid_res_Tm = 1;%GRID RESOLUTION
              
         end
         
      elseif ~isempty(grid_Tm_1) & grid_res_Tm_1 == 1 & ~isempty(TmModel_1)
             
             if strcmpi(TmModel_1,'GPT2w')% IF Tm MODEL IS GPT2w
                grid_Tm     = grid_Tm_1; %Gridfile for 1° resolution radio button
                grid_res_Tm = 1;  %Grid resolution for 1° resolution radio button
                
             else    
                 grid_Tm     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                 grid_res_Tm = 1;%GRID RESOLUTION
                 
             end
         
      else
          grid_Tm     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
          grid_res_Tm = 1;%GRID RESOLUTION
      end
      
   elseif any([strncmpi(TmMODEL,'GPT3 model [Landskron & Böhm, 2018]',35),TmVAL == 12])
    
          TmModel = 'GPT3';%INDICATE Tm model type
          
          if ~isempty(grid_Tm1) & grid_res_Tm1 == 1 & ~isempty(Tm_model1)
          
             if strcmpi(Tm_model1,'GPT3')% IF Tm MODEL IS GPT2w
                grid_Tm     = grid_Tm1; %Gridfile for 1° resolution radio button
                grid_res_Tm = 1;  %Grid resolution for 1° resolution radio button 
            
             else  
                 grid_Tm     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                 grid_res_Tm = 1;%GRID RESOLUTION
              
             end
             
          elseif ~isempty(grid_Tm_1) & grid_res_Tm_1 == 1 & ~isempty(TmModel_1)
              
                 if strcmpi(TmModel_1,'GPT3')% IF Tm MODEL IS GPT3
                    grid_Tm     = grid_Tm_1; %Gridfile for 1° resolution radio button
                    grid_res_Tm = 1;  %Grid resolution for 1° resolution radio button
                    
                 else  
                     grid_Tm     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                     grid_res_Tm = 1;%GRID RESOLUTION
                     
                 end 
         
          else 
              grid_Tm     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
              grid_res_Tm = 1;%GRID RESOLUTION
          
          end
          
   end %//if any([strncmpi(TmMODEL,'GPT2w model [Böhm et al 2014]',28),TmVAL == 11])
   
   %SAVE GRID FILE & RESOLUTION
   if any([exist('grid_Tm','var') exist('grid_res_Tm','var')])
      setappdata(0,'grid_Tm_1',grid_Tm); %gridvalues
      setappdata(0,'grid_res_Tm_1',grid_res_Tm);%grid resolution   
      setappdata(0,'TmModel_1',TmModel);%GPT MODEL TYPE FOR Tm    
   end
         
end %//if get(hObject,'Value') ==1

%===============================END OF GridResolution_1_rb_pwv_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%                   ------------------
%(6.2)************ 5° GRID RESOLUTION 
%                   ------------------
% --- EXECUTES ON BUTTON PRESS IN rb_grid_resolution_1.
function GridResolution_5_rb_pwv_Callback(hObject,handles)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%*******FOR READING GLOBAL PRESSURE & TEMPERATURE(GPT) GRIDFILES
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%GET Tm popup menu contents
getconts = cellstr(get(handles.pop_Tm,'String'));%get popup menu contents as cell array
TmMODEL  = getconts{get(handles.pop_Tm,'Value')};%get selected item from pop_Tm
TmVAL    = get(handles.pop_Tm,'Value');%GET SELECTED VALUE

%GET STORED GRIDVALUES & RESOLUTION FROM Tm MODELS popup menu(I.E. GPT MODELS)
grid_Tm1     = getappdata(0,'grid_Tm');%Previous gridvalues
grid_res_Tm1 = getappdata(0,'grid_res_Tm');%Previous grid resolution
Tm_model1    = getappdata(0,'Tm_model');% Tm MODEL IF MODEL IS GPT2w or GPT3

%GET STORED GRIDVALUES & RESOLUTION FROM GPT MODELS 5° GRIDRESOLUTION RADIO BUTTON
grid_Tm_5     = getappdata(0,'grid_Tm_5');%Previous gridvalues
grid_res_Tm_5 = getappdata(0,'grid_res_Tm_5');%Previous grid resolution
TmModel_5     = getappdata(0,'TmModel_5');% Tm MODEL IF MODEL IS GPT2w or GPT3

%SET 5° GRID RESOLUTION RADIO BUTTON TO 1
set(hObject,'Value',1,'ForegroundColor',[0,0.75,0.75]) 

%SET 1° GRID RESOLUTION RADIO BUTTON TO 0
set(handles.rb_grid_resolution_1_pwv,'Value',0,'ForegroundColor','black')

if get(hObject,'Value') ==1 %IF USER SELECTS 5° resolution radio button
    
   %CHECK IF SELECTED Tm MODEL IS GPT2w MODEL  
   if any([strncmpi(TmMODEL,'GPT2w model [Böhm et al 2014]',28),TmVAL == 11])
       
      TmModel = 'GPT2w';%INDICATE Tm model type
      
      if ~isempty(grid_Tm1) & grid_res_Tm1 == 5 & ~isempty(Tm_model1)
          
         if strcmpi(Tm_model1,'GPT2w')% IF Tm MODEL IS GPT2w
            grid_Tm     = grid_Tm1; %Gridfile for 5° resolution radio button
            grid_res_Tm = 5;  %Grid resolution for 5° resolution radio button 
  
         else
              grid_Tm     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
              grid_res_Tm = 5;%GRID RESOLUTION
              
         end
         
      elseif ~isempty(grid_Tm_5) & grid_res_Tm_5 == 5 & ~isempty(TmModel_5)
             
             if strcmpi(TmModel_5,'GPT2w')% IF Tm MODEL IS GPT2w
                grid_Tm     = grid_Tm_5; %Gridfile for 5° resolution radio button
                grid_res_Tm = 5;  %Grid resolution for 5° resolution radio button
                
             else    
                 grid_Tm     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                 grid_res_Tm = 5;%GRID RESOLUTION
                 
             end
      else    
          grid_Tm     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
          grid_res_Tm = 5;%GRID RESOLUTION
      end
      
   elseif any([strncmpi(TmMODEL,'GPT3 model [Landskron & Böhm, 2018]',35),TmVAL == 12])
    
          TmModel = 'GPT3';%INDICATE Tm model type
          
          if ~isempty(grid_Tm1) & grid_res_Tm1 == 5 & ~isempty(Tm_model1)
          
             if strcmpi(Tm_model1,'GPT3')% IF Tm MODEL IS GPT3
                grid_Tm     = grid_Tm1; %Gridfile for 5° resolution radio button
                grid_res_Tm = 5;  %Grid resolution for 5° resolution radio button 
                         
             else  
                 grid_Tm     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                 grid_res_Tm = 5;%GRID RESOLUTION
              
             end 
         
          elseif ~isempty(grid_Tm_5) & grid_res_Tm_5 == 5 & ~isempty(TmModel_5)
              
                 if strcmpi(TmModel_5,'GPT3')% IF Tm MODEL IS GPT3
                    grid_Tm     = grid_Tm_5; %Gridfile for 5° resolution radio button
                    grid_res_Tm = 5;  %Grid resolution for 5° resolution radio button
                    
                 else  
                     grid_Tm     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                     grid_res_Tm = 5;%GRID RESOLUTION
                     
                 end
         
          else  
              grid_Tm     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
              grid_res_Tm = 5;%GRID RESOLUTION
          
          end
          
   end %//if any([strncmpi(TmMODEL,'GPT2w model [Böhm et al 2014]',28),TmVAL == 11])
   
   
   %SAVE GRID FILE & RESOLUTION
   if any([exist('grid_Tm','var') exist('grid_res_Tm','var')])
      setappdata(0,'grid_Tm_5',grid_Tm); %gridvalues
      setappdata(0,'grid_res_Tm_5',grid_res_Tm);%grid resolution   
      setappdata(0,'TmModel_5',TmModel);%GPT MODEL TYPE FOR Tm 
   end
         
end %//if get(hObject,'Value') ==1
%===============================END OF GridResolution_5_rb_pwv_Callback.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%                    ------------------------               
%7.******************THE go ! PUSH BUTTON(pb)
%                    ------------------------

% --- EXECUTES ON BUTTON PRESS IN go_button.
function go_pb_Callback(hObject,handles) %#ok<*INUSL>
%--------------------------------------------------------------------------    
%*******************GET ALL ATMOSPHERIC TOOL COMPONENTs 
%                   =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%--------------------------------------------------------------------------+

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%1.**********************IONOSPHERIC MODELLING
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
getIONOconts = cellstr(get(handles.pop_lIono_correction,'String'));%get IIono popup menu contents as cell array
IonoModel=getIONOconts{get(handles.pop_lIono_correction,'Value')};%get selected item from pop_lIono_correction

%SAVE SELECTED IONO MODEL
setappdata(0,'Iono_model',IonoModel)
%============================================IONOSPHERIC MODELLING         +
%--------------------------------------------------------------------------+


%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%2.***************TROPOSPHERIC MODELLING
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%GET SELECTED MODEL TYPE
tropmodel_c=get(handles.rb_Tropo_combine,'Value');%Get combine tropmodel rb
tropmodel_s=get(handles.rb_Tropo_separate,'Value');%Get separate tropmodel rb

%SAVE OPTIONs FOR COMBINE(DRY+WET) & SEPARATE(DRY,WET) TROPO MODELLING
setappdata(0,'option_Cmodels',get(handles.rb_Tropo_combine,'Value'))% option button for combine tropo models
setappdata(0,'option_Smodels',get(handles.rb_Tropo_separate,'Value'))%option button for separate tropo models

%--------------------------------------------------------------------------+
if tropmodel_c == 1 %if USER SELECTs Combine(DRY+WET)TROPOSPHERIC MODEL
%--------------------------------------------------------------------------+ 

   %=======================================================================+
   %2.1***GET SELECTED COMBINE(DRY+WET) TROPO MODEL
   %-----------------------------------------------------------------------+
   getTROPOconts = cellstr(get(handles.pop_lTropo_combine,'String'));%get Tropo_combine popup menu contents as cell array
   tropoModel=getTROPOconts{get(handles.pop_lTropo_combine,'Value')};%get selected item from pop_lTropo_combine
   Modelv=get(handles.pop_lTropo_combine,'Value');%get selected item value from pop_lTropo_combine

   %-----------------------------------------------------------------------+
   %*****BLIND(LATITUDE,LONGITUDE,HEIGHT & TIME/DOY DEPENDENT) MODELS
   %-----------------------------------------------------------------------+
   
   %(1)*************WORKING ON GPT SERIES MODELS(GPT2,GPT2w,GPT3)
   
   %**********RETRIEVE STORED GRID VALUES/FILES FROM SELECTED MODELS
   %                    ---------------------------------
   
   %GET & SAVE GRID VALUES IF USER SELECTS ANY OF THE GPT MODELs 
   if any([strncmpi(tropoModel,'GPT2 (5° x 5°)',14),strncmpi(tropoModel,'GPT2w (1° x 1°)',14),strncmpi(tropoModel,'GPT2w (5° x 5°)',14),...
           strncmpi(tropoModel,'GPT3 (1° x 1°)',14),strncmpi(tropoModel,'GPT3 (5° x 5°)',14),Modelv == [12,13,14,15,16]])
      
       %*********CHECK IF STORED GRID VALUES FROM Tropo_combine popup menu
       if ~isempty(getappdata(0,'grid_c'))%if not empty([])
          
          %GET & SAVE GRID VALUES & GRID RESOLUTION
          setappdata(0,'grid_c',getappdata(0,'grid_c'))
          setappdata(0,'grid_res_c',getappdata(0,'grid_res_c')) 
          
       else %IF EMPTY([]),READ NEW GRID VALUES USING "readGPTgrid" fxn
           
           %CHECK GPT TROPO MODEL SELECTION TYPE
           if any([strncmpi(tropoModel,'GPT2 (5° x 5°)',14),Modelv == 12])
              grid_c     = readGPTgrid('gpt2_5.mat','GPT2',5);%GRID VALUES
              grid_res_c = 5;%GRID RESOLUTION
              
           elseif any([strncmpi(tropoModel,'GPT2w (1° x 1°)',14),Modelv == 13])
                  grid_c     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                  grid_res_c = 1;%GRID RESOLUTION
            
           elseif any([strncmpi(tropoModel,'GPT2w (5° x 5°)',14),Modelv == 14])
                  grid_c     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                  grid_res_c = 5;%GRID RESOLUTION
                  
           elseif any([strncmpi(tropoModel,'GPT3 (1° x 1°)',14),Modelv == 15])
                  grid_c     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                  grid_res_c = 1;%GRID RESOLUTION
                  
           elseif any([strncmpi(tropoModel,'GPT3 (5° x 5°)',14),Modelv == 16])
                  grid_c     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                  grid_res_c = 5;%GRID RESOLUTION
                  
           end %//if any([strncmpi(tropoModel,'GPT2 (5° x 5°)',14),Modelv == 12])
           
           %SAVE GRID VALUES & GRID RESOLUTION
           setappdata(0,'grid_c',grid_c)
           setappdata(0,'grid_res_c',grid_res_c)
           
       end %//if ~isempty(getappdata(0,'grid_c'))
        
       %****CHECK TIME VARIATION STATUS
       if ~isempty(getappdata(0,'Timevar_c'))
           
          %GET & SAVE TIME VARIATION STATUS
          setappdata(0,'Timevar_c',getappdata(0,'Timevar_c'))%time variation
          
       else
           %SET DEFAULT OF 1 & SAVE 
           setappdata(0,'Timevar_c',0)%time variation
           
       end %if ~isempty(getappdata(0,'Timevar_c'))
       
       %(2)*************WORKING ON VMF GRIDDED ZENITH DELAYS MODEL
        
       %GET & SAVE GRIDFILES IF USER SELECTS VMF ZTD MODEL
   elseif any([strncmpi(tropoModel,'VMF gridded ZTD',15),Modelv == 17])
          
          %GET & SAVE VMF GRID  FILES
          if any([~isempty(getappdata(0,'VMFgrids')),~isempty(getappdata(0,'VMF_grid_found'))]) 
        
             %*********SAVE STRUCT GRID FILES
             setappdata(0,'VMFgrids',getappdata(0,'VMFgrids')) 
           
             %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
             setappdata(0,'VMF_grid_found',getappdata(0,'VMF_grid_found'))
             
             %SAVE GRID FILE VERSION & TYPE
             setappdata(0,'VMFgrid_type',getappdata(0,'VMFgrid_type'))
             setappdata(0,'VMFgrid_res',getappdata(0,'VMFgrid_res'))  
           
          else 
            
              %CHECK IF VMF GRID FILEs ARE AVAILABLE IN goGPS
              %Call for the 'SearchVMFgrids.m' fxn
              [VMF_grid_found,VMFgrids] = SearchVMFgrids() ;
                        
              %*********SAVE STRUCT GRID FILES
              setappdata(0,'VMFgrids',VMFgrids)
         
              %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
              setappdata(0,'VMF_grid_found',VMF_grid_found)
              
              if VMF_grid_found == 1 %IF VMF grid FILES FOUND,LET USER CHOOSE GRID TYPE & RESOLUTION
                                     %I.E.['VMF1(2°x2.5°)','VMF3(1°x1°)','VMF3(5°x5°)']
                               
                 [VMF_type,grid_res] = VMFgrid_type(); 
                 
                 %SAVE GRID FILE VERSION & TYPE
                 setappdata(0,'VMFgrid_type',VMF_type)
                 setappdata(0,'VMFgrid_res',grid_res)
           
                 %*********CHECK IF VMF_type,grid_res ARE BOTH EMPTY
                 %NOTE:
                 %     VMF_type,grid_res ARE ASSIGNED EMPTY WHEN USER CANCEL SELECTION
                 %     OR CLOSE DIALOGUE BOX FOR THE SELECTION OF VMF GRID TYPE &
                 %     RESOLUTION.IF IT HAPPENS SO, THE SELECTED TROPO DELAY MODEL,
                 %     'VMF gridded ZTD' WOULD BE REPLACED WITH SAASTEMOIN MODEL
                 %     OR THE POPUP MENU WILL BE SET TO SAASTAMOINEN AS DEFAULT
        
                if all([isempty(VMF_type),isempty(grid_res)])%IF USER CANCEL SELECTION/CLOSE DIALOGUE BOX
            
                   %GET SOURCE OF MET PARAMETERS popup menu contents
                   getcontss = cellstr(get(handles.pop_source_metpara,'String'));%get popup menu contents as cell array
                   METpara=getcontss{get(handles.pop_source_metpara,'Value')};%get selected item from pop_lTropo_combine
                   metVAL=get(handles.pop_source_metpara,'Value');%Get selected item from pop_source_metpara
                   
                   beep %Give a beep sound
                   errmsg700{1}=sprintf('No VMF grid file(s) Type & Resolution selected .\n');
                   errmsg700{2}=sprintf('Try Again by selecting VMF grid file(s) Type & Resolution.\n');
                   errmsg700{3}=sprintf('i.e. Select the VMF gridded ZTD model again and choose from the popup\n');
                   warndlg(errmsg700,'VMF grid type Error','modal')
                    
                    
                   %SET COMBINE TROPO MODEL POPUP MENU VALUE TO 2(I.E.SAASTAMOINEN)
                   set(handles.pop_lTropo_combine,'Value',2)
              
                   %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
                   set(handles.cb_extract_VMF_ZTDs,'enable','on','value',0)
                   
                   %******SET MET PARAMETER SOURCE UICONTROLS VISIBLE/ENABLE ON
                   set(handles.pop_source_metpara,'enable','on')
                   set(handles.text_source_metpara,'enable','on')
    
                   if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
                                   any(metVAL==[7,8])])
        
                      if strcmpi(get(handles.pop_met_manual,'visible'),'Off') 
                         set(handles.pop_met_manual,'visible','On')
                         set(handles.text_metmanual_source,'visible','On')
                      end     
     
                   elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
                                       any(metVAL==[4,5])])
          
                          if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'off'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'off')])
                             set(handles.rb_grid_resolution_1,'Visible','on')
                             set(handles.rb_grid_resolution_5,'Visible','on')
                             set(handles.text_grid_resolution,'Visible','on')  
                          end     
          
                   end %//if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
                       %          any(metVAL==[7,8])])  
            
                  return
                  
                end  %//if all([isempty(VMF_type),isempty(grid_res)])
          

              end %if VMF_grid_found == 1
             
          end %//if any([~isempty(getappdata(0,'VMFgrids')),~isempty(getappdata(0,'VMF_grid_found'))]) 
          
          
        %GET & SAVE GRID VALUES IF USER SELECTS GTrop MODEL  
   elseif any([strncmpi(tropoModel,'GTrop [Sun et al 2019]',22),Modelv == 18])
                        
          %GET & SAVE GTrop GRID VALUES
          if ~isempty(getappdata(0,'gridV_GTrop'))
            
             setappdata(0,'gridV_GTrop',getappdata(0,'gridV_GTrop')) 
           
          else  
               %SEARCH FOR GTrop model GRID FILE(GTropCoefficient.mat) & DIRECTORY & LOAD VALUES               
               %Call the "SearchGTropgrid.m" fxn
               [~,gridV_GTrop] = SearchGTropgrid();%****LOADED GTrop MODEL COEFFICIENTS
       
               %SAVE GRID VALUES
               setappdata(0,'gridV_GTrop',gridV_GTrop)
          end        
                    
   end %//if any([strncmpi(tropoModel,'GPT2 (5° x 5°)',14),strncmpi(tropoModel,'GPT2w (1° x 1°)',14),strncmpi(tropoModel,'GPT2w (5° x 5°)',14),...
       %          strncmpi(tropoModel,'GPT3 (1° x 1°)',14),strncmpi(tropoModel,'GPT3 (5° x 5°)',14),Modelv == [12,13,14,15,16]])
 
       
%*******GET GOEID MODEL FOR MODELS THAT REQUIRE +VE ORTHOMETRIC HEIGHTS        
%--------------------------------------------------------------------------   
%NOTE:
%SAASTAMOINEN,MARINI,Askne & Nordius,UNB3m,EGNOS & MOPS MODELS REQUIRE +VE
%ORTHOMETRIC HEIGHT--------------------------------------------------------
%********GET goGPS DEFAULT GEOID MODEL (EGM2008)    
%CHECK SELECTED TROPO MODEL 

if any([strncmpi(tropoModel,'Saastamoinen',12),strncmpi(tropoModel,'Saastamoinen(Refined)',21),...
        strncmpi(tropoModel,'Askne & Nordius',15),strncmpi(tropoModel,'Marini',6),...
        strncmpi(tropoModel,'UNB3m',5),strncmpi(tropoModel,'EGNOS',5),strncmpi(tropoModel,'MOPS',4)])        

   if ~isempty(getappdata(0,'geoid'))
               
      %SAVE FOR FUTURE USE
      setappdata(0,'geoid',getappdata(0,'geoid')) 
              
   else  
       try %goGPS v0.5.2 beta1
          gs = Go_State.getInstance;
          geoid = gs.getRefGeoid(); 
      
       catch  %goGPS v0.4.3 
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
         
       end  %try 

    %SAVE FOR FUTURE USE
    setappdata(0,'geoid',geoid)
   
   end %//if ~isempty(getappdata(0,'geoid'))
   
end %//if any([strncmpi(tropoModel,'Saastamoinen',12),strncmpi(tropoModel,'Saastamoinen(Refined)',21),...
    %          strncmpi(tropoModel,'Askne & Nordius',15),strncmpi(tropoModel,'Marini',6),...
    %          strncmpi(tropoModel,'UNB3m',5),strncmpi(tropoModel,'EGNOS',5),strncmpi(tropoModel,'MOPS',4)])
 
   %SAVE SELECTED COMBINE TROPO MODEL
   setappdata(0,'Tropo_Cmodel',tropoModel)
   setappdata(0,'TROPO_Cmodel',tropoModel)
   setappdata(0,'Tropo_Cmodel_v',Modelv)
   
   %===================================COMBINE(DRY+WET) TROPOSPHERIC MODELS+
   %-----------------------------------------------------------------------

%--------------------------------------------------------------------------+   
elseif tropmodel_s == 1 %if USER SELECTs SEPARATE(DRY,WET)TROPOSPHERIC MODEL
%--------------------------------------------------------------------------+

       %===================================================================+
       %2.2*****************GET SELECTED SEPARATE TROPO MODEL
       %-------------------------------------------------------------------+
      
       %HYDROSTATIC model
       getcontsH = cellstr(get(handles.pop_ITropo_hydrostatic,'String'));%get Tropo_hydrostatic popup menu contents as cell array
       dryModel=getcontsH{get(handles.pop_ITropo_hydrostatic,'Value')};%get selected item from pop_Tropo_hydrostatic
       Modelv_H=get(handles.pop_ITropo_hydrostatic,'Value');%get selected item value from pop_ITropo_hydrostatic
       
       %WET model
       getcontsW = cellstr(get(handles.pop_ITropo_wet,'String'));%get Tropo_wet popup menu contents as cell array
       wetModel=getcontsW{get(handles.pop_ITropo_wet,'Value')};%get selected item from pop_Tropo_wet
       Modelv_W=get(handles.pop_ITropo_wet,'Value');%get selected item value from pop_ITropo_wet

       %-------------------------------------------------------------------+
       %*****BLIND(LATITUDE,LONGITUDE,HEIGHT & TIME/DOY DEPENDENT) MODELS
       %-------------------------------------------------------------------+
   
       %(1)*************WORKING ON GPT SERIES MODELS(GPT2,GPT2w,GPT3)
       
       %GET & SAVE GRIDVALUES IF USER SELECTS ANY OF THE HYDROSTATIC GPT MODELs 
        if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),Modelv_H==[10,11,12,13,14]]) & ...   
           ~any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                 strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),Modelv_W==[15,16,17,18,19]])  
          
          
          %****CHECK STORED GRID VALUES FROM HYDROSTATIC MODEL POPUP MENU
          if ~isempty(getappdata(0,'grid_h'))%if not empty([])
          
             %GET & SAVE GRID VALUES & GRID RESOLUTION
             setappdata(0,'grid_h',getappdata(0,'grid_h'))
             setappdata(0,'grid_res_h',getappdata(0,'grid_res_h'))
             
          else %IF EMPTY([]),READ NEW GRID VALUES USING "readGPTgrid" fxn 
               
               %CHECK GPT TROPO MODEL SELECTION TYPE
               if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),Modelv_H == 10])
                  grid_h     = readGPTgrid('gpt2_5.mat','GPT2',5);%GRID VALUES
                  grid_res_h = 5;%GRID RESOLUTION
              
           elseif any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),Modelv_H == 11])
                  grid_h     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                  grid_res_h = 1;%GRID RESOLUTION
            
           elseif any([strncmpi(dryModel,'GPT2w (5° x 5°)',14),Modelv_H == 12])
                  grid_h     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                  grid_res_h = 5;%GRID RESOLUTION
                  
           elseif any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),Modelv_H == 13])
                  grid_h     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                  grid_res_h = 1;%GRID RESOLUTION
                  
           elseif any([strncmpi(dryModel,'GPT3 (5° x 5°)',14),Modelv_H == 14])
                  grid_h     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                  grid_res_h = 5;%GRID RESOLUTION
                  
               end %//if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),Modelv_H == 10])
           
           %SAVE GRID VALUES & GRID RESOLUTION
           setappdata(0,'grid_h',grid_h)
           setappdata(0,'grid_res_h',grid_res_h)
           
          end  %//if ~isempty(getappdata(0,'grid_h'))
          
          %SET WET GRID FILE & RESOLUTION TO EMPTY([])
          setappdata(0,'grid_w',[])
          setappdata(0,'grid_res_w',[])
          
          %****CHECK TIME VARIATION STATUS
          if ~isempty(getappdata(0,'Timevar_dry'))
           
             %GET & SAVE TIME VARIATION STATUS
             setappdata(0,'Timevar_dry',getappdata(0,'Timevar_dry'))%time variation
          
          else  
              %SET DEFAULT OF 1 & SAVE 
              setappdata(0,'Timevar_dry',0)%time variation
           
          end  %if ~isempty(getappdata(0,'Timevar_dry'))    
             
          setappdata(0,'Timevar_wet',1)%set time variation for wet model to static
          
              
       %GET & SAVE GRIDFILE IF USER SELECTS ANY OF THE WET GPT MODELs 
       elseif any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                   strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),Modelv_W==[15,16,17,18,19]]) & ...
              ~any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                    strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),Modelv_H==[10,11,12,13,14]])  
    
              %****CHECK STORED GRID VALUES FROM WET MODEL POPUP MENU
              if ~isempty(getappdata(0,'grid_w'))%if not empty([])
          
                 %GET & SAVE GRID VALUES & GRID RESOLUTION
                 setappdata(0,'grid_w',getappdata(0,'grid_w'))
                 setappdata(0,'grid_res_w',getappdata(0,'grid_res_w'))
             
              else %IF EMPTY([]),READ NEW GRID VALUES USING "readGPTgrid" fxn 
               
                   %CHECK GPT TROPO MODEL SELECTION TYPE
                   if any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),Modelv_W == 15])
                      grid_w     = readGPTgrid('gpt2_5.mat','GPT2',5);%GRID VALUES
                      grid_res_w = 5;%GRID RESOLUTION
              
                   elseif any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),Modelv_W == 16])
                          grid_w     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                          grid_res_w = 1;%GRID RESOLUTION
            
                   elseif any([strncmpi(wetModel,'GPT2w (5° x 5°)',14),Modelv_W == 17])
                          grid_w     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                          grid_res_w = 5;%GRID RESOLUTION
                  
                   elseif any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),Modelv_W == 18])
                          grid_w     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                          grid_res_w = 1;%GRID RESOLUTION
                  
                   elseif any([strncmpi(wetModel,'GPT3 (5° x 5°)',14),Modelv_W == 19])
                          grid_w     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                          grid_res_w = 5;%GRID RESOLUTION
                  
                   end  %//if any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),Modelv_W == 15])
           
                   %SAVE GRID VALUES & GRID RESOLUTION
                   setappdata(0,'grid_w',grid_w)
                   setappdata(0,'grid_res_w',grid_res_w)
           
              end    %//if ~isempty(getappdata(0,'grid_w'))
          
              %SET HYDROSTATIC GRID VALUES & RESOLUTION TO EMPTY([])
              setappdata(0,'grid_h',[])
              setappdata(0,'grid_res_h',[])
          
              %****CHECK TIME VARIATION STATUS
              if ~isempty(getappdata(0,'Timevar_wet'))
           
                 %GET & SAVE TIME VARIATION STATUS
                 setappdata(0,'Timevar_wet',getappdata(0,'Timevar_wet'))%time variation
          
              else   
                  %SET DEFAULT OF 1 & SAVE 
                  setappdata(0,'Timevar_wet',0)%time variation
           
              end  %if ~isempty(getappdata(0,'Timevar_w'))    
             
              setappdata(0,'Timevar_dry',1)%set time variation for dry model to static
  
        %GET & SAVE GRIDFILES IF USER SELECTS BOTH THE DRY & WET GPT MODELs 
        elseif any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                   strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),Modelv_H==[10,11,12,13,14]]) & ...   
              any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                   strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),Modelv_W==[15,16,17,18,19]])  
            
                %IF DRY & WET MODELs ARE THE SAME & ARE GPT MODELs
               if strncmpi(dryModel,wetModel,14)%Compare 1st 14 characters
                   
                 %GET ANY OF THE SAVED GRIDs FILE(e.g.HYDROSTATIC)& RESOLUTION
                 if ~isempty(getappdata(0,'grid_h'))%if not empty([])
          
                    %GET & SAVE GRID VALUES & GRID RESOLUTION
                    setappdata(0,'grid_h',getappdata(0,'grid_h'))
                    setappdata(0,'grid_res_h',getappdata(0,'grid_res_h'))
             
                    setappdata(0,'grid_w',getappdata(0,'grid_h'))
                    setappdata(0,'grid_res_w',getappdata(0,'grid_res_h'))
                    
                 else %IF EMPTY([]),READ NEW GRID VALUES USING "readGPTgrid" fxn 
               
                      %CHECK GPT TROPO MODEL SELECTION TYPE
                      if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),Modelv_H == 10])
                         grid_h     = readGPTgrid('gpt2_5.mat','GPT2',5);%GRID VALUES
                         grid_res_h = 5;%GRID RESOLUTION
              
                      elseif any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),Modelv_H == 11])
                             grid_h     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                             grid_res_h = 1;%GRID RESOLUTION
            
                      elseif any([strncmpi(dryModel,'GPT2w (5° x 5°)',14),Modelv_H == 12])
                             grid_h     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                             grid_res_h = 5;%GRID RESOLUTION
                  
                      elseif any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),Modelv_H == 13])
                             grid_h     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                             grid_res_h = 1;%GRID RESOLUTION
                  
                      elseif any([strncmpi(dryModel,'GPT3 (5° x 5°)',14),Modelv_H == 14])
                             grid_h     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                             grid_res_h = 5;%GRID RESOLUTION
                  
                      end  %//if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),Modelv_H == 10])
           
                      %SAVE GRID VALUES & GRID RESOLUTION 
                      %HYDROSTATIC 
                      setappdata(0,'grid_h',grid_h)
                      setappdata(0,'grid_res_h',grid_res_h)
                      %WET
                      setappdata(0,'grid_w',grid_h)
                      setappdata(0,'grid_res_w',grid_res_h)
           
                 end   %//if ~isempty(getappdata(0,'grid_h'))
          
     
                %****CHECK TIME VARIATION STATUS
                if ~isempty(getappdata(0,'Timevar_dry'))
           
                   %GET & SAVE TIME VARIATION STATUS
                   setappdata(0,'Timevar_dry',getappdata(0,'Timevar_dry'))%time variation
                   setappdata(0,'Timevar_wet',getappdata(0,'Timevar_wet'))%time variation for wet delay
                else   
                    %SET DEFAULT VALUE OF 1 & SAVE 
                    setappdata(0,'Timevar_dry',0)%time variation For dry model 
                    setappdata(0,'Timevar_wet',0)%time variation for wet model 
           
                end %if ~isempty(getappdata(0,'Timevar_dry'))    
             
            
               else %IF DRY & WET MODELs ARE NOT THE SAME & ARE GPT MODELs 
                  
                   %****CHECK STORED GRID VALUES FROM HYDROSTATIC MODEL POPUP MENU
                   if ~isempty(getappdata(0,'grid_h'))%if not empty([])
          
                      %GET & SAVE GRID VALUES & GRID RESOLUTION
                      setappdata(0,'grid_h',getappdata(0,'grid_h'))
                      setappdata(0,'grid_res_h',getappdata(0,'grid_res_h'))
             
                   else %IF EMPTY([]),READ NEW GRID VALUES USING "readGPTgrid" fxn 
               
                       %CHECK GPT TROPO MODEL SELECTION TYPE
                       if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),Modelv_H == 10])
                          grid_h     = readGPTgrid('gpt2_5.mat','GPT2',5);%GRID VALUES
                          grid_res_h = 5;%GRID RESOLUTION
              
                       elseif any([strncmpi(dryModel,'GPT2w (1° x 1°)',14),Modelv_H == 11])
                              grid_h     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                              grid_res_h = 1;%GRID RESOLUTION
            
                       elseif any([strncmpi(dryModel,'GPT2w (5° x 5°)',14),Modelv_H == 12])
                              grid_h     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                              grid_res_h = 5;%GRID RESOLUTION
                  
                       elseif any([strncmpi(dryModel,'GPT3 (1° x 1°)',14),Modelv_H == 13])
                              grid_h     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                              grid_res_h = 1;%GRID RESOLUTION
                  
                       elseif any([strncmpi(dryModel,'GPT3 (5° x 5°)',14),Modelv_H == 14])
                              grid_h     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                              grid_res_h = 5;%GRID RESOLUTION
                  
                       end %//if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),Modelv_H == 10])
           
                       %SAVE GRID VALUES & GRID RESOLUTION
                       setappdata(0,'grid_h',grid_h)
                       setappdata(0,'grid_res_h',grid_res_h)
           
                   end %//if ~isempty(getappdata(0,'grid_h'))
          
                  %****CHECK TIME VARIATION STATUS
                  if ~isempty(getappdata(0,'Timevar_dry'))
           
                     %GET & SAVE TIME VARIATION STATUS
                     setappdata(0,'Timevar_dry',getappdata(0,'Timevar_dry'))%time variation
          
                  else   
                      %SET DEFAULT OF 1 & SAVE 
                       setappdata(0,'Timevar_dry',0)%time variation
           
                  end %if ~isempty(getappdata(0,'Timevar_dry'))    
                  
                  
                 %****CHECK STORED GRID VALUES FROM WET MODEL POPUP MENU
                 if ~isempty(getappdata(0,'grid_w'))%if not empty([])
          
                    %GET & SAVE GRID VALUES & GRID RESOLUTION
                    setappdata(0,'grid_w',getappdata(0,'grid_w'))
                    setappdata(0,'grid_res_w',getappdata(0,'grid_res_w'))
             
                 else %IF EMPTY([]),READ NEW GRID VALUES USING "readGPTgrid" fxn 
               
                      %CHECK GPT TROPO MODEL SELECTION TYPE
                      if any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),Modelv_W == 15])
                         grid_w     = readGPTgrid('gpt2_5.mat','GPT2',5);%GRID VALUES
                         grid_res_w = 5;%GRID RESOLUTION
              
                      elseif any([strncmpi(wetModel,'GPT2w (1° x 1°)',14),Modelv_W == 16])
                             grid_w     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                             grid_res_w = 1;%GRID RESOLUTION
            
                      elseif any([strncmpi(wetModel,'GPT2w (5° x 5°)',14),Modelv_W == 17])
                             grid_w     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                             grid_res_w = 5;%GRID RESOLUTION
                  
                      elseif any([strncmpi(wetModel,'GPT3 (1° x 1°)',14),Modelv_W == 18])
                             grid_w     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                             grid_res_w = 1;%GRID RESOLUTION
                  
                      elseif any([strncmpi(wetModel,'GPT3 (5° x 5°)',14),Modelv_W == 19])
                             grid_w     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                             grid_res_w = 5;%GRID RESOLUTION
                  
                      end %//if any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),Modelv_W == 15])
           
                      %SAVE GRID VALUES & GRID RESOLUTION
                      setappdata(0,'grid_w',grid_w)
                      setappdata(0,'grid_res_w',grid_res_w)
           
                 end %//if ~isempty(getappdata(0,'grid_w'))
          
                 %****CHECK TIME VARIATION STATUS
                 if ~isempty(getappdata(0,'Timevar_wet'))
           
                    %GET & SAVE TIME VARIATION STATUS
                    setappdata(0,'Timevar_wet',getappdata(0,'Timevar_wet'))%time variation
          
                 else    
                     %SET DEFAULT VALUE OF 1 & SAVE 
                     setappdata(0,'Timevar_wet',0)%time variation
           
                 end %if ~isempty(getappdata(0,'Timevar_w')) 
                         

               end%//if strncmpi(dryModel,wetModel,14) 
              
               
               %(2)*************WORKING ON VMF GRIDDED ZENITH DELAYS MODEL  
               
              %IF VMF ZPD MODEL IS SELECTED 
        elseif any([any([any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv_H == 15]),...
                         any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv_W == 20])]),...
                    all([any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv_H == 15]),...
                         any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv_W == 20])])])
                     
               %GET & SAVE VMF GRID  FILES
               if all([~isempty(getappdata(0,'VMFgrids')),~isempty(getappdata(0,'VMF_grid_found'))]) 
        
                  %*********SAVE STRUCT GRID FILES
                  setappdata(0,'VMFgrids',getappdata(0,'VMFgrids')) 
           
                  %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
                  setappdata(0,'VMF_grid_found',getappdata(0,'VMF_grid_found'))
                  
                  %SAVE GRID FILE VERSION & TYPE
                  setappdata(0,'VMFgrid_type',getappdata(0,'VMFgrid_type'))
                  setappdata(0,'VMFgrid_res',getappdata(0,'VMFgrid_res'))  
           
               else 
            
                   %CHECK IF VMF GRID FILEs ARE AVAILABLE IN goGPS
                   %Call for the 'SearchVMFgrids.m' fxn
                   [VMF_grid_found,VMFgrids] = SearchVMFgrids(handles) ;
                        
                   %*********SAVE STRUCT GRID FILES
                   setappdata(0,'VMFgrids',VMFgrids)
         
                   %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
                   setappdata(0,'VMF_grid_found',VMF_grid_found)
                   
                   if VMF_grid_found == 1 %IF VMF grid FILES FOUND,LET USER CHOOSE GRID TYPE & RESOLUTION
                                          %I.E.['VMF1(2°x2.5°)','VMF3(1°x1°)','VMF3(5°x5°)']
                               
                      [VMF_type,grid_res] = VMFgrid_type(); 
           
                      %SAVE GRID FILE VERSION & TYPE
                      setappdata(0,'VMFgrid_type',VMF_type)
                      setappdata(0,'VMFgrid_res',grid_res)
                      
                      %*********CHECK IF VMF_type,grid_res ARE BOTH EMPTY
                      %NOTE:
                      %     VMF_type,grid_res ARE ASSIGNED EMPTY WHEN USER CANCEL SELECTION
                      %     OR CLOSE DIALOGUE BOX FOR THE SELECTION OF VMF GRID TYPE &
                      %     RESOLUTION.IF IT HAPPENS SO, THE SELECTED TROPO DELAY MODEL,
                      %     'VMF gridded ZTD' WOULD BE REPLACED WITH SAASTEMOIN MODEL
                      %     OR THE POPUP MENU WILL BE SET TO SAASTAMOINEN AS DEFAULT
        
                     if all([isempty(VMF_type),isempty(grid_res)])%IF USER CANCEL SELECTION/CLOSE DIALOGUE BOX
                   
                        %GET SOURCE OF MET PARAMETERS popup menu contents
                        getcontss = cellstr(get(handles.pop_source_metpara,'String'));%get popup menu contents as cell array
                        METpara=getcontss{get(handles.pop_source_metpara,'Value')};%get selected item from pop_lTropo_combine
                        metVAL=get(handles.pop_source_metpara,'Value');%Get selected item from pop_source_metpara
                    
                        beep %Give a beep sound
                        errmsg701{1}=sprintf('No VMF grid file(s) Type & Resolution selected .\n');
                        errmsg701{2}=sprintf('Try Again by selecting VMF grid file(s) Type & Resolution.\n');
                        errmsg701{3}=sprintf('i.e. Select the VMF gridded ZTD model again and choose from the popup\n');
                        warndlg(errmsg701,'VMF grid type Error','modal')
                        
                        if any([strncmpi(dryModel,'VMF ZHD',7),Modelv_H == 15])
                            
                           %SET DRY TROPO MODEL POPUP MENU VALUE TO 2(I.E.SAASTAMOINEN)
                           set(handles.pop_ITropo_hydrostatic,'Value',1)
                        end  
                        
                        if any([strncmpi(wetModel,'VMF ZWD',7),Modelv_W == 20])
                            
                           %SET WET TROPO MODEL POPUP MENU VALUE TO 1(I.E.SAASTAMOINEN)
                           set(handles.pop_ITropo_wet,'Value',1)
                        end  
              
                        %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
                        set(handles.cb_extract_VMF_ZTDs,'enable','on','value',0)
                   
                        %******SET MET PARAMETER SOURCE UICONTROLS VISIBLE/ENABLE ON
                        set(handles.pop_source_metpara,'enable','on')
                        set(handles.text_source_metpara,'enable','on')
    
                        if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
                                        any(metVAL==[7,8])])
        
                           if strcmpi(get(handles.pop_met_manual,'visible'),'Off') 
                              set(handles.pop_met_manual,'visible','On')
                              set(handles.text_metmanual_source,'visible','On')
                           end      
     
                        elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
                                            any(metVAL==[4,5])])
          
                               if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'off'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'off')])
                                  set(handles.rb_grid_resolution_1,'Visible','on')
                                  set(handles.rb_grid_resolution_5,'Visible','on')
                                  set(handles.text_grid_resolution,'Visible','on')  
                               end      
          
                        end  %//if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
                             %          any(metVAL==[7,8])])  
                  
                        return
                  
                     end 
              

                   end  %if VMF_grid_found == 1
             
               end %//if all([~isempty(getappdata(0,'VMFgrids')),~isempty(getappdata(0,'VMF_grid_found'))])       
                     
              %IF GTrop MODEL IS SELECTED               
        elseif any([any([any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),Modelv_H == 16]),...
                         any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),Modelv_W == 21])]),...
                    all([any([strncmpi(dryModel,'GTrop [Sun et al 2019]',22),Modelv_H == 16]),...
                         any([strncmpi(wetModel,'GTrop [Sun et al 2019]',22),Modelv_W == 21])])])
               
               %GET & SAVE GTrop GRID VALUES
               if ~isempty(getappdata(0,'gridV_GTrop'))
            
                  setappdata(0,'gridV_GTrop',getappdata(0,'gridV_GTrop')) 
           
               else 
                   %SEARCH FOR GTrop model GRID FILE(GTropCoefficient.mat) & DIRECTORY & LOAD VALUES               
                   %Call the "SearchGTropgrid.m" fxn
                   [~,gridV_GTrop] = SearchGTropgrid();%****LOADED GTrop MODEL COEFFICIENTS
       
                   %SAVE GRID VALUES
                   setappdata(0,'gridV_GTrop',gridV_GTrop)
               end       
                                                                 
       end %//if any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
           %          strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),Modelv_H==[10,11,12,13,14]]) & ...   
           %    ~any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
           %          strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),Modelv_W==[15,16,17,18,19]])
        
           
       %*******GET GOEID MODEL FOR MODELS THAT REQUIRE +VE ORTHOMETRIC HEIGHTS  
       %-------------------------------------------------------------------  
       %NOTE"
       %SAASTAMOINEN,DAVIS ET AL,Askne & Nordius,UNB3m,EGNOS & MOPS MODELS 
       %REQUIRE +VE ORTHOMETRIC HEIGHT-------------------------------------
       %********GET goGPS DEFAULT GEOID MODEL (EGM2008)    
       %CHECK SELECTED TROPO MODEL
       if any([strncmpi(dryModel,'Saastamoinen',12),strncmpi(dryModel,'Davis et al)',12),strncmpi(dryModel,'Askne & Nordius',15),...
               strncmpi(dryModel,'UNB3m',5),strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),strncmpi(wetModel,'Saastamoinen',12),...
               strncmpi(wetModel,'Askne & Nordius',15),strncmpi(wetModel,'UNB3m',5),strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4)])          
          
           %LOAD GEOID MODEL[%USING goGPS DEFAULT GEOID MODEL (EGM2008)] 
           if ~isempty(getappdata(0,'geoid'))
               
              %SAVE FOR FUTURE USE
              setappdata(0,'geoid',getappdata(0,'geoid')) 
              
           else
               
               try %goGPS v0.5.2 beta1
                  gs = Go_State.getInstance;
                  geoid = gs.getRefGeoid(); 
      
               catch  %goGPS v0.4.3 
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
         
               end  %try

              %SAVE FOR FUTURE USE
             setappdata(0,'geoid',geoid)
             
           end %/if ~isempty(getappdata(0,'geoid'))    
       
       end %//if any([strncmpi(dryModel,'Saastamoinen',12),strncmpi(dryModel,'Davis et al)',12),strncmpi(dryModel,'Askne & Nordius',15),...
           %          strncmpi(dryModel,'UNB3m',5),strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),strncmpi(wetModel,'Saastamoinen',12),...
           %          strncmpi(wetModel,'Askne & Nordius',15),strncmpi(wetModel,'UNB3m',5),strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4)])
           
        %SAVE SELECTED SEPARATE TROPO MODELs
        setappdata(0,'dryModel',dryModel)
        setappdata(0,'wetModel',wetModel)
        setappdata(0,'DryModel',dryModel)
        setappdata(0,'WetModel',wetModel)
        setappdata(0,'dryModel_v',Modelv_H)
        setappdata(0,'wetModel_v',Modelv_W)
        
       %=================================SEPARATE TROPOSPHERIC MODELS(DRY, WET)
       %-----------------------------------------------------------------------

       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       %(2.3)****************GET SELECTED MAPPING FUNCTION(MF)
       %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
       
       if get(handles.rb_different_MF,'Value') == 1 %if Use other Mapping Models is selected
           
          %HYDROSTATIC MF
          getcontsHMF = cellstr(get(handles.pop_IHydrostic_MF,'String'));%get Hydrostatic MF popup menu contents as cell array
          MFh_model=getcontsHMF{get(handles.pop_IHydrostic_MF,'Value')};%get selected item from pop_IHydrostic_MF
          
          %WET MF
          getcontsWMF = cellstr(get(handles.pop_IWet_MF,'String'));%get Wet MF popup menu contents as cell array
          MFw_model=getcontsWMF{get(handles.pop_IWet_MF,'Value')};%get selected item from pop_IWet_MF
          
           %----------------------------------------------------------------+
          %****WORKING ON VIENNA MAPPING FUNCTION(VMF) MODELS (VMF1/VMF3)
          %----------------------------------------------------------------+
          
          %GET & SAVE GRIDFILE IF USER SELECTS ANY OF THE HYDROSTATIC VMF MODELs 
          if any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),...
                  strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12)]) & ...
             ~any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),...
                   strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12)])
              
              %****CHECK STORED GRID VALUES FROM HYDROSTATIC MODEL POPUP MENU
              if ~isempty(getappdata(0,'grid_MFh'))%if not empty([])
          
                 %GET & SAVE GRID VALUES & GRID RESOLUTION
                 setappdata(0,'grid_MFh',getappdata(0,'grid_MFh'))
                 setappdata(0,'gridRES_MFh',getappdata(0,'gridRES_MFh'))
             
              else %IF EMPTY([]),READ NEW GRID VALUES USING "readGPTgrid" fxn  
                  
                   %CHECK VMF MODEL SELECTION TYPE
                   if strncmpi(MFh_model,'VMF1(1° x 1°)',12) 
                      grid_MFh     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                      gridRES_MFh  = 1;%GRID RESOLUTION
                      
                   elseif strncmpi(MFh_model,'VMF1(5° x 5°)',12)
                          grid_MFh     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                          gridRES_MFh  = 5;%GRID RESOLUTION 
               
                   elseif strncmpi(MFh_model,'VMF3(1° x 1°)',12)
                          grid_MFh     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                          gridRES_MFh  = 1;%GRID RESOLUTION
                
                   elseif strncmpi(MFh_model,'VMF3(5° x 5°)',12)
                          grid_MFh     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                          gridRES_MFh  = 5;%GRID RESOLUTION     
                   end 
               
                   %SAVED HYDROSTATIC VMF GRID VALUES & RESOLUTION                                            
                   setappdata(0,'grid_MFh',grid_MFh)
                   setappdata(0,'gridRES_MFh',gridRES_MFh)
                   
              end %//if ~isempty(getappdata(0,'grid_MFh'))
             
              %SET WET VMF GRID FILE & RESOLUTION TO EMPTY([])
              setappdata(0,'grid_MFw',[])
              setappdata(0,'gridRES_MFw',[])
             
          %GET & SAVE GRIDFILE IF USER SELECTS ANY OF THE WET VMF MODELs   
          elseif any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),...
                      strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12)]) & ...
                 ~any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),...
                       strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12)])
                
                %****CHECK STORED GRID VALUES FROM HYDROSTATIC MODEL POPUP MENU
                if ~isempty(getappdata(0,'grid_MFw'))%if not empty([])
          
                   %GET & SAVE GRID VALUES & GRID RESOLUTION
                   setappdata(0,'grid_MFw',getappdata(0,'grid_MFw'))
                   setappdata(0,'gridRES_MFw',getappdata(0,'gridRES_MFw'))
             
                else %IF EMPTY([]),READ NEW GRID VALUES USING "readGPTgrid" fxn  
                  
                     %CHECK VMF MODEL SELECTION TYPE
                     if strncmpi(MFw_model,'VMF1(1° x 1°)',12) 
                        grid_MFw     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                        gridRES_MFw  = 1;%GRID RESOLUTION
                      
                     elseif strncmpi(MFw_model,'VMF1(5° x 5°)',12)
                            grid_MFw     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                            gridRES_MFw  = 5;%GRID RESOLUTION 
               
                     elseif strncmpi(MFw_model,'VMF3(1° x 1°)',12)
                            grid_MFw     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                            gridRES_MFw  = 1;%GRID RESOLUTION
                
                     elseif strncmpi(MFw_model,'VMF3(5° x 5°)',12)
                            grid_MFw     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                            gridRES_MFw  = 5;%GRID RESOLUTION     
                     end  
               
                   %SAVED WET VMF GRID VALUES & RESOLUTION                                            
                   setappdata(0,'grid_MFw',grid_MFw)
                   setappdata(0,'gridRES_MFw',gridRES_MFw)
                   
                end %//if ~isempty(getappdata(0,'grid_MFw'))
             
               %SET HYDROSTATIC VMF GRID FILE & RESOLUTION TO EMPTY([])
               setappdata(0,'grid_MFh',[])
               setappdata(0,'gridRES_MFh',[])  
                         
          %GET & SAVE GRIDFILES IF USER SELECTS BOTH HYDROSTATIC & WET VMF MODELs        
          elseif any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),...
                      strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12)]) & ...
                 any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),...
                      strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12)])
           
                  
                 %IF HYDROSTATIC & WET MODELs ARE THE SAME & ARE VMF MODELs
               if strncmpi(MFh_model,MFw_model,12)%Compare 1st 12 characters
                 
                  %GET ANY OF THE SAVED GRIDs FILE(e.g.HYDROSTATIC)& RESOLUTION 
                  %HYDROSTATIC MF
                  %****CHECK STORED GRID VALUES FROM HYDROSTATIC MODEL POPUP MENU
                  if ~isempty(getappdata(0,'grid_MFh'))%if not empty([])
          
                     %GET & SAVE GRID VALUES & GRID RESOLUTION
                     setappdata(0,'grid_MFh',getappdata(0,'grid_MFh'))
                     setappdata(0,'gridRES_MFh',getappdata(0,'gridRES_MFh'))
             
                  else %IF EMPTY([]),READ NEW GRID VALUES USING "readGPTgrid" fxn  
                  
                       %CHECK VMF MODEL SELECTION TYPE
                       if strncmpi(MFh_model,'VMF1(1° x 1°)',12) 
                          grid_MFh     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                          gridRES_MFh  = 1;%GRID RESOLUTION
                      
                       elseif strncmpi(MFh_model,'VMF1(5° x 5°)',12)
                              grid_MFh     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                              gridRES_MFh  = 5;%GRID RESOLUTION 
               
                       elseif strncmpi(MFh_model,'VMF3(1° x 1°)',12)
                              grid_MFh     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                              gridRES_MFh  = 1;%GRID RESOLUTION
                
                       elseif strncmpi(MFh_model,'VMF3(5° x 5°)',12)
                              grid_MFh     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                              gridRES_MFh  = 5;%GRID RESOLUTION     
                       end  
               
                       %SAVE GRID VALUES & GRID RESOLUTION 
                       %HYDROSTATIC                                            
                       setappdata(0,'grid_MFh',grid_MFh)
                       setappdata(0,'gridRES_MFh',gridRES_MFh)
                       %WET MF
                       setappdata(0,'grid_MFw',grid_MFh)
                       setappdata(0,'gridRES_MFw',gridRES_MFh)
                   
                  end %//if ~isempty(getappdata(0,'grid_MFh'))
                     
               else %IF HYDROSTATIC & WET MF MODELs ARE NOT THE SAME & ARE VMF MODELs 
                 
                    %****CHECK STORED GRID VALUES FROM HYDROSTATIC MODEL POPUP MENU
                    if ~isempty(getappdata(0,'grid_MFh'))%if not empty([])
                       
                       %GET & SAVE GRID VALUES & GRID RESOLUTION
                       setappdata(0,'grid_MFh',getappdata(0,'grid_MFh'))
                       setappdata(0,'gridRES_MFh',getappdata(0,'gridRES_MFh'))
             
                    else %IF EMPTY([]),READ NEW GRID VALUES USING "readGPTgrid" fxn 
                        
                         %CHECK VMF MODEL SELECTION TYPE
                         if strncmpi(MFh_model,'VMF1(1° x 1°)',12) 
                            grid_MFh     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                            gridRES_MFh  = 1;%GRID RESOLUTION
                      
                         elseif strncmpi(MFh_model,'VMF1(5° x 5°)',12)
                                grid_MFh     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                                gridRES_MFh  = 5;%GRID RESOLUTION 
               
                         elseif strncmpi(MFh_model,'VMF3(1° x 1°)',12)
                                grid_MFh     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                                gridRES_MFh  = 1;%GRID RESOLUTION
                
                         elseif strncmpi(MFh_model,'VMF3(5° x 5°)',12)
                                grid_MFh     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                                gridRES_MFh  = 5;%GRID RESOLUTION     
                         end    
               
                         %SAVE GRID VALUES & GRID RESOLUTION                                                
                         setappdata(0,'grid_MFh',grid_MFh)
                         setappdata(0,'gridRES_MFh',gridRES_MFh)
                   
                    end %//if ~isempty(getappdata(0,'grid_MFh'))
               
                   
                   %****CHECK STORED GRID VALUES FROM WET MODEL POPUP MENU
                   if ~isempty(getappdata(0,'grid_MFw'))%if not empty([])
          
                      %GET & SAVE GRID VALUES & GRID RESOLUTION
                      setappdata(0,'grid_MFw',getappdata(0,'grid_MFw'))
                      setappdata(0,'gridRES_MFw',getappdata(0,'gridRES_MFw'))
             
                   else %IF EMPTY([]),READ NEW GRID VALUES USING "readGPTgrid" fxn  
                  
                        %CHECK VMF MODEL SELECTION TYPE
                        if strncmpi(MFw_model,'VMF1(1° x 1°)',12) 
                           grid_MFw     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                           gridRES_MFw  = 1;%GRID RESOLUTION
                      
                        elseif strncmpi(MFw_model,'VMF1(5° x 5°)',12)
                               grid_MFw     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                               gridRES_MFw  = 5;%GRID RESOLUTION 
               
                        elseif strncmpi(MFw_model,'VMF3(1° x 1°)',12)
                               grid_MFw     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                               gridRES_MFw  = 1;%GRID RESOLUTION
                
                        elseif strncmpi(MFw_model,'VMF3(5° x 5°)',12)
                               grid_MFw     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                               gridRES_MFw  = 5;%GRID RESOLUTION     
                        end   
               
                        %SAVED WET VMF GRID VALUES & RESOLUTION                                            
                        setappdata(0,'grid_MFw',grid_MFw)
                        setappdata(0,'gridRES_MFw',gridRES_MFw)
                   
                   end %//if ~isempty(getappdata(0,'grid_MFw'))
                   

               end %//if strncmpi(MFh_model,MFw_model,12)
               
 
          end %//if any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),...
              %          strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12)]) & ...
              %     ~any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),...
              %           strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12)])
           
              
          %GET & SAVE VMF GRID FILES [eg:'VMFG_20180101.H00','VMFG_20180101.H06','VMFG_20180101.H12','VMFG_20180101.H18']
          %----------------------------------------------------------------
          %NOTE: 
          %     WHEN VMF GRID FILES ARE NOT FOUND, THE SYSTEM RESORT TO USING
          %     THE GPT SERIES MODELS GRID VALUES FOR THE COMPUTATION OF VMF
          %----------------------------------------------------------------
          if any([any([any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),...
                            strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12)]),...
                       any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),...
                            strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12)])]),...
                  all([any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),...
                            strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12)]),...
                       any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),...
                            strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12)])])])
                   
                        %GET & SAVE VMF GRID  FILES
             if all([~isempty(getappdata(0,'VMFgrids')),~isempty(getappdata(0,'VMF_grid_found'))]) 
        
                %*********SAVE STRUCT GRID FILES
                setappdata(0,'VMFgrids',getappdata(0,'VMFgrids')) 
           
                %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
                setappdata(0,'VMF_grid_found',getappdata(0,'VMF_grid_found'))
                
                %SAVE GRID FILE VERSION & TYPE
                setappdata(0,'VMFgrid_type_mf',getappdata(0,'VMFgrid_type_mf'))
                setappdata(0,'VMFgrid_res_mf',getappdata(0,'VMFgrid_res_mf'))      
           
             else 
            
                 %CHECK IF VMF GRID FILEs ARE AVAILABLE IN goGPS
                 %Call for the 'SearchVMFgrids.m' fxn
                 [VMF_grid_found,VMFgrids] = SearchVMFgrids(handles) ;
                        
                 %*********SAVE STRUCT GRID FILES
                 setappdata(0,'VMFgrids',VMFgrids)
         
                 %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
                 setappdata(0,'VMF_grid_found',VMF_grid_found)
                 
                 if VMF_grid_found == 1 %IF VMF grid FILES FOUND
                           
                    if any([any([any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12)]),...
                                 any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12)])]),...
                            all([any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12)]),...
                                 any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12)])])])
                            
                       VMF_type_mf = 'VMF1';
                       grid_res_mf = 2; 
          
                    elseif  any([any([strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(1° x 1°)',12)]),...
                            all([strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(1° x 1°)',12)])])
              
                            VMF_type_mf = 'VMF3';
                             grid_res_mf = 1;
              
                    elseif  any([any([strncmpi(MFh_model,'VMF3(5° x 5°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12)]),...
                            all([strncmpi(MFh_model,'VMF3(5° x 5°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12)])])
                      
                            VMF_type_mf = 'VMF3';
                            grid_res_mf = 5;
              
                    end  %//if any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12)])
             
                    
                    %SAVE GRID FILE VERSION & TYPE
                    setappdata(0,'VMFgrid_type_mf',VMF_type_mf)
                    setappdata(0,'VMFgrid_res_mf',grid_res_mf)
                       
                 end %//if VMF_grid_found == 1 
             
             end %//if all([~isempty(getappdata(0,'VMFgrids')),~isempty(getappdata(0,'VMF_grid_found'))]) 
              
          end
                
          %SAVE SELECTED MAPPING FUNCTION MODELs
          setappdata(0,'MFh_model',MFh_model)
          setappdata(0,'MFw_model',MFw_model)
          setappdata(0,'MFh_Model',MFh_model)
          setappdata(0,'MFw_Model',MFw_model)
      %===========================================MAPPING FUNCTIONS
      %--------------------------------------------------------------------+
      
       end  %//if get(handles.rb_different_MF,'Value') == 1
       
       %SAVE STATE OF OPTIONs FOR USE OF MODELs / OTHER MAPPING FUNCTIONS
       setappdata(0,'option_model_MF',get(handles.rb_model_MF,'Value'))%save option button for models MF
       setappdata(0,'option_different_MF',get(handles.rb_different_MF,'Value'))%save option button for different MF models
            
end  %//if tropmodel_c == 1

% % %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
% % %*************GET MAPPING FUNCTION FOR PPP TROPO MODELLING
% % %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
% % if all([get(handles.flag_tropo,'value') == 1, strcmpi(get(handles.flag_tropo,'enable'),'on')])
% %    
% %    %RETRIEVE & RESAVED ALL SAVED DATA FROM PPP MAPPING FUNCTION(MF) SELECTION DIALOGUE
% %    PPPmf_POPval = getappdata(0,'PPPmf_POPval');%MF INDEX FROM POPUP MENU
% %   
% %    %GET SAVE MAPPING FUNCTION MODELS
% %    if all([~isempty(getappdata(0,'MFh_model_ppp')),~isempty(getappdata(0,'MFw_model_ppp'))])
% %       
% %       setappdata(0,'MFh_model_ppp',getappdata(0,'MFh_model_ppp'))
% %       setappdata(0,'MFw_model_ppp',getappdata(0,'MFw_model_ppp'))
% %       
% %       if PPPmf_POPval == 3 %IF USER SELECTS VMF MODEL
% %      
% %          %GET & SAVE STRUCT GRID FILES
% %          setappdata(0,'VMFgrids_ppp',getappdata(0,'VMFgrids_ppp'))
% %          
% %          %GET & SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
% %          setappdata(0,'VMF_grid_found_ppp',getappdata(0,'VMF_grid_found_ppp')) 
% %         
% %          %SAVE GRID FILE VERSION & TYPE WITH EMPTY MATRIX([])
% %          setappdata(0,'VMFgrid_type_ppp',getappdata(0,'VMFgrid_type_ppp'))
% %          setappdata(0,'VMFgrid_res_ppp',getappdata(0,'VMFgrid_res_ppp')) 
% %        
% %          %GET & SAVED GPT2w/GPT3 GRID VALUES & RESOLUTION
% %          %HYDROSTATIC
% %          setappdata(0,'grid_MFh_ppp',getappdata(0,'grid_MFh_ppp'))
% %          setappdata(0,'gridRES_MFh_ppp',getappdata(0,'gridRES_MFh_ppp'))
% %       
% %          %WET
% %          setappdata(0,'grid_MFw_ppp',getappdata(0,'grid_MFw_ppp'))
% %          setappdata(0,'gridRES_MFw_ppp',getappdata(0,'gridRES_MFw_ppp'))
% %          
% %          %GPT MODEL
% %          setappdata(0,'GPTmodel_ppp',getappdata(0,'GPTmodel_ppp'))
% %       
% %       end  %//if PPPmf_POPval == 3
% %       
% %    else  
% %        if PPPmf_POPval == 3 %IF USER SELECTS VMF MODEL
% %     
% %           %CHECK IF VMF GRID FILEs ARE AVAILABLE IN goGPS
% %           %Call for the 'SearchVMFgrids.m' fxn
% %           [VMF_grid_found,VMFgrids] = SearchVMFgrids(handles); 
% %                         
% %           %*********SAVE STRUCT GRID FILES
% %           setappdata(0,'VMFgrids_ppp',VMFgrids)
% %          
% %           %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
% %           setappdata(0,'VMF_grid_found_ppp',VMF_grid_found)
% %    
% %           if VMF_grid_found == 1 %IF VMF grid FILES FOUND,LET USER CHOOSE GRID TYPE & RESOLUTION
% %                                %I.E.['VMF1(2°x2.5°)','VMF3(1°x1°)','VMF3(5°x5°)']
% %                                
% %              [VMF_type,grid_res] = VMFgrid_type(); 
% %       
% %              %SAVE GRID FILE VERSION & TYPE
% %              setappdata(0,'VMFgrid_type_ppp',VMF_type)
% %              setappdata(0,'VMFgrid_res_ppp',grid_res)
% %       
% %              grid_MF_ppp = [];%ASSIGN EMPTY([]) MATRIX TO GPT series grid values
% %              gridRES_MF_ppp = [];
% %              GPTmodel_ppp = []; %GPT MODEL WITH EMPTY([]) MATRIX
% %       
% %              %*********CHECK IF VMF_type,grid_res ARE BOTH EMPTY
% %              %NOTE:
% %              %     VMF_type,grid_res ARE ASSIGNED EMPTY WHEN USER CANCEL SELECTION
% %              %     OR CLOSE DIALOGUE BOX FOR THE SELECTION OF VMF GRID TYPE &
% %              %     RESOLUTION.IF IT HAPPENS SO, THE SELECTED TROPO DELAY MODEL,
% %              %     'VMF gridded ZTD' WOULD BE REPLACED WITH SAASTEMOIN MODEL
% %              %     OR THE POPUP MENU WILL BE SET TO SAASTAMOINEN AS DEFAULT
% %       
% %              if all([isempty(VMF_type),isempty(grid_res)])%IF USER CANCEL SELECTION/CLOSE DIALOGUE BOX
% %           
% %                 beep %Give a beep sound
% %                 warnmsg1{1}=sprintf('VMF grid type & resolution has not been selected .\n');
% %                 warnmsg1{2}=sprintf('Try Again and select VMF grid type & resolution.\n');        
% %                 warndlg(warnmsg1,'VMF grid file Error','modal')
% %          
% %                 %READ FROM GPT3 GRID FILE(using GPT3 1° x 1° model)
% %                 grid_MF_ppp=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDVALUES  
% %                 gridRES_MF_ppp = 1;%GRID RESOLUTION
% %          
% %                 MFh_ppp = 'VMF3(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
% %                 MFw_ppp = 'VMF3(1° x 1°)';%WET MAPPING FUNCTION MODEL 
% %                 GPTmodel_ppp = 'GPT3'; %GPT3 MODEL 
% %                 
% %                 if isempty(grid_MF_ppp) %if GPT3 model grid is unavailable
% %             
% %                    %READ FROM GPT2w GRID FILE(using GPT2w 1° x 1° model)
% %                    grid_MF_ppp = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
% %                    gridRES_MF_ppp = 1;%GRID RESOLUTION
% %          
% %                    MFh_ppp = 'VMF1(1° x 1°)'; %HYDROSTATIC MAPPING FUNCTION MODEL 
% %                    MFw_ppp = 'VMF1(1° x 1°)'; %WET MAPPING FUNCTION MODEL 
% %                    GPTmodel_ppp = 'GPT2w'; %GPT2w MODEL 
% %                    
% %                 end          
% %             
% %              else %IF USER DOES NOT CANCEL SELECTION/CLOSE DIALOGUE BOX
% %                   if any([strcmpi(VMF_type,'VMF1'),strfind(VMF_type,'VMF1')])
% %                
% %                      MFh_ppp = 'VMF1(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
% %                      MFw_ppp = 'VMF1(1° x 1°)';%WET MAPPING FUNCTION MODEL
% %               
% %                   elseif any([strcmpi(VMF_type,'VMF3'),strfind(VMF_type,'VMF3')])
% %                   
% %                          if grid_res == 1
% %                       
% %                             MFh_ppp = 'VMF3(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
% %                             MFw_ppp = 'VMF3(1° x 1°)';%WET MAPPING FUNCTION MODEL
% %                      
% %                          elseif grid_res == 5
% %                       
% %                                 MFh_ppp = 'VMF3(5° x 5°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
% %                                 MFw_ppp = 'VMF3(5° x 5°)';%WET MAPPING FUNCTION MODEL
% %                        
% %                          end   
% %                       
% %                   end %//if any([strcmpi(VMF_type,'VMF1'),strfind(VMF_type,'VMF1')])
% %            
% %              end  %//if all([isempty(VMF_type),isempty(grid_res)])
% %       
% %           else %IF VMF grid FILES NOT FOUND
% %             
% %                %SAVE GRID FILE VERSION & TYPE WITH EMPTY MATRIX([])
% %                setappdata(0,'VMFgrid_type_ppp',[])
% %                setappdata(0,'VMFgrid_res_ppp',[])
% %         
% %                beep %Give a beep sound
% %                errmsg3{1}=sprintf('No VMF grid file(s) found in goGPS directory.\n');
% %                errmsg3{2}=sprintf('System will resort to GPT2w / GPT3 grid file(s) instead.\n');
% %                errmsg3{3}=sprintf('Mean while, VMF grid(s) can as well be downloaded from : http://vmf.geo.tuwien.ac.at/trop_products/GRID/\n');
% %                warndlg(errmsg3,'VMF grid file Error','modal')
% %         
% %               %READ FROM GPT3 GRID FILE(using GPT3 1° x 1° model)
% %               grid_MF_ppp=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE 
% %               gridRES_MF_ppp = 1;%GRID RESOLUTION
% %         
% %               MFh_ppp = 'VMF3(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
% %               MFw_ppp = 'VMF3(1° x 1°)';%WET MAPPING FUNCTION MODEL 
% %               GPTmodel_ppp = 'GPT3'; %GPT3 MODEL 
% %               
% %               if isempty(grid_MF_ppp) %if GPT3 model grid is unavailable
% %             
% %                  %READ FROM GPT2w GRID FILE(using GPT2w 1° x 1° model)
% %                  grid_MF_ppp = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
% %          
% %                  MFh_ppp = 'VMF1(1° x 1°)'; %HYDROSTATIC MAPPING FUNCTION MODEL 
% %                  MFw_ppp = 'VMF1(1° x 1°)'; %WET MAPPING FUNCTION MODEL 
% %                  GPTmodel_ppp = 'GPT2w'; %GPT2w MODEL 
% %           
% %               end                
% %    
% %           end   %//if VMF_grid_found == 1
% %    
% %           %SAVED GPT2w/GPT3 GRID VALUES & RESOLUTION
% %           %HYDROSTATIC
% %           setappdata(0,'grid_MFh_ppp',grid_MF_ppp)
% %           setappdata(0,'gridRES_MFh_ppp',gridRES_MF_ppp)
% %       
% %           %WET
% %           setappdata(0,'grid_MFw_ppp',grid_MF_ppp)
% %           setappdata(0,'gridRES_MFw_ppp',gridRES_MF_ppp)
% %           
% %           %GPT MODEL
% %           setappdata(0,'GPTmodel_ppp',GPTmodel_ppp)
% %      
% %        elseif PPPmf_POPval == 2 %IF USER SELECTS GMF
% %               MFh_ppp = 'GMF';
% %               MFw_ppp = 'GMF';
% %        
% %        elseif PPPmf_POPval == 1 %IF USER SELECTS NMF
% %               MFh_ppp = 'NMF';
% %               MFw_ppp = 'NMF';   
% %     
% %        end %//if PPPmf_POPval == 3  
% %      
% %        %SAVE MAPPING FUNCTION MODELS
% %        setappdata(0,'MFh_model_ppp',MFh_ppp)
% %        setappdata(0,'MFw_model_ppp',MFw_ppp)
% %          
% %    end %//if all([~isempty(getappdata(0,'MFh_model_ppp')),~isempty(getappdata(0,'MFw_model_ppp'))]) 
% %      
% % end %//if all([get(handles.flag_tropo,'value') == 1, strcmpi(get(handles.flag_tropo,'enable'),'on')])    
% %       

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%(3)****************GET SELECTED SOURCE MET PARAMETERS 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%GET ENABILITY STATE OF MET PARAMETERS popup menu
metENABLE = get(handles.pop_source_metpara,'Enable');%Enable state (on or off)

if strcmpi(metENABLE,'on')%if enable 'on'
    
 %GET SOURCE OF MET PARAMETERS popup menu contents
 getconts = cellstr(get(handles.pop_source_metpara,'String'));%get popup menu contents as cell array
 metPARA=getconts{get(handles.pop_source_metpara,'Value')};%get selected item from pop_lTropo_combine
 metVAL=get(handles.pop_source_metpara,'Value');%Get selected item from pop_source_metpara

 %GET GRID RESOLUTION FOR GPT2w/GPT3 MODELS
 grid_res1=get(handles.rb_grid_resolution_1,'value');%for 1° resolution
 grid_res5=get(handles.rb_grid_resolution_5,'value');%for 5° resolution
   
 %CHECK WHICH MODEL IS SELECTED FOR THE COMPUTATION OF MET PARAMETERS
 if any([strcmpi(metPARA,'1.Standard Meteorological Parameters') isequal(metVAL,1)] ) 
    METpara='Standard';
   
 elseif any([strcmpi(metPARA,'2.Generate from GPT model') isequal(metVAL,2)] ) 
       METpara='GPT';
   
       %****GET SUPPLEMENTARY FILES FOR SUPPLEMENTARY MET PARAMETERS
       %-------------------------------------------------------------------
       %GPT MODEL PROVIDES ONLY PRESSURE(P) & TEMPERATURE(T).WE NEED OTHER
       %MET PARAMETERS LIKE WATER VAPOUR PARTIAL PRESSURE(e) TO SUPPLEMENT
       %WHAT THE GPT MODEL PROVIDES----------------------------------------
       
       %FIRST CHECK IF ALREADY STORED GRID FILE IS EMPTY([]) OR NOT
       if ~isempty(getappdata(0,'grid_met_s'))
          
          %GET & STORE GRID VALUES,GRID RESOLUTION & GPT MODEL TYPE(GPT2w/GPT3) 
          setappdata(0,'grid_met_s',getappdata(0,'grid_met_s'))
          setappdata(0,'grid_res_met_s',getappdata(0,'grid_res_met_s'))
          setappdata(0,'grid_METmodel',getappdata(0,'grid_METmodel'))
          
       else   
           %READ FROM GRID FILE(using GPT2w model)
           grid_met_ss = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
           grid_res_met_ss = 5;%GRID RESOLUTION 
           grid_MET_model = 'GPT2w';
       
           if isempty(grid_met_ss) %if GPT2w model grid is unavailable
              %READ FROM GRID FILE(using GPT3 model)
              grid_met_ss=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE 
              grid_MET_model = 'GPT3';
           end   
   
           setappdata(0,'grid_met_s',grid_met_ss) %SAVE GRIDFILE FOR MET PARAMETER
           setappdata(0,'grid_res_met_s',grid_res_met_ss)%SAVE GRID RESOLUTION FOR MET PARAMETER
           setappdata(0,'grid_METmodel',grid_MET_model)%SAVE GRID MODEL
           
       end %if ~isempty(getappdata(0,'grid_met_s'))
 
       %FOR GPT2,GPT2w, & GPT3 MODELs, RETRIEVE A CELL ARRAY OF GRID FILES & RESOLUTION
       %-------------------------------------------------------------------------------
 elseif any([strcmpi(metPARA,'3.Generate from GPT2 model') isequal(metVAL,3)] ) 
        METpara='GPT2';
        
       grid_Met=getappdata(0,'grid_met'); %Gridfile
       grid_res_Met=5; %Grid Resolution
       
 elseif any([strcmpi(metPARA,'4.Generate from GPT2w model') isequal(metVAL,4)] ) %USER SELECTS GPT2w MODEL
        METpara='GPT2w';
       if get(handles.rb_grid_resolution_1,'Value') == 1 & get(handles.rb_grid_resolution_5,'Value') == 0
           
         if any([~isempty(getappdata(0,'grid_met')) ~isempty(getappdata(0,'grid_res_met'))]) &  getappdata(0,'grid_res_met') == 1   
            grid_Met=getappdata(0,'grid_met'); %Gridfile
            grid_res_Met=1; %Grid Resolution
         else
             grid_Met=getappdata(0,'grid_met_rb1'); %Gridfile
             grid_res_Met=1;%Grid Resolution
         end
         
       elseif get(handles.rb_grid_resolution_1,'Value') == 0 & get(handles.rb_grid_resolution_5,'Value') == 1
             if any([~isempty(getappdata(0,'grid_met')) ~isempty(getappdata(0,'grid_res_met'))]) &  getappdata(0,'grid_res_met') == 5   
                grid_Met=getappdata(0,'grid_met'); %Gridfile
                grid_res_Met=5; %Grid Resolution
             else 
                 grid_Met=getappdata(0,'grid_met_rb5'); %Gridfile
                 grid_res_Met=5;%Grid Resolution
             end 
       end 
       
 elseif any([strcmpi(metPARA,'5.Generate from GPT3 model') isequal(metVAL,5)] )%USER SELECTS GPT3 MODEL
        METpara='GPT3';
       
        if get(handles.rb_grid_resolution_1,'Value') == 1 & get(handles.rb_grid_resolution_5,'Value') == 0
           
         if any([~isempty(getappdata(0,'grid_met')) ~isempty(getappdata(0,'grid_res_met'))]) &  getappdata(0,'grid_res_met') == 1   
            
             grid_Met=getappdata(0,'grid_met'); %Gridfile
             grid_res_Met=1; %Grid Resolution
         else
             grid_Met=getappdata(0,'grid_met_rb1'); %Gridfile
             grid_res_Met=1;%Grid Resolution
         end
         
       elseif get(handles.rb_grid_resolution_1,'Value') == 0 & get(handles.rb_grid_resolution_5,'Value') == 1
             
              if any([~isempty(getappdata(0,'grid_met')) ~isempty(getappdata(0,'grid_res_met'))]) &  getappdata(0,'grid_res_met') == 5   
                
                 grid_Met=getappdata(0,'grid_met'); %Gridfile
                 grid_res_Met=5; %Grid Resolution
             else 
                 grid_Met=getappdata(0,'grid_met_rb5'); %Gridfile
                 grid_res_Met=5;%Grid Resolution
             end 
        end
       
 elseif any([strcmpi(metPARA,'6.Generate from UNB3m model') isequal(metVAL,6)] )%USER SELECTS UNB3m MODEL
        METpara='UNB3m';
       
 elseif any([strcmpi(metPARA,'7.Enter Meteorological Parameters') isequal(metVAL,7)] )
      
     
     if ~isempty(getappdata(0,'metPARA_csv'))
         
        METpara=getappdata(0,'metPARA_csv');
        
     else
         %OPEN A DIALOGUE FOR USER FILE IMPORT
         importMETPARA_csv();
         METpara=getappdata(0,'metPARA_csv');
     end
       
 elseif any([strcmpi(metPARA,'8.Import from csv file [T P RH]') isequal(metVAL,8)] )
        METpara=getappdata(0,'metPARA_table');
        
   
 elseif any([strcmpi(metPARA,'9.MET file'),isequal(metVAL,9)])%if user opt for use of MET file
                                                                  %MET file contains Met parameters
                                                                  %observed/measured at GNSS site
        METpara = 'MET file';
     
        %GET OBSERVED METEO DATA @ GNSS SITE
        gs = Go_State.getInstance();
        state = gs.getCurrentSettings();
        
        if exist('ini_settings_file', 'var')
           state.importIniFile(ini_settings_file);
        end 
        
        md = Meteo_Data(state.getMetFile());%METEO DATA
        
        if any([~(md.isValid()),md.isValid()]) %IF FILE IS INVALID, %USE GPT2W 1° x 1° MET DATA AS ALTERNATIVE 
        
           if ~isempty(getappdata(0,'grid_METfile'))%IF STORED GRID VALUES IS NOT EMPTY
              grid_METfile = getappdata(0,'grid_METfile');
              grid_res_METfile = getappdata(0,'grid_res_METfile');
          
           else %*****READ GRID FILE
                grid_METfile =readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                grid_res_METfile = 1;%GRID RESOLUTION
            
           end
           
           %SAVE GRID VALUES & GRID RESOLUTION
           setappdata(0,'grid_METfile',grid_METfile)
           setappdata(0,'grid_res_METfile',grid_res_METfile)   
           
        end %//if any([~(md.isValid()),md.isValid()])      
        
        
 end %//if any([strcmpi(metPARA,'1.Standard Meteorological Parameters') isequal(metVAL,1)] ) 
 
 
if any([isequal(metVAL,7) isequal(metVAL,8)]) 
    
   %GET SOURCE OF USER ENTERED METEOROLOGICAL PARAMETERS
   metSource=get(handles.pop_met_manual,'Value'); 
    
  if metSource == 1
    
    if tropmodel_c == 1
    
       if ~any([strncmpi(tropoModel,'GPT2 (5° x 5°)',14),strncmpi(tropoModel,'GPT2w (1° x 1°)',14),strncmpi(tropoModel,'GPT2w (5° x 5°)',14),...
                strncmpi(tropoModel,'GPT3 (1° x 1°)',14),strncmpi(tropoModel,'GPT3 (5° x 5°)',14),strncmpi(tropoModel,'UNB3m',5),...
                strncmpi(tropoModel,'EGNOS',5),strncmpi(tropoModel,'MOPS',4),strncmpi(tropoModel,'VMF gridded ZTD',15),...
                strncmpi(tropoModel,'GTrop [Sun et al 2019]',22)])
    
           beep %Give a beep sound
           msg=sprintf('Select Source of Entered / Imported Meteorological Parameters');
           warndlg(msg,'Entered / Imported Meteorological Source', 'modal');
          return
          
       end
       
    elseif tropmodel_s == 1
   
           if any([~any([strncmpi(dryModel,'GPT2 (5° x 5°)',14),strncmpi(dryModel,'GPT2w (1° x 1°)',14),strncmpi(dryModel,'GPT2w (5° x 5°)',14),...
                    strncmpi(dryModel,'GPT3 (1° x 1°)',14),strncmpi(dryModel,'GPT3 (5° x 5°)',14),strncmpi(dryModel,'UNB3m',5),...
                    strncmpi(dryModel,'EGNOS',5),strncmpi(dryModel,'MOPS',4),strncmpi(dryModel,'VMF gridded ZHD',15),strncmpi(dryModel,'GTrop [Sun et al 2019]',22)]),...   
              ~any([strncmpi(wetModel,'GPT2 (5° x 5°)',14),strncmpi(wetModel,'GPT2w (1° x 1°)',14),strncmpi(wetModel,'GPT2w (5° x 5°)',14),...
                    strncmpi(wetModel,'GPT3 (1° x 1°)',14),strncmpi(wetModel,'GPT3 (5° x 5°)',14),strncmpi(wetModel,'UNB3m',5),...
                    strncmpi(wetModel,'EGNOS',5),strncmpi(wetModel,'MOPS',4),strncmpi(wetModel,'VMF gridded ZWD',15),strncmpi(wetModel,'GTrop [Sun et al 2019]',22)])])
            
             beep %Give a beep sound
             msg=sprintf('Select Source of Entered / Imported Meteorological Parameters');
             warndlg(msg,'Entered / Imported Meteorological Source', 'modal');
            return
   
           end 
           
    end
   
else
     if metSource == 2
    
        METsource = 'WEATHER STATION';
   
     elseif  metSource == 3
    
           METsource = 'NWP';
                
     elseif  metSource == 4
    
           METsource = 'IN-SITU';    
     end
     
     %SAVE SOURCE OF USER INPUT MET PARAMETERS
     setappdata(0,'METsource',METsource)
     
   end %if metSource == 1
     
end %//if get(handles.rb_source_ENTERmetpara,'Value')==1 & metSource == 1

%SAVE SELECTED MET PARAMETER & GRID RESOLUTION OPTIONs
setappdata(0,'METpara',METpara)%SAVE SELECTED MET PARAMETER OPTION
setappdata(0,'metVAL',metVAL)%SAVE SELECTED MET PARAMETER OPTION(Numeric)

if any([any(strcmpi(METpara,{'GPT2','GPT2w','GPT3'})),any(metVAL==[3,4,5])])

   if all([~isempty(grid_Met),~isempty(grid_res_Met)])
       
      setappdata(0,'grid_Met',grid_Met) %SAVE GRIDFILE FOR MET PARAMETER
      setappdata(0,'grid_res_Met',grid_res_Met)%SAVE GRID RESOLUTION FOR MET PARAMETER 
      
   else %IF EMPTY([]),READ NEW GRID VALUES USING "readGPTgrid" fxn 
       
       %IF MET MODEL IS GPT2
       if any([strcmpi(METpara,'GPT2'),metVAL==3])
          grid_Met     = readGPTgrid('gpt2_5.mat','GPT2',5);%GRID VALUES
          grid_res_Met = 5;%GRID RESOLUTION
          
          %IF MET MODEL IS GPT2w
       elseif any([strcmpi(METpara,'GPT2w'),metVAL==4])
              
              %**********CHECK GRID REOLUTION TYPE(1° OR 5°)
              %IF 1° GRID RESOLUTION IS SELECTED
              if all([grid_res1 == 1,grid_res5 == 0])
                 grid_Met     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
                 grid_res_Met = 1;%GRID RESOLUTION 
                 
              %IF 5° GRID RESOLUTION IS SELECTED    
              elseif all([grid_res1 == 0,grid_res5 == 1])
                     grid_Met     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                     grid_res_Met = 5;%GRID RESOLUTION
                     
              end
            
             %IF MET MODEL IS GPT3
       elseif any([strcmpi(METpara,'GPT3'),metVAL==5])
              
              %***********CHECK GRID REOLUTION TYPE(1° OR 5°)
              %IF 1° GRID RESOLUTION IS SELECTED
              if all([grid_res1 == 1,grid_res5 == 0])
                 grid_Met     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                 grid_res_Met = 1;%GRID RESOLUTION 
                 
                %IF 5° GRID RESOLUTION IS SELECTED  
              elseif all([grid_res1 == 0,grid_res5 == 1])
                     grid_Met     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                     grid_res_Met = 5;%GRID RESOLUTION
                     
              end
              
       end %//if any([strcmpi(METpara,'GPT2'),metVAL==3])
       
       setappdata(0,'grid_Met',grid_Met) %SAVE GRIDFILE FOR MET PARAMETER
       setappdata(0,'grid_res_Met',grid_res_Met)%SAVE GRID RESOLUTION FOR MET PARAMETER
       
       
   end %//if all([~isempty(grid_Met),~isempty(grid_res_Met)])
          
   %****CHECK TIME VARIATION STATUS
   if ~isempty(getappdata(0,'Timevar_met'))
           
      %GET & SAVE TIME VARIATION STATUS
      setappdata(0,'Timevar_met',getappdata(0,'Timevar_met'))%time variation
        
   else   
       %SET DEFAULT VALUE OF 1 & SAVE 
       setappdata(0,'Timevar_met',0)%time variation
           
   end %if ~isempty(getappdata(0,'Timevar_Met'))      
    
end %//if any([any(strcmpi(METpara,{'GPT2','GPT2w','GPT3'})),any(metVAL==[3,4,5])])

%**************GET SUPPLEMENTARY MET PARAMETERS
%NOTE:---------------------------------------------------------------------
%GPT,GPT2 MODELS & USER INPUT MET PARAMETERS DOESN'T PROVIDE ALL THE NEEDED
%METEOROLOGICAL PARAMETERS(e.g.:Tm,lambda).SOME OF THE TROPO
%MODELS(e.g.:'Askne & Nordius' WET MODEL) REQUIRE Mean temperature of water 
%vapor(Tm) & water vapour `lapse rate' or decrease factor(lambda). WE NEED
%TO SUPPLEMENT THE MET PARAMETERS ANYTIME ANY OF THESE MODELS IS SELECTED
%--------------------------------------------------------------------------
if tropmodel_c == 1
    
   if strncmpi(tropoModel,'Askne & Nordius',15)%IF SELECTED COMBINE TROP MODEL IS 'Askne & Nordius'
  
      if any([strcmpi(metPARA,'1.Standard Meteorological Parameters'),...
              strcmpi(metPARA,'3.Generate from GPT2 model'),strcmpi(metPARA,'7.Enter Meteorological Parameters'),...
              strcmpi(metPARA,'8.Import from csv file [T P RH]'),metVAL==[1,3,7,8]])
        
          
         %FIRST CHECK IF ALREADY STORED GRID FILE IS EMPTY([]) OR NOT
         if ~isempty(getappdata(0,'grid_met_s'))
          
            %GET & STORE GRID VALUES,GRID RESOLUTION & GPT MODEL TYPE(GPT2w/GPT3) 
            setappdata(0,'grid_met_s',getappdata(0,'grid_met_s'))
            setappdata(0,'grid_res_met_s',getappdata(0,'grid_res_met_s'))
            setappdata(0,'grid_METmodel',getappdata(0,'grid_METmodel'))
          
         else    
             %READ FROM GRID FILE(using GPT2w model)
             grid_met_ss = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
             grid_res_met_ss = 5;%GRID RESOLUTION 
             grid_MET_model = 'GPT2w';
       
             if isempty(grid_met_ss) %if GPT2w model grid is unavailable
                %READ FROM GRID FILE(using GPT3 model)
                grid_met_ss=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE 
                grid_MET_model = 'GPT3';
             end    
   
                setappdata(0,'grid_met_s',grid_met_ss) %SAVE GRIDFILE FOR MET PARAMETER
                setappdata(0,'grid_res_met_s',grid_res_met_ss)%SAVE GRID RESOLUTION FOR MET PARAMETER
                setappdata(0,'grid_METmodel',grid_MET_model)%SAVE GRID MODEL
           
         end  %if ~isempty(getappdata(0,'grid_met_s')) 
          
      end %//if any([strcmpi(metPARA,'1.Standard Meteorological Parameters'),...
          %          strcmpi(metPARA,'3.Generate from GPT2 model'),strcmpi(metPARA,'7.Enter Meteorological Parameters'),...
          %          strcmpi(metPARA,'8.Import from csv file [T P RH]'),metVAL==[1,3,7,8]])
      
   end %//if strncmpi(tropoModel,'Askne & Nordius',15)
     
elseif tropmodel_s == 1
    
       if strncmpi(wetModel,'Askne & Nordius',15)
  
          if any([strcmpi(metPARA,'1.Standard Meteorological Parameters'),...
                  strcmpi(metPARA,'3.Generate from GPT2 model'),strcmpi(metPARA,'7.Enter Meteorological Parameters'),...
                  strcmpi(metPARA,'8.Import from csv file [T P RH]'),metVAL==[1,3,7,8]])
        
             %FIRST CHECK IF ALREADY STORED GRID FILE IS EMPTY([]) OR NOT
             if ~isempty(getappdata(0,'grid_met_s'))
          
                %GET & STORE GRID VALUES,GRID RESOLUTION & GPT MODEL TYPE(GPT2w/GPT3) 
                setappdata(0,'grid_met_s',getappdata(0,'grid_met_s'))
                setappdata(0,'grid_res_met_s',getappdata(0,'grid_res_met_s'))
                setappdata(0,'grid_METmodel',getappdata(0,'grid_METmodel'))
          
             else     
                 %READ FROM GRID FILE(using GPT2w model)
                 grid_met_ss = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRIDFILE
                 grid_res_met_ss = 5;%GRID RESOLUTION 
                 grid_MET_model = 'GPT2w';
       
                 if isempty(grid_met_ss) %if GPT2w model grid is unavailable
                    %READ FROM GRID FILE(using GPT3 model)
                    grid_met_ss=readGPTgrid('gpt3_5.mat','GPT3',5);%GRIDFILE 
                    grid_MET_model = 'GPT3';
                 end     
   
                 setappdata(0,'grid_met_s',grid_met_ss) %SAVE GRIDFILE FOR MET PARAMETER
                 setappdata(0,'grid_res_met_s',grid_res_met_ss)%SAVE GRID RESOLUTION FOR MET PARAMETER
                 setappdata(0,'grid_METmodel',grid_MET_model)%SAVE GRID MODEL
           
             end   %if ~isempty(getappdata(0,'grid_met_s')) 
             
          end %//if any([strcmpi(metPARA,'1.Standard Meteorological Parameters'),...
              %          strcmpi(metPARA,'3.Generate from GPT2 model'),strcmpi(metPARA,'7.Enter Meteorological Parameters'),...
              %          strcmpi(metPARA,'8.Import from csv file [T P RH]'),metVAL==[1,3,7,8]])
          
       end %if strncmpi(wetModel,'Askne & Nordius',15)
       
end %//if tropmodel_c == 1


end %//if strcmpi(metENABLE,'on')

%SAVE ENABLE STATE(on or off)
setappdata(0,'metENABLE',metENABLE)
    
%================================================METEOROLOGICAL PARAMETERS+
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-    
 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%ADD VMF GRIDDED ZENITH DELAYS TO OUTPUT FILE IF USER SELECTED MODEL TYPE 
%IS NOT  THE VMF GRIDDED ZENITH DELAYS(I.E VMF ZTD,VMF ZHD,VMF ZWD
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
getvmfZTDs_state = get(handles.cb_extract_VMF_ZTDs,'Value');%Get state of checkbox(1 or 0)
getvmfZTDs_enable = get(handles.cb_extract_VMF_ZTDs,'enable');%Get state of checkbox ENABILITY

if all([getvmfZTDs_state == 1,strcmpi(getvmfZTDs_enable,'on')])

   %GET & SAVE VMF GRID  FILES
   if all([~isempty(getappdata(0,'VMFgrids')),~isempty(getappdata(0,'VMF_grid_found'))]) 
        
      %*********SAVE STRUCT GRID FILES
      setappdata(0,'VMFgrids',getappdata(0,'VMFgrids')) 
           
      %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
      setappdata(0,'VMF_grid_found',getappdata(0,'VMF_grid_found'))
      
      %SAVE GRID FILE VERSION & TYPE
      setappdata(0,'VMFgrid_type',getappdata(0,'VMFgrid_type'))
      setappdata(0,'VMFgrid_res',getappdata(0,'VMFgrid_res'))
           
   else  
            
       %CHECK IF VMF GRID FILEs ARE AVAILABLE IN goGPS
       %Call for the 'SearchVMFgrids.m' fxn
       [VMF_grid_found,VMFgrids] = SearchVMFgrids(handles) ;
                        
       %*********SAVE STRUCT GRID FILES
       setappdata(0,'VMFgrids',VMFgrids)
         
      %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
      setappdata(0,'VMF_grid_found',VMF_grid_found)
      
      if VMF_grid_found == 0 %IF VMF grid FILES NOT FOUND
           
         beep %Give a beep sound
         
         %PRINT SOME MESSAGE TO USER
         fprintf('\n\nNo VMF grid file(s) found in goGPS directory. VMF gridded Zenith Delays will not be included in the outputfile.\nOtherwise you can provide the required VMF grid files & try again\n\n') 

         %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
         set(handles.cb_extract_VMF_ZTDs,'value',0)
        
      else %IF VMF grid FILES FOUND,LET USER CHOOSE GRID TYPE & RESOLUTION
           %I.E.['VMF1(2°x2.5°)','VMF3(1°x1°)','VMF3(5°x5°)']
         
           [VMF_type,grid_res] = VMFgrid_type(); 
           
           %SAVE GRID FILE VERSION & TYPE
           setappdata(0,'VMFgrid_type',VMF_type)
           setappdata(0,'VMFgrid_res',grid_res)
           
           %*********CHECK IF VMF_type,grid_res ARE BOTH EMPTY
           %NOTE:
           %     VMF_type,grid_res ARE ASSIGNED EMPTY WHEN USER CANCEL SELECTION
           %     OR CLOSE DIALOGUE BOX FOR THE SELECTION OF VMF GRID TYPE &
           %     RESOLUTION.IF IT HAPPENS SO, THE SELECTED TROPO DELAY MODEL,
           %     'VMF gridded ZTD' WOULD BE REPLACED WITH SAASTEMOIN MODEL
           %     OR THE POPUP MENU WILL BE SET TO SAASTAMOINEN AS DEFAULT
        
           if all([isempty(VMF_type),isempty(grid_res)])%IF USER CANCEL SELECTION/CLOSE DIALOGUE BOX
            
              beep %Give a beep sound
              errmsg702{1}=sprintf('No VMF grid file(s) Type & Resolution selected .\n');
              errmsg702{2}=sprintf('Try Again by selecting VMF grid file(s) Type & Resolution.\n');
              errmsg702{3}=sprintf('i.e. Select the VMF gridded ZTD model again and choose from the popup\n');
              warndlg(errmsg702,'VMF grid type Error','modal') 
              
              %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
              set(handles.cb_extract_VMF_ZTDs,'value',0)
              
             return
             
           end %// if all([isempty(VMF_type),isempty(grid_res)])
           
               
      end %//if VMF_grid_found == 0 
      
       
   end %//if all([~isempty(getappdata(0,'VMFgrids')),~isempty(getappdata(0,'VMF_grid_found'))])
      
end %//if all([getvmfZTDs_state == 1,strcmpi(getvmfZTDs_enable,'on')])

%SAVE CHECKBOX STATE(1 or 0)
setappdata(0,'getvmfZTDs_state',get(handles.cb_extract_VMF_ZTDs,'value'));

%==========================================================extract_VMF_ZTDs+
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%******PRECIPITABLE WATER VAPOUR(PWV) RETRIEVAL FROM TROPOPHERIC MODELS 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%1.GET & SAVE PWV CHECKBOX STATE(1 OR 0)
pwv_state = get(handles.cb_retrievePWV,'Value');%Get state of checkbox(1 or 0)

Tmpop_enable = get(handles.pop_Tm,'enable');%GET Tm POPUP MENU ENABILITY

if any([pwv_state == 1,strcmpi(Tmpop_enable,'on')])

%2.GET WEIGHTED MEAN TEMPERATURE(Tm) MODEL
%GET SOURCE OF MET PARAMETERS popup menu contents
getconts_Tm = cellstr(get(handles.pop_Tm,'String'));%get popup menu contents as cell array
Tm_model    = getconts_Tm{get(handles.pop_Tm,'Value')};%get selected item from pop_Tm
TmVAL       = get(handles.pop_Tm,'Value');%Get value of selected item from pop_Tm

%GET GRID RESOLUTION FOR GPT2w/GPT3 MODELS
grid_res1_pwv = get(handles.rb_grid_resolution_1_pwv,'value');%for 1° resolution
grid_res5_pwv = get(handles.rb_grid_resolution_5_pwv,'value');%for 5° resolution

%CHECK IF SELECTED Tm MODEL IS ANY OF THE GPT MODELS & GET GRID VALUES & RESOLUTION
if any([strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),...
        strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35),any(TmVAL==[11,12])])
    
   %IF 1° GRID RESOLUTION IS SELECTED 
   if all([grid_res1_pwv == 1,grid_res5_pwv == 0])
         
        %CHECK GRID VALUES & RESOLUTION STORED FROM THE Tm MODEL POPUP
         if all([any([~isempty(getappdata(0,'grid_Tm')) ~isempty(getappdata(0,'grid_res_Tm'))]),...
                 getappdata(0,'grid_res_Tm') == 1])  
             
            grid_TM     = getappdata(0,'grid_Tm'); %Grid Values
            grid_res_TM = 1; %Grid Resolution
         else
             grid_TM     = getappdata(0,'grid_Tm_1'); %Grid Values
             grid_res_TM = 1;%Grid Resolution
         end

   %IF 5° GRID RESOLUTION IS SELECTED 
   elseif all([grid_res1_pwv == 0,grid_res5_pwv == 1])
         
        %CHECK GRID VALUES & RESOLUTION STORED FROM THE Tm MODEL POPUP
         if all([any([~isempty(getappdata(0,'grid_Tm')) ~isempty(getappdata(0,'grid_res_Tm'))]),...
                 getappdata(0,'grid_res_Tm') == 5 ]) 
             
            grid_TM     = getappdata(0,'grid_Tm'); %Grid Values
            grid_res_TM = 5; %Grid Resolution
         else
             grid_TM     = getappdata(0,'grid_Tm_5'); %Grid Values
             grid_res_TM = 5;%Grid Resolution
         end
         
   end %//if all([grid_res1_pwv == 1,grid_res5_pwv == 0])
   
   %*****SAVE TM GRID VALUES & RESOLUTION 
  
   if all([~isempty(grid_TM),~isempty(grid_res_TM)])%IF THEY ARE NOT EMPTY([])
       
      setappdata(0,'grid_TM',grid_TM) %SAVE GRIDFILE FOR MET PARAMETER
      setappdata(0,'grid_res_TM',grid_res_TM)%SAVE GRID RESOLUTION FOR MET PARAMETER
      
   else %IF ANY/BOTH ARE EMPTY([]),READ NEW GRID VALUES USING "readGPTgrid" fxn
   
       %IF TM MODEL IS GPT2w
       if any([strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),TmVAL==11])
              
          %CHECK GRID REOLUTION TYPE(1° OR 5°)
          if all([grid_res1_pwv == 1,grid_res5_pwv == 0])
             grid_TM     = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRID VALUES
             grid_res_TM = 1;%GRID RESOLUTION 
                  
          elseif all([grid_res1_pwv == 0,grid_res5_pwv == 1])
                 grid_TM     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
                 grid_res_TM = 5;%GRID RESOLUTION
                     
          end 
            
             %IF TM MODEL IS GPT3
       elseif any([strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35),TmVAL==12])
              
              %CHECK GRID REOLUTION TYPE(1° OR 5°)
              if all([grid_res1_pwv == 1,grid_res5_pwv == 0])
                 grid_TM     = readGPTgrid('gpt3_1.mat','GPT3',1);%GRID VALUES
                 grid_res_TM = 1;%GRID RESOLUTION 
                  
              elseif all([grid_res1_pwv == 0,grid_res5_pwv == 1])
                     grid_TM     = readGPTgrid('gpt3_5.mat','GPT3',5);%GRID VALUES
                     grid_res_TM = 5;%GRID RESOLUTION
                     
              end  
              
       end  %//if if any([strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),TmVAL==11])
       
       setappdata(0,'grid_TM',grid_TM) %SAVE GRIDFILE FOR MET PARAMETER
       setappdata(0,'grid_res_TM',grid_res_TM)%SAVE GRID RESOLUTION FOR MET PARAMETER 
   
   end %//if all([~isempty(grid_TM),~isempty(grid_res_TM)])
   
   %****CHECK TIME VARIATION STATUS
   if ~isempty(getappdata(0,'Timevar_Tm'))
           
      %GET & SAVE TIME VARIATION STATUS
      setappdata(0,'Timevar_Tm',getappdata(0,'Timevar_Tm'))%time variation
        
   else   
       %SET DEFAULT VALUE OF 0 & SAVE 
       setappdata(0,'Timevar_Tm',0)%time variation
           
   end %if ~isempty(getappdata(0,'Timevar_Tm'))  
   
  switch getappdata(0,'Timevar_Tm')
      
      case 0
            Timevar_gpt = 'time variation (annual and semiannual terms)';
                         
      case  1
             Timevar_gpt = 'no time variation but static quantities';
  end                
    
elseif any([strncmpi(Tm_model,'GTVR-Tm model [Yao et al 2018]',30),TmVAL == 9]) 
    
       %GET & SAVE GTVR-Tm model GRID VALUES
        if ~isempty(getappdata(0,'gridV_GTVR'))
            
           setappdata(0,'gridV_GTVR',getappdata(0,'gridV_GTVR'))
           
        else    
            %SEARCH FOR GTVR-Tm model GRID FILE(GTVR-Tm.txt) & DIRECTORY & IMPORT VALUES               
            %Call the "searchGTVRgrid.m" function
            [~,gridV_GTVR] = SearchGTVRgrid();%****IMPORTED TVGG-Tm MODEL COEFFICIENTS
       
            %SAVE GRID VALUES
            setappdata(0,'gridV_GTVR',gridV_GTVR)
        end
       
       %GET & SAVE GPT2w model GRID VALUES
       if ~isempty(getappdata(0,'GPT2w_gridV'))
           
          setappdata(0,'GPT2w_gridV',getappdata(0,'GPT2w_gridV'))
          
       else    
           %READ GPT2w GRID FILE
           GPT2w_gridV     = readGPTgrid('gpt2_5w.mat','GPT2w',5);%GRID VALUES
       
          %SAVE GPT2w GRID VALUES
          setappdata(0,'GPT2w_gridV',GPT2w_gridV)
          
       end
    
elseif any([strncmpi(Tm_model,'TVGG-Tm model [Jiang et al 2018]',32),TmVAL == 10])
    
       %SET TOOLTIP STRING
       set(hObject,'TooltipString',strjoin({'PWV Retrieval Using','Time-Varying Global Gridded (TVGG) Tm model'}))%change Tooltip any time user select Tm Model type 
       
       %GET & SAVE TVGG-Tm model GRID VALUES
       if ~isempty(getappdata(0,'gridV_TVGG'))
       
          setappdata(0,'gridV_TVGG',getappdata(0,'gridV_TVGG'))
          
       else          
           %SEARCH FOR TVGG-Tm model GRID FILE(Coeff_TVGG_ERAI.mat) & DIRECTORY & LOAD VALUES               
           %Call the "SearchTVGGgrid.m" fxn
           [~,gridV_TVGG] = SearchTVGGgrid();%****LOADED TVGG-Tm MODEL COEFFICIENTS
       
           %SAVE GRID VALUES
           setappdata(0,'gridV_TVGG',gridV_TVGG)
           
       end
       
elseif any([strncmpi(Tm_model,'GTrop model [Sun et al 2019]',28),TmVAL == 15])
    
       %GET & SAVE GTrop GRID VALUES
        if ~isempty(getappdata(0,'gridV_GTrop'))
            
           setappdata(0,'gridV_GTrop',getappdata(0,'gridV_GTrop')) 
           
        else
            %SEARCH FOR GTrop model GRID FILE(GTropCoefficient.mat) & DIRECTORY & LOAD VALUES               
            %Call the "SearchGTropgrid.m" fxn
            [~,gridV_GTrop] = SearchGTropgrid();%****LOADED GTrop MODEL COEFFICIENTS
       
            %SAVE GRID VALUES
            setappdata(0,'gridV_GTrop',gridV_GTrop)
        end      
    
end %//any([strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),...
    %       strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35),any(TmVAL==[11,12])])

     
 %********GET & SAVE TM  Tm MODEL TYPE STRING FOR goGPS OUTPUT FILE
 %CHECK IF Tm MODEL TYPE STRING IS EMPTY([]) OR NOT
 if ~isempty(getappdata(0,'TmModel_type'))
     
    %SAVE TmModel_type
    setappdata(0,'TmModel_type',getappdata(0,'TmModel_type'))
    
 else
     if any([strncmpi(Tm_model,'Bevis et al (1992) model',24),TmVAL == 1])
 
        %CREATE OUTPUT STRING FOR Tm MODEL TYPE
        TmModel_type = 'Bevis et al (1992) [i.e. Tm = 70.2 + 0.72.*T]' ;
        
     elseif  any([strncmpi(Tm_model,'Bevis et al (1995) model',24),TmVAL == 2])
       
             %CREATE OUTPUT STRING FOR Tm MODEL TYPE
             TmModel_type = 'Bevis et al (1995) [i.e. Tm = 86.63 + 0.668.*T]' ;  
             
    elseif any([strncmpi(Tm_model,'Mendes et al (2000) model',25),TmVAL == 3])
       
           %CREATE OUTPUT STRING FOR Tm MODEL TYPE
           TmModel_type = 'Mendes et al (2000) [i.e. Tm = 50.4 + 0.789.*T]' ;
       
     elseif any([strncmpi(Tm_model,'Schueler et al (2001) model',27),TmVAL == 4])
       
            %CREATE OUTPUT STRING FOR Tm MODEL TYPE
            TmModel_type = 'Schueler et al (2001) [i.e. Tm = 86.9 + 0.647.*T]' ;
       
     elseif any([strncmpi(Tm_model,'Yao et al (2014) model',22),TmVAL == 5])
       
            %CREATE OUTPUT STRING FOR Tm MODEL TYPE
            TmModel_type = 'Yao et al (2014) [i.e. Tm = 43.69 + 0.8116.*T]' ;
       
     elseif any([strncmpi(Tm_model,'GTm-I model [Yao et al 2012]',28),TmVAL == 6])
       
            %CREATE OUTPUT STRING FOR Tm MODEL TYPE
            TmModel_type = 'Global Weighted Mean Temperature I (GWMT-I/GTm-I) [Yao et al 2012]';
       
     elseif any([strncmpi(Tm_model,'GTm-II model [Yao et al 2013]',29),TmVAL == 7])
       
            %CREATE OUTPUT STRING FOR Tm MODEL TYPE
            TmModel_type = 'Global Weighted Mean Temperature II (GWMT-II/GTm-II) [Yao et al 2013]';
       
     elseif any([strncmpi(Tm_model,'GTm-III model [Yao et al 2014]',30),TmVAL == 8])
       
            %CREATE OUTPUT STRING FOR Tm MODEL TYPE
            TmModel_type = 'Global Weighted Mean Temperature III (GWMT-III/GTm-III) [Yao et al 2014]';
       
     elseif any([strncmpi(Tm_model,'GTVR-Tm model [Yao et al 2018]',30),TmVAL == 9])
            
            %CREATE OUTPUT STRING FOR Tm MODEL TYPE
            TmModel_type = 'Global Weighted Mean Temperature Vertical Correction (GTVR-Tm) model [Yao et al 2018]';
    
     elseif any([strncmpi(Tm_model,'TVGG-Tm model [Jiang et al 2018]',32),TmVAL == 10])
            
            %CREATE OUTPUT STRING FOR Tm MODEL TYPE
            TmModel_type = 'Time-Varying Global Gridded (TVGG) model [Jiang et al 2018]';
            
            
     elseif any([strncmpi(Tm_model,'GPT2w model [Böhm et al 2014]',28),TmVAL == 11])
         
            %*************GRID RESOLUTION TYPE(1° OR 5°)
            %IF RESOLUTION IS 1°
            if grid_res1_pwv == 1 & grid_res5_pwv == 0
          
              %CREATE OUTPUT STRING FOR Tm MODEL TYPE
              TmModel_type =  strjoin({'GPT2w based on a 1° x 1° grid file with',Timevar_gpt,'[Böhm et al 2014]'});
         
            %IF RESOLUTION IS 5°
            elseif grid_res1_pwv == 0 & grid_res5_pwv == 1
              
                   %CREATE OUTPUT STRING FOR Tm MODEL TYPE
                   TmModel_type =  strjoin({'GPT2w based on a 5° x 5° grid file with',Timevar_gpt,'[Böhm et al 2014]'});
            end
            
     elseif any([strncmpi(Tm_model,'GPT3 model [Landskron & Böhm, 2018]',35),TmVAL == 12])
       
            %*************GRID RESOLUTION TYPE(1° OR 5°)
            %IF RESOLUTION IS 1°
           if grid_res1_pwv == 1 & grid_res5_pwv == 0
          
              %CREATE OUTPUT STRING FOR Tm MODEL TYPE
              TmModel_type =  strjoin({'GPT3 based on a 1° x 1° grid file with',Timevar_gpt,'[Landskron & Böhm, 2018]'});      
            
           %IF RESOLUTION IS 5°   
           elseif grid_res1_pwv == 0 & grid_res5_pwv == 1
              
                  %CREATE OUTPUT STRING FOR Tm MODEL TYPE
                  TmModel_type =  strjoin({'GPT3 based on a 5° x 5° grid file with',Timevar_gpt,'[Landskron & Böhm, 2018]'});      
           end
           
     elseif any([strncmpi(Tm_model,'UNB3m model [Leandro et al 2006]',32),TmVAL == 13])
       
            %CREATE OUTPUT STRING FOR Tm MODEL TYPE
            TmModel_type = 'University of New Brunswick v3 modified(UNB3m) [Leandro et al 2006]';
       
     elseif any([strncmpi(Tm_model,'Simplified PI  model [Manandhar et al 2017]',43),TmVAL == 14])
       
            %CREATE OUTPUT STRING FOR Tm MODEL TYPE
            TmModel_type = 'Simplified PI  model [Manandhar et al 2017]';
       
     elseif any([strncmpi(Tm_model,'GTrop model [Sun et al 2019]',28),TmVAL == 15])       
           
            %CREATE OUTPUT STRING FOR Tm MODEL TYPE
            TmModel_type = 'Global Tropospheric model(GTrop) [Sun et al 2019]';
            
     end %//if any([strncmpi(TmMODEL,'Bevis et al (1992) model',24),TmVAL == 1])
     
     %SAVE TmModel_type
     setappdata(0,'TmModel_type',TmModel_type) 
     
 end %//if ~isempty(getappdata(0,'TmModel_type'))
                       
%SAVE USER SELECTION OF Tm MODEL
setappdata(0,'Tm_model',Tm_model);%Save selected model as string
setappdata(0,'TmVAL',TmVAL);%Save selected model as value/numeric    

end %//if any([pwv_state == 1,strcmpi(Tmpop_enable,'on')])

%SAVE CHECKBOX STATE(1 or 0)
setappdata(0,'pwv_state',pwv_state);

%=============================================================PWV RETRIEVAL+
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%GET & ASVE STATE OF OPTION FOR OTHER REPORTS IN goGPS OUTPUTS
setappdata(0,'report_option',get(handles.cb_report_option,'value'))

%*********SAVE ATMOSPHERIC MEODELLING GUI COMPONENTS
%Call the "saveGUIState.m" fxn
saveGUIState(handles)
%                               -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+
%===============================END OF go_pb_Callback.m
%                               -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+   
 

%                          -------------------------------
%(2)***********************NON GUI COMPONENT SUB-ROUTINES
%                          -------------------------------

%                 ----------------------------------------------
%(1)**************FUNCTION TO IMPORT CSV FILE WITH MET PARAMTERS
%                 ----------------------------------------------
function [] = importMETPARA_csv()

%*********IMPORT FILE
%ALLOW USER TO SELECT csv FILE
[filename, pathname]=uigetfile({'*.xls;*.xlsx;*.csv;*.txt'},'Select File'); 
 
%***IF USER CHOOSE TO CANCEL IMPORT
if isequal (pathname,0) || isequal(filename, 0);%user pressed cancel  
     
   return
else    %if user does not cancel 
    
    %*******LET USER CHOOSE INPUT MET FORMAT
    %CALL THE 'metPARAformat.m' FUNCTION
    metPARAformat()
    
    %BLOCK PROGRAM EXECUTION AND WAIT TO RESUME WHEN USER CLICKS ON OK FROM
    %THE 'metPARAformat'
    %SUB-ROUTINE=======================================+
    uiwait(gcf)  
    
    %CONSTRUCT FULL PATH OF DATA
    fname= fullfile(pathname,filename);%FULL FILE PATH 
        
% end   

%****READ DATA FROM FILE
Data=importdata(fname);

%***GET NUMERIC PART OF DATA
 metpara=Data.data;
 
%*************DATA ARRANGEMENT / FORMAT
%DATA IS SHOULD BE IN 3 COLUMNS:
%COLUMN ARRANGEMENT IS AS FOLLOWS:
%1.[PRESSURE(P) TEMPERATURE(T) RELATIVE HUMIDITY (RH)]
%2.[PRESSURE(P) TEMPERATURE(T) WATER VAPOUR PRESSURE(WVP)(e)]
%3.[TEMPERATURE(T) PRESSURE(P) RELATIVE HUMIDITY (RH)]
%4.[TEMPERATURE(T) PRESSURE(P) WATER VAPOUR PRESSURE(WVP)(e)]
 
%RETRIEVE STORED MET FORMAT FROM THE 'metPARAformat' SUB-ROUTINE
METformat=getappdata(0,'metPARAformat');
 
%CHECK SELECTED FORMAT 
 
if strcmpi(METformat,'P,T,RH')
 %***GET SEPARATE DATA
 P=metpara(:,1);
 T=metpara(:,2);
 WV=metpara(:,3);%REPRESENTING EITHER RH or e WITH WV(WATER VAPOUR)
 
elseif strcmpi(METformat,'P,T,e')
       P=metpara(:,1);
       T=metpara(:,2);
       WV=metpara(:,3);
       
elseif strcmpi(METformat,'T,P,RH')
       T=metpara(:,1);
       P=metpara(:,2);
       WV=metpara(:,3);
    
elseif strcmpi(METformat,'T,P,e')
       T=metpara(:,1);
       P=metpara(:,2);
       WV=metpara(:,3);
       
end
%***REMOVING EMPTY CELLS(NAN) IF ANY
     
if any([isnan(T) isnan(P) isnan(WV)])%if contains nan
    
    T(isnan(T))=[];%deleting rows and columns with not a number(NAN)
    
    P(isnan(P))=[];%deleting rows and columns with not a number(NAN)
    
    WV(isnan(WV))=[];%deleting rows and columns with not a number(NAN)
end

%SAVE SELECTED MET PARAMETERS
 setappdata(0,'metPARA_csv',[P,T,WV])
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF importMETPARA_csv.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%              ------------------------------------------
%(2)***********SUB-ROUTINE TO SELECT IMPORTED DATA FORMAT
%              ------------------------------------------
function metPARAformat()

%FIGURE
S.fh = figure('units','pixels','position',[500 500 378 90], 'menubar','none','Color',[0.9,0.9,0.9],...
              'numbertitle','off','name','Import Format','resize','off');
          
  fc = get(S.fh,'color');
       
%UIPANEL        
S.p = uipanel('units','pix','pos',[15 40 349 45],'Title','Select Input Format','Fontsize',10,'Fontname'...
                ,'Cambria','Fontweight','bold','ForegroundColor','red','BackgroundColor',fc);                               
%PUSH BUTTON 
S.pb = uicontrol('style','push','units','pixels','position',[290 8 70 30],...
                 'fonts',14,'str','OK','callback',{@ok,S});
             
%RADIO BUTTONS         
rtext = {'P,T,RH';'P,T,e';'T,P,RH';'T,P,e'};     
rpos = {[30 45 75 20];[120 45 75 20];[190 45 75 20];[275 45 87 20]};
toolTipString = {'Select if format is [Pressure(P) , Temp(T) , Relative Humidity (RH) ]';...
                 'Select if format is [Pressure(P) , Temp(T) , Water Vapour Pressure (e) ]';...
                 'Select if format is [Temp(T) , Pressure(P) , Relative Humidity (RH)]';...
                 'Select if format is [Temp(T) , Pressure(P) , Water Vapour Pressure (e) ]'};         
for i=1:4
     
    S.rb(i)= uicontrol('style','radio','unit','pix','position',rpos{i},'string',rtext{i}...
                       ,'Fontsize',12,'Fontname','Cambria','Fontweight','bold','BackgroundColor',fc,...
                       'tooltip',toolTipString{i}); 
end 

%SET DEFAULT FORMAT
set((S.rb(1)),'value',1);%set radiobutton1 on
     
%SETTING CALLBACKS
set(S.rb(1),'callback',{@ptrh,S});
set(S.rb(2),'callback',{@pte,S});
set(S.rb(3),'callback',{@tprh,S});
set(S.rb(4),'callback',{@tpe,S});
 
%STORING RADIOBUTTONS
setappdata(0,'rb1',S.rb(1)) 
setappdata(0,'rb2',S.rb(2)) 
setappdata(0,'rb3',S.rb(3)) 
setappdata(0,'rb4',S.rb(4))                  
  
function[]=ptrh(varargin)
            
S = varargin{3};  % Get the structure.   
        
set(S.rb(1),'value',1);
set(S.rb(2),'value',0);
set(S.rb(3),'value',0);
set(S.rb(4),'value',0);

function[]=pte(varargin)
        
S = varargin{3};  % Get the structure.   
        
set(S.rb(1),'value',0);
set(S.rb(2),'value',1);
set(S.rb(3),'value',0);
set(S.rb(4),'value',0);

function[]=tprh(varargin)
            
S = varargin{3};  % Get the structure.   
        
set(S.rb(1),'value',0);
set(S.rb(2),'value',0);
set(S.rb(3),'value',1);
set(S.rb(4),'value',0);


function[]=tpe(varargin)
            
S = varargin{3};  % Get the structure.   
        
set(S.rb(1),'value',0);
set(S.rb(2),'value',0);
set(S.rb(3),'value',0);
set(S.rb(4),'value',1);
 
function[]=ok(varargin)
        
S = varargin{3};  % Get the structure

rb1=getappdata(0,'rb1');       
rb2=getappdata(0,'rb2');  
rb3=getappdata(0,'rb3');  
rb4=getappdata(0,'rb4');  
 
%CHECK WHICH RADIO BUTTON(FORMAT) IS SELECTED
if get(rb1,'value') == 1 
   format='P,T,RH';%Met format[Pressure(P),Temp(T),Relative Humidity(RH)] 
  
elseif get(rb2,'value') == 1
       format='P,T,e';%Met format [Pressure(T),Temp(T),Water Vapour Pressure(e)]
       
elseif get(rb3,'value') == 1
    format='T,P,RH';%Met format [Temp(T),Pressure(P),Relative Humidity(RH)]
    
elseif get(rb4,'value') == 1 
       format='T,P,e';%Met format [Temp(T),Pressure(P),Water Vapour Pressure(e)]
       
end

%SAVE MET FORMAT
setappdata(0,'metPARAformat',format);   
 
 close(S.fh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF metPARAformat.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%        ----------------------------------------------------------------
%(3)*****SUB-ROUTINE TO CREATE UITABLE FOR USER ENTRY/INPUT MET PARAMTERS
%        -----------------------------------------------------------------
function metpara_table()

%****FIGURE
S.fh = figure('units','pixels','Position', [1000 120 480 180],'menubar','none',...
              'numbertitle','off', 'name','goGPS : USER DEFINED METEOROLOGICAL PARAMETERS','resize','off');

fc=get(S.fh,'color');

%****UITABLE
S.t = uitable('Parent', S.fh, 'Position', [2 109 930 70]);

%COLUMN WIDTH
cw={157 157 162} ;

%COLUMN NAMEs
cn={'PRESSURE [ hpa ]','TEMPERATURE [°C ]','RH [ % ] OR  WVP(e) [ hpa ]'};


cData={[] [] []};%Initialize with empty([]) matrix/data
set(S.t,'ColumnName',cn,'ColumnWidth',cw,'RowName',[],'ColumnEditable',true,'SelectionHighlight','on','CellSelectionCallback',@cellSelect,'data',cData,'Fontsize',12) ;     

%UIPANEL        
S.p = uipanel('units','pix','pos',[3 55 270 45],'Title','Select Input Format','Fontsize',10,'Fontname'...
                ,'Cambria','Fontweight','bold','ForegroundColor','red','BackgroundColor',fc);                               


%RADIO BUTTONS         
rtext = {'P,T,RH';'P,T,e'};     
rpos = {[30 60 75 20];[160 60 75 20]};
toolTipString = {'Select if format is [Pressure(P) , Temp(T) , Relative Humidity (RH) ]';...
                 'Select if format is [Pressure(P) , Temp(T) , Water Vapour Pressure (e) ]'};         
for i=1:2
     
    S.rb(i)= uicontrol('style','radio','unit','pix','position',rpos{i},'string',rtext{i}...
                       ,'Fontsize',12,'Fontname','Cambria','Fontweight','bold','BackgroundColor',fc,...
                       'tooltip',toolTipString{i}); 
end

%SET DEFAULT FORMAT
set((S.rb(1)),'value',1);%set radiobutton1 on  

%SETTING CALLBACKS
set(S.rb(1),'callback',{@PTRH,S});
set(S.rb(2),'callback',{@PTe,S});

%STORING RADIOBUTTONS
setappdata(0,'rb1',S.rb(1)) 
setappdata(0,'rb2',S.rb(2)) 
        
              
%PUSH BUTTON FOR ADDING NEW ROW
S.pb1 = uicontrol('style','pushbutton','units','pix','position',[5 10 100 30],...
                  'string','Add Row','ForegroundColor','black','Fontsize',12,'Fontweight','bold',...
                  'tooltip','Click to Add Row','callback',{@AddRow,S});
              
%PUSH BUTTON FOR DELETING SELECTED ROW
S.pb2 = uicontrol('style','pushbutton','units','pix','position',[130 10 100 30],...
                  'string','Delete Row','ForegroundColor','black','Fontsize',12,'Fontweight','bold',...
                  'tooltip','Click to Delete Row','callback',{@DeleteRow,S});
                          
%PUSH BUTTON TO SAVING INPUTs
S.pb3 = uicontrol('style','pushbutton','units','pix','position',[255 10 100 30],...
                  'string','OK','ForegroundColor','BLUE','Fontsize',12,'Fontweight','bold',...
                  'tooltip','Click to Save inputs','callback',{@OK,S});
                                   
%PUSH BUTTON TO EXIT FIGURE            
S.pb4 = uicontrol('style','pushbutton','units','pix','position',[377 10 100 30],...                                
                  'string','Exit','ForegroundColor','red','Fontsize',12,'Fontweight','bold',...
                  'tooltip','Click to Exit','callback',{@exit,S});                       
%--------------------------------------------------------------------------
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-              
              
%*************WORKING ON VARIOUS UI CONTROLS
%             ------------------------------
%1.RADIO BUTTONS
function[]=PTRH(varargin)
            
S = varargin{3};  % Get the structure.   
        
set(S.rb(1),'value',1);
set(S.rb(2),'value',0);

function[]=PTe(varargin)
        
S = varargin{3};  % Get the structure.   
        
set(S.rb(1),'value',0);
set(S.rb(2),'value',1);

%2.PUSH BUTTONS              
%2.1.****ADD NEW ROW              
function[]=AddRow(varargin) 
%GET STRUCTUREs
S = varargin{3};      
              
%GET UITABLE 
uitable=S.t;

%GET DATA FROM TABLE
tableData=get(uitable, 'data'); % Gets the actual data for use

%ADD ONE NEW ROW
if iscell(tableData)
   tableData(end+1,:)={[]};% Add new row of []'s as a cell type
 else
     tableData(end+1,:)=[];
 end

set(uitable, 'data', tableData); %Update uitable with a rows with string type and other rows as cell type

function cellSelect(src,evt)
%GET INDICES OF SELECTED CELLs AND MAKE THEM AVAILABLE FOR OTHER CALLBACKS    
currentCell=evt.Indices;%GET CURRENT CELL  FROM IUTABLE

if any(currentCell) %Loop necessary to surpress unimportant errors
    rows=currentCell(:,1);
    set(src,'UserData',rows)    
end

%2.2****DELETE SELECTED ROW              
function[]=DeleteRow(varargin) 
%GET STRUCTUREs
S = varargin{3};      
              
%GET UITABLE 
uitable=S.t;

%GET DATA FROM TABLE
tableData=get(uitable, 'data'); % Gets the actual data for use

if isempty(tableData)
             
   tableData={''};%CREATE EMPTY CELL ARRAY       
end
%GET INDICES OF SELECTED ROWS
index=get(uitable,'UserData');

%CREATE MASK CONTAINING ROWS TO KEEP
mask=(1:size(tableData,1))';

mask(index)=[];

%DELETE SELECTED ROWS AND RE-WRITE DATA
tableData=tableData(mask,:);
set(uitable,'data',tableData)


%2.3.SAVE INPUTDATA
function[]=OK(varargin) 
 
%GET STRUCTUREs
S = varargin{3};  

%GET UITABLE 
uitable = S.t;

%GET DATA FROM TABLE
data_str = get(uitable, 'data'); % Gets the actual data in string
data = str2double(data_str); % Change all strings inside the uitable to a double type
 
%CHECK IF TABLE IS EMPTY([])
if all(isnan(data))
   beep %Give a beep sound
   errmsg{1}=sprintf('No input for meteorological parameters.\n');
   errmsg{2}='Please Provide meteorological parameters & Try Again.';
   warndlg(errmsg,'Meteorological Input Error','modal')
   return 
   
elseif any(isnan(data))
       beep %Give a beep sound
       errmsg{1}=sprintf('Missing input(s) for meteorological parameters.\n');
       errmsg{2}='Please Provide all meteorological parameters & Try Again.';
       warndlg(errmsg,'Meteorological Input Error','modal')
       return 
end
      
%GET RADIO BUTTONS
rb1 = getappdata(0,'rb1'); % FOR 'P,T,RH  format     
rb2 = getappdata(0,'rb2'); % FOR 'P,T,e  format    

%CHECK WHICH RADIO BUTTON(FORMAT) IS SELECTED
if get(rb1,'value') == 1 
   format = 'P,T,RH';%Met format[Pressure(P),Temp(T),Relative Humidity(RH)] 
  
elseif get(rb2,'value') == 1
       format = 'P,T,e';%Met format [Pressure(T),Temp(T),Water Vapour Pressure(e)]       
end

%*************DATA ARRANGEMENT / FORMAT
%DATA IS SHOULD BE IN 3 COLUMNS:
%COLUMN ARRANGEMENT IS AS FOLLOWS:
%1.[PRESSURE(P) TEMPERATURE(T) RELATIVE HUMIDITY (RH)]
%2.[PRESSURE(P) TEMPERATURE(T) WATER VAPOUR PRESSURE(WVP)(e)]

 %***GET SEPARATE DATA
 P  = data(:,1);%Pressure
 T  = data(:,2);%Temperature
 WV = data(:,3);%REPRESENTING EITHER RH or e WITH WV(WATER VAPOUR)
  
%***REMOVING EMPTY CELLS(NAN) IF ANY    
if any([isnan(T) isnan(P) isnan(WV)])%if contains nan
    
    T(isnan(T))   = [];%deleting rows and columns with not a number(NAN)
    
    P(isnan(P))   = [];%deleting rows and columns with not a number(NAN)
    
    WV(isnan(WV)) = [];%deleting rows and columns with not a number(NAN)
end

%SAVE MET DATA
setappdata(0,'metPARA_table',[P,T,WV])

%SAVE MET FORMAT
setappdata(0,'metPARAformat',format);   

%CLOSE FIGURE
close(S.fh)
              
%3.EXIT PROGRAMME
 function [] = exit(varargin)
        
S = varargin{3};  % Get the structure.
        
 close(S.fh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF metpara_table.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%              -----------------------------------------------------
%(4)***********SUB-ROUTINE FOR TIME VARIATION APPLICATION(GPT MODELS)
%              ------------------------------------------------------
function TIMEvariation()

%FIGURE
S.fh = figure('units','pixels','position',[1220 500 250 90], 'menubar','none','Color',[0.9,0.9,0.9],...
              'numbertitle','off','name','Time Variation Status','resize','off');
          
  fc = get(S.fh,'color');
       
%UIPANEL       
S.p = uipanel('units','pix','pos',[15 40 225 45],'Title','Apply Time Variation or Static','Fontsize',10,'Fontname'...
              ,'Cambria','Fontweight','bold','ForegroundColor','red','BackgroundColor',fc);                               
%PUSH BUTTON[290 8 70 30] 
S.pb = uicontrol('style','push','units','pixels','position',[170 8 70 30],...
                 'fonts',14,'str','OK','callback',{@Ok,S});
             
%CHECK BOX         
rtext = {'Apply Time Variation'};     
rpos = {[30 45 180 20]};
toolTipString = {'Select to apply time variation ( i.e. annual and semiannual terms )';};         
for i=1
     
    S.cb(i)= uicontrol('style','checkbox','unit','pix','position',rpos{i},'string',rtext{i}...
                       ,'Fontsize',12,'Fontname','Cambria','Fontweight','bold','BackgroundColor',fc,...
                       'Value',0,'TooltipString',toolTipString{i}); 
end 

%SETTING CALLBACKS
set(S.cb,'callback',{@timeVAR,S});

%STORING CHECKBOX
setappdata(0,'cb',S.cb) 

function[]=timeVAR(varargin)
%NOTE:
%    GLOBAL PRESSURE & TEMPERATURE(GPT) MODELS SUCH AS GPT2,GPT2w & GPT3 
%    COMPUTES METEOROLOGICAL PARAMETERS EITHER IN STATIC MODE OR IN TIME
%    VARYING MODE:
%    case 1: no time variation but static quantities
%    case 0: with time variation (annual and semiannual terms)

S = varargin{3};  % Get the structure.   

%GET CHECK BOX STATE(0,or 1)
Timevar = get(S.cb,'value');

%******TOGGLE CHECKBOX
if get(S.cb,'value') == 1 
   Timevar = Timevar - 1; %set Timevar to zero
   set(S.cb,'TooltipString','UnCheck to apply static terms');
else
     Timevar = Timevar + 1;
     set(S.cb,'TooltipString','Select to apply time variation ( i.e. annual and semiannual terms )');
end 
%SAVE Timevar
setappdata(0,'Timevar',Timevar)

function[]=Ok(varargin)
        
S = varargin{3}; % Get the structure

%GET CHECK BOX HANDLE
 cb=getappdata(0,'cb');
 
%GET CHECK BOX STATE (0, OR 1) 
 Timevar = get(cb,'value');
 
 %GET CHECK BOX STATE (0, OR 1) 
 if get(cb,'value') == 1 
    Timevar = Timevar - 1; %set Timevar to zero
    
 else 
     Timevar = Timevar + 1;
     
 end
 
%SAVE Timevar(TIME VARIATION)
setappdata(0,'Timevar',Timevar);   
 
 close(S.fh)
%==========================END OF TIMEvariation.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%              --------------------------------------------------------
%(5)***********SUB-ROUTINE FOR CHECKING THE EXISTENCE OF VMF GRID FILES
%              --------------------------------------------------------
function [VMF_grid_found,VMFgrids,VMF_type,grid_res] = SearchVMFgrids(handles)

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
%      [VMF_grid_found,VMFgrids] = SearchVMFgrids()                      +
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
%INPUT                                                                     *
%1.handles : goGPS GUI FIG HANDLES. NOTE :INPUT IS OPTIONAL                *

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
%9.    Orography: Orography file                                           *                     *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
% WRITTEN BY: 
%            OSAH SAMUEL, MSC GEOMATIC ENGINEERING (PhD STUDENT)           +
%            Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com                +
%            Phone:+233(0)246137410 / +233(0)509438484                     +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%**************************************************************************+
if nargin == 0
   handles = [];
end


%GET TROPOSPHERIC MODELING GUI COMPONENTS
if ~isempty(handles)
                  
   %GET TROPOSPHERIC MODELLING GUI COMPONENTS FROM goGPS GUI FIG
   %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
   %***************TROPOSPHERIC MODELLING
   %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
              
   %GET SELECTED MODEL TYPE [OPTION FOR TROPO MODELLING]
   tropmodel_c=get(handles.rb_Tropo_combine,'Value');%Get combine tropmodel rb [option button for combine tropo models]
   tropmodel_s=get(handles.rb_Tropo_separate,'Value');%Get separate tropmodel rb[option button for separate tropo models]
               
   %============================================================+
   %***GET SELECTED COMBINE(DRY+WET) TROPO MODEL
   %------------------------------------------------------------+
   getTROPOconts = cellstr(get(handles.pop_lTropo_combine,'String'));%get Tropo_combine popup menu contents as cell array
   tropoModel=getTROPOconts{get(handles.pop_lTropo_combine,'Value')};%get selected item from pop_lTropo_combine
   Modelv=get(handles.pop_lTropo_combine,'Value');%get selected item value from pop_lTropo_combine
              
   %HYDROSTATIC model
   getcontsH = cellstr(get(handles.pop_ITropo_hydrostatic,'String'));%get Tropo_hydrostatic popup menu contents as cell array
   dryModel=getcontsH{get(handles.pop_ITropo_hydrostatic,'Value')};%get selected item from pop_Tropo_hydrostatic
   Modelv_H=get(handles.pop_ITropo_hydrostatic,'Value');%get selected item value from pop_ITropo_hydrostatic
       
   %WET model
   getcontsW = cellstr(get(handles.pop_ITropo_wet,'String'));%get Tropo_wet popup menu contents as cell array
   wetModel=getcontsW{get(handles.pop_ITropo_wet,'Value')};%get selected item from pop_Tropo_wet
   Modelv_W=get(handles.pop_ITropo_wet,'Value');%get selected item value from pop_ITropo_wet
              
   %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
   %****************GET SELECTED MAPPING FUNCTION(MF)
   %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
       
   %HYDROSTATIC MF
   getcontsHMF = cellstr(get(handles.pop_IHydrostic_MF,'String'));%get Hydrostatic MF popup menu contents as cell array
   MFh_model=getcontsHMF{get(handles.pop_IHydrostic_MF,'Value')};%get selected item from pop_IHydrostic_MF
          
   %WET MF
   getcontsWMF = cellstr(get(handles.pop_IWet_MF,'String'));%get Wet MF popup menu contents as cell array
   MFw_model=getcontsWMF{get(handles.pop_IWet_MF,'Value')};%get selected item from pop_IWet_MF
    
   %GET SELECTED METEOROLOGICAL PARAMETER SOURCE
   getconts_MET = cellstr(get(handles.pop_source_metpara,'String'));%get MET popup menu contents as cell array
   METpara = getconts_MET{get(handles.pop_source_metpara,'Value')};%get selected item from MET POPUP MENU
   metVAL=get(handles.pop_source_metpara,'Value');%GET SELECTED VALUE
   
end

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
          
             %DISPLAY ERROR/WARNING MESSAGE BASED ON SOME PARTICULAR SELECTED TROPO MODELS
             if ~isempty(handles)  
                  
              if all([tropmodel_c == 1,tropmodel_s == 0])
                   
                 if any([Modelv == 17,strncmpi(tropoModel,'VMF gridded ZTD',15)])
              
                    beep %Give a beep sound
                    errmsg1{1}=sprintf('No VMF grid file(s) found in goGPS directory.\n');
                    errmsg1{2}=sprintf('Please Provide VMF grid file(s) & Try Again.\n');
                    errmsg1{3}=sprintf('VMF grid(s) can as well be downloaded from : http://vmf.geo.tuwien.ac.at/trop_products/GRID/\n');
                    warndlg(errmsg1,'VMF grid file Error','modal')
                    
                    %SET COMBINE TROPO MODEL POPUP MENU VALUE TO 2(I.E.SAASTAMOINEN)
                    set(handles.pop_lTropo_combine,'Value',2)
                    
                    %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
                    set(handles.cb_extract_VMF_ZTDs,'value',0,'enable','on')
                    
                   %******SET MET PARAMETER SOURCE UICONTROLS VISIBLE/ENABLE ON
                    set(handles.pop_source_metpara,'enable','on')
                    set(handles.text_source_metpara,'enable','on')
    
                    if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
                            any(metVAL==[7,8])])
        
                       if strcmpi(get(handles.pop_met_manual,'visible'),'Off') 
                          set(handles.pop_met_manual,'visible','On')
                          set(handles.text_metmanual_source,'visible','On')
                       end 
     
                    elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
                                any(metVAL==[4,5])])
          
                           if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'off'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'off')])
                              set(handles.rb_grid_resolution_1,'Visible','on')
                              set(handles.rb_grid_resolution_5,'Visible','on')
                              set(handles.text_grid_resolution,'Visible','on')  
                           end 
          
                    end %//if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
                        %          any(metVAL==[7,8])]) 
                    
                 end %//if any([Modelv == 17,strncmpi(tropoModel,'VMF ZTD',7)]) 
                  
              elseif all([tropmodel_c == 0,tropmodel_s == 1])
                   
                     if any([any([any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv_H == 15]),any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv_W == 20])]),...
                             all([any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv_H == 15]),any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv_W == 20])])])
                          
                        beep %Give a beep sound
                        errmsg2{1}=sprintf('No VMF grid file(s) found in goGPS directory.\n');
                        errmsg2{2}=sprintf('Please Provide VMF grid file(s) & Try Again.\n');
                        errmsg2{3}=sprintf('VMF grid(s) can as well be downloaded from : http://vmf.geo.tuwien.ac.at/trop_products/GRID/\n');
                        warndlg(errmsg2,'VMF grid file Error','modal')
                        
                        if any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv_H == 15])
                            
                           %SET DRY TROPO MODEL POPUP MENU VALUE TO 2(I.E.SAASTAMOINEN)
                           set(handles.pop_ITropo_hydrostatic,'Value',1)
                           
                           %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
                           set(handles.cb_extract_VMF_ZTDs,'value',0,'enable','on')
                           
                        end %//if any([strncmpi(dryModel,'VMF gridded ZHD',15),Modelv_H == 15])
                        
                        
                        if any([strncmpi(wetModel,'VMF gridded ZWD',15),Modelv_W == 20])
                            
                           %SET WET TROPO MODEL POPUP MENU VALUE TO 1(I.E.SAASTAMOINEN)
                           set(handles.pop_ITropo_wet,'Value',1)
                           
                           %SET extract_VMF_ZTDs CHECKBOX VALUE TO 0,enable on
                           set(handles.cb_extract_VMF_ZTDs,'value',0,'enable','on')

                        end
                        
                        %******SET MET PARAMETER SOURCE UICONTROLS VISIBLE/ENABLE ON
                        set(handles.pop_source_metpara,'enable','on')
                        set(handles.text_source_metpara,'enable','on')
    
                        if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
                                        any(metVAL==[7,8])])
        
                           if strcmpi(get(handles.pop_met_manual,'visible'),'Off') 
                              set(handles.pop_met_manual,'visible','On')
                              set(handles.text_metmanual_source,'visible','On')
                           end  
     
                        elseif any([strcmpi(METpara,'4.Generate from GPT2w model'),strcmpi(METpara,'5.Generate from GPT3 model'),...
                                            any(metVAL==[4,5])])
          
                               if any([strcmpi(get(handles.rb_grid_resolution_1,'Visible'),'off'),strcmpi(get(handles.rb_grid_resolution_5,'Visible'),'off')])
                                  set(handles.rb_grid_resolution_1,'Visible','on')
                                  set(handles.rb_grid_resolution_5,'Visible','on')
                                  set(handles.text_grid_resolution,'Visible','on')  
                               end  
          
                        end  %//if any([strcmpi(METpara,'7.Enter Meteorological Parameters'),strcmpi(METpara,'8.Import from csv file[T P RH]'),...
                               %                  any(metVAL==[7,8])])
    
                     end
                     

                     if get(handles.rb_different_MF,'Value') == 1 %if Use other Mapping Models is selected
                      
                        %GET & SAVE GRIDFILE IF USER SELECTS ANY OF THE HYDROSTATIC VMF MODELs 
                        if any([any([any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),...
                                          strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12)]),...
                                     any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),...
                                          strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12)])]),...
                                all([any([strncmpi(MFh_model,'VMF1(1° x 1°)',12),strncmpi(MFh_model,'VMF1(5° x 5°)',12),...
                                          strncmpi(MFh_model,'VMF3(1° x 1°)',12),strncmpi(MFh_model,'VMF3(5° x 5°)',12)]),...
                                     any([strncmpi(MFw_model,'VMF1(1° x 1°)',12),strncmpi(MFw_model,'VMF1(5° x 5°)',12),...
                                          strncmpi(MFw_model,'VMF3(1° x 1°)',12),strncmpi(MFw_model,'VMF3(5° x 5°)',12)])])])        
             
                           %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR           
                           setappdata(0,'VMF_grid_found_mf',VMF_grid_found)%SAVED FOR MF MODEL
                           
                           errmsg3{1}=sprintf('No VMF grid file(s) found in goGPS directory.\n');
                           errmsg3{2}=sprintf('System will resort to GPT2w / GPT3 grid file(s) instead.\n');
                           errmsg3{3}=sprintf('Mean while, VMF grid(s) can as well be downloaded from : http://vmf.geo.tuwien.ac.at/trop_products/GRID/\n');
                           warndlg(errmsg3,'VMF grid file Error','modal')
                           
                        end
                        
                     end %//if get(handles.rb_different_MF,'Value') == 1
                     
              end %//if all([tropmodel_c == 1,tropmodel_s == 0])
             
               
              %MAPPING FUCTION FOR PPP TROPO MODELLING
              if all([get(handles.flag_tropo,'value') == 1, strcmpi(get(handles.flag_tropo,'enable'),'on')])
                  
                 msg1{1}=sprintf('No VMF grid file(s) found in goGPS directory.\n');
                 msg1{2}=sprintf('System will resort to GPT2w / GPT3 grid file(s) instead.\n');
                 msg1{3}=sprintf('Mean while, VMF grid(s) can as well be downloaded from : http://vmf.geo.tuwien.ac.at/trop_products/GRID/\n');
                 warndlg(errmsg3,'VMF grid file Error','modal')
                 
              end
              

              else
                  beep %Give a beep sound
                  errmsg4{1}=sprintf('No VMF grid file(s) found in goGPS directory.\n');
                  errmsg4{2}=sprintf('Please Provide VMF grid file(s) & Try Again.\n');
                  errmsg4{3}=sprintf('VMF grid(s) can as well be downloaded from : http://vmf.geo.tuwien.ac.at/trop_products/GRID/\n');
                  warndlg(errmsg4,'VMF grid file Error','modal')
                  
                  return
             end
              
           end
            
         end 
         
      end 
      
   end
   
end
                          
% SAVE VMF grid FILES & GRIDFILE FOUND/NOT FOUND INDICATOR
setappdata(0,'VMF_grid_found',VMF_grid_found)
setappdata(0,'VMFgrids',VMFgrids)

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
        
      if all([~exist('h0Files','var'),~exist('h6Files','var'),~exist('h12Files','var'),~exist('h18Files','var'),...
              ~exist('h0_1Files','var'),~exist('h6_1Files','var'),~exist('h12_1Files','var'),~exist('h18_1Files','var')]) 
         
         %SET VMF GRID FILE & OROGRAPHY FILE NAMES & PATHS TO EMPTY([])  
         h0files  = []; h6files  = []; h12files = []; h18files = [];
         H0files  = []; H6files  = []; H12files = []; H18files = [];
         orogFILE = [];
         orogPATH = [];
         
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
   
else
     orofile = []; %SET orography FILE TO EMPTY([]) IF IT DOES NOT EXIST IN FOLDER   
     oropath = [];%SET orography FILE PATH TO EMPTY([])
    
   
end      

%***********CHECK IF ALL GRID FILES ARE EMPTY([]) 
if any([all([isempty(H0files),isempty(H6files),isempty(H12files),isempty(H18files),isempty(orofile)]),...
        all([isempty(H0files),isempty(H6files),isempty(H12files),isempty(H18files),~isempty(orofile)])])
   
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
%                           END OF  SearchVMFgrids.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


function [VMF_type,grid_res] = VMFgrid_type()
    
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *                                              
%            "VMFgrid_type" subroutine allows user to determine/Select the +   
%            Vienna Mapping Function(VMF) Grid File Version & Type         +    
    
choice = questdlg('SELECT VMF GRID TYPE :','Vienna Mapping Function(VMF) Grid Files Type',...
                  'VMF1(2°x2.5°) Grid','VMF3(1°x1°) Grid','VMF3(5°x5°) Grid','VMF3(5°x5°) Grid');                
               
%************GET USER CHOICE OF SELECTION 
         
if ischar(choice)
   userResponse = {'VMF1(2°x2.5°) Grid','VMF3(1°x1°) Grid','VMF3(5°x5°) Grid'};
   
else   
     userResponse = {1,2,3};
    
end   

%HANDLE RESPONSE 
switch choice
          
    case userResponse 
         
         if any([strcmpi(choice,'VMF1(2°x2.5°) Grid'),strfind(choice,'VMF1(2°x2.5°)')])
             
            VMF_type = 'VMF1'; 
            
            grid_res = 2;
            
             %IF VMF3 FILE TYPE, SPECIFY GRID RESOLUTION
         else  
              VMF_type = 'VMF3';
             
              if any([strcmpi(choice,'VMF3(1°x1°) Grid'),strfind(choice,'VMF3(1°x1°)')])
          
                 grid_res = 1 ;
           
              elseif any([strcmpi(choice,'VMF3(5°x5°) Grid'),strfind(choice,'VMF3(5°x5°)')])
          
                     grid_res = 5;     
              end    
              
         end   
        
    otherwise %WHEN USER CANCEL SELECTION OR CLOSE DIALOGUE BOX FOR THE 
              %SELECTION OF VMF GRID TYPE & RESOLUTION.
          
              %ASSIGN EMPTY([]) MATRIX
              VMF_type = [];
              grid_res = [];
                 
              return
                
end %//switch choice
  
%=========================================END OF VMFgrid_type.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-        


%=====FOR READING SOME WEIGHTED MEAN TEMPERATURE(Tm) MODEL'S GRID FILES====
%     -----------------------------------------------------------------

%              -----------------------------------------------------------
%(6)***********SUB-ROUTINE FOR CHECKING THE EXISTENCE OF GTVR-Tm GRID FILE
%              ------------------------------------------------------------
function [GTVR_grid,gridV_GTVR,folder] = SearchGTVRgrid()
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *                                              
%            The "SearchGTVRgrid" subroutine searches Global Weighted...   +
%            Mean Temperature Vertical Correction model(GTVR-Tm)Grid file  +
%            in a given Directory/Folder as indicated in the sub-routine.  + 
%            Situation where GTVR-Tm file is not found in the default ...  +   
%            directory/folder,it searches recursively through all ...      +  
%            Sub-Directories/Folders of the given Directory.GTVR-Tm file... + 
%            is extracted by looping through all the listed files in the...+
%            provided folder or sub-folders.Finally,if GTVR-Tm file is ... + 
%            still not found in the default directory/folder and its ...   +
%            sub-folders,the search for GTVR-Tm file is extended to the ...+
%            current directory.
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%%USAGE:                                                                   +
%       [GTVR_grid,folder] = SearchGTVRgrid()                              +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%NOTE:                                                                     *
%     THE FF FOLDERS/DIRECTORIES ARE SEARCHED FOR GTVR GRID FILE           *
%1.   'GTVR-Tm file' folder; main folder for GTVR GRIDs                 *
%2.   'TropoGRIDS' folder; main folder for Tropospheric & Tm GRIDs         *
%3.   'data' folder; folder for goGPS data                                 *
%4.   'pwd'; Current folder/Directory. In goGPS, it is the 'goGPS' folder  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%--------------------------------------------------------------------------
%OUTPUTs                                                                   *
%1.   GTVR_grid : Flag to indicate presence of GTVR grid file.1 to indicate*
%                 there is file in the folder & 0 to mean absence of file  *
%2.  gridV_GTVR :  Extracted GTVR grid Values                              *
%3.     folder  : The found directory or folder(File Path) for GTVR grid...* 
%                 file.e.g.:'C:...data\TropoGRIDS\GTVR-Tm file\GTVR-Tm.txt *       
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
% WRITTEN BY: 
%            OSAH SAMUEL, MSC GEOMATIC ENGINEERING (PhD STUDENT)           +
%            Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com                +
%            Phone:+233(0)246137410 / +233(0)509438484                     +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%**************************************************************************+

%         -----------------------------------------------------
%*********SEARCH FOR GTVR GRID FILES FROM THE VARIOUS DIRECTORY
%         -----------------------------------------------------

%1.SEARCH GTVR GRIDs from 'GTVR-Tm' folder
[GTVR_grid,GTVRgDir] = SearchGTVRGrid('../data/TropoGRIDS/GTVR-Tm file');

if isequal(GTVR_grid,0)%if GTVR grids are not in 'GTVR-Tm file' folder, try the 'TropoGRIDS'
    
   %2.SEARCH GTVR-Tm GRIDs from 'TropoGRIDS' folder 
   [GTVR_grid,GTVRgDir] = SearchGTVRGrid('../data/TropoGRIDS');

   if isequal(GTVR_grid,0)%if GTVR grids are not in 'TropoGRIDS' folder, try the 'data' folder
    
      %3.SEARCH GTVR-Tm GRIDs from 'data' folder 
      [GTVR_grid,GTVRgDir] = SearchGTVRGrid('../data');%searching the 'data' folder
   
      if isequal(GTVR_grid,0) %if GTVR grids are not in the 'data' folder, try the 'current' directory
       
         %4.SEARCH GTVR-Tm GRIDs from 'current' folder 
         [GTVR_grid,GTVRgDir] = SearchGTVRGrid(pwd);
      
         if isequal(GTVR_grid,0)%if GTVR grids are not in the 'current' directory, then issue a message
          
            beep %Give a beep sound
            errmsg2{1}=sprintf('No GTVR grid file found in goGPS directory.\n');
            errmsg2{2}=sprintf('Please Provide GTVR grid file & Try Again.\n');
            warndlg(errmsg2,'GTVR grid file Error','modal')
            
            %ASSIGN EMPTY([]) MATRIX
            gridV_GTVR = [];
            folder     = [];
             
            return
      
         else  
             folder = GTVRgDir;%GTVR grid directory/full path [EX.C:\Users\...data\TropoGRIDS\GTVR-Tm file\GTVR-Tm.txt]
      
         end  
      
      else  
          folder = GTVRgDir;%GTVR grid directory/full path 
       
      end  
   
   else  
       folder = GTVRgDir;%GTVR grid directory/full path 
    
   end 
   
else
    folder = GTVRgDir;%GTVR grid directory/full path 
    
end

%            ---------------------------------------------
%************IF A DIRECTORY/FOLDER WITH GTVR GRID IS FOUND....
%            ---------------------------------------------
 if ~isempty(folder) %If Folder is not Empty
    
    %**IMPORT GTVR GRID file & ASSIGN TO gridV_GTVR 
    gridV_GTVR = importdata(folder);%GRID VALUES 
    
 end
 
  %SAVE GTVR GRID DIRECTORY & VALUES
  setappdata(0,'gridV_GTVR',gridV_GTVR) 
  setappdata(0,'GTVRgDir',folder) 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 

%            -----------------------------------------------------------
%***********MAIN SUB-ROUTINE FOR CHECKING THE EXISTENCE OF GTVR GRID FILE
%            ------------------------------------------------------------

function [GTVR_grid,Folder] = SearchGTVRGrid(directory)

%ASSIGNMENT
folder = directory;  %folder/Directory suppose to be containing GTVR FILES in goGPS

%****FIRST CHECK IF FOLDER IS EMPTY OR NOT
SizeListing=size(dir(folder)); %Size of folder content

%********CHECK IF FOLDER IS EMPTY OR NOT
if any([SizeListing(1,1)==2,SizeListing(1,1)==0])% size of 2 means folder is empty & 0 means folder doesn't exist
                  
   GTVR_grid = 0; %Flag to indicate the absence of GTVR grid file in directory
   Folder = []; %File directory is empty([]) if GTVR-Tm.txt is not found
   
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
          
       %*****GET INDEXEs FOR THE VARIOUS GTVR GRID FILES(i.e.Integer Val /[])
       GTVRIndex = regexp(fileList,'\w*GTVR-Tm.txt'); %matches any words with GTVR-Tm.txt.  
          
       %************GET THE GRID FILES
       GTVRFiles  = fileList(~cellfun(@isempty,GTVRIndex));%GTVR files(eg:'GTVR-Tm.txt')
       
    else
        GTVRFiles = [];
               
    end  %if ~isempty(fileList)
    
    %IF NO GTVR FILES EXIST IN THE CURRENT FOLDER, CHECK SUB-FOLDERS/              
    %SUB-DIRECTORIES FOR GTVR GRID FILEs
    %                              +                        ==========
     if isempty(fileList)
              
        checkSubfolders = 1; %Flag to search Sub-Directories/folders
           
     elseif isempty(GTVRFiles)
                   
            checkSubfolders = 1;%Flag to search Sub-Directories/folders 
                         
     else
          if ~isempty(GTVRFiles)
          
             checkSubfolders=0;%Flag not to search Sub-Directories/folders
             
             GTVR_grid = 1; %Flag to indicate the presence of GTVR grid file in directory
             
             %***CONSTRUCT THE FULL PATH [File directory [if GTVR-Tm.txt file is found)]      
             Folder = cellfun(@(x) fullfile(folder,x),GTVRFiles,'UniformOutput',false);%Prepend path to files  
             
             %SORT FILES & REMOVE DUPLICATE FILES
             [GTVRFiles,i_GTVR] = unique(GTVRFiles);
             Folder =(Folder(i_GTVR));
             
             %CONVERT CELL ARRAY TO CHARACTER
             if iscell(Folder)
                Folder = char(Folder);
             end    
             
             if iscell(GTVRFiles)
                GTVRFiles = char(GTVRFiles);
             end    
             
          end
           
     end  %\\if isempty(fileList)
                
     
     if checkSubfolders==1 %if this condition is true,then either files in the current
                           %Directory were not GTVR grid file or there were no files at all.        
                           %System has to look through Sub-Folders or Sub-directories for GTVR Grid File                                                              
                                                                                                                                                                                                                                                                                  
        %***OPEN ALL SUB-FOLDERS & GET GTVR GRID FILES
                
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
            
                     j=1;%Loop index for the VMF grid file types
                     
                      for i = 1:length(filelists)%****LOOP OVER THE FILES & CHECK FOR GTVR  GRID FILE
                                                 %Identify the needed file and extract it separately
                                                 %into different file name(GTVRfiles)
                                                 %Note that dir also lists the directories,so you have to check for them. 
                                                         
                           %***GET NAME OF FILES FROM THE LIST OF FILES FOR EACH LOOP               
                           
                           FileName = filelists{i,1};
                           filepath = filepaths{i,1};

                           %********MATCH ANY WORD STARTING GTVR-Tm.txt
                           if regexpi(FileName,'\w*GTVR-Tm.txt') 

                             %When Tested & the file is a GTVR file
                             GTVRfiles{j,1}=FileName; %#ok<*NASGU>
                             GTVRpath{j,1}=filepath;%GTVR file file path in cells
                             
                             j=j+1;%Update index
                                                                           
                           end  %//if regexpi(FileName,'\w*GTVR-Tm.txt')
                           
                      end %//for i = 1:length(fileLists) 
                                                                                                               
                  end   %//if ~isempty(fileLists)
                      
            else                          
                 %WRITE ERROR MESSAGE
                 erMsg{q,1}=sprintf('The Folder %s and  It''s Sub-Folder, %s , are empty / Contain No GTVR Grid file .\n',folder,subDirs{iDir});
                         
                 q=q+1; %update index
                        
            end   %if Sizesublist(1,1)~=2
                   
        end    %for iDir = find(validIndex)
        
      if ~exist('GTVRfiles','var')
         
         GTVR_grid = 0;%Flag to indicate the absence of GTVR grid file in directory
         
         Folder = []; %File directory is empty([]) if GTVR-Tm.txt file is not found
         
         return
     
      else 
          GTVR_grid = 1; %Flag to indicate the presence of GTVR grid file in directory
          
          Folder = GTVRpath; %if GTVR-Tm.txt file is found
          
          %SORT FILES & REMOVE DUPLICATE FILES
          [GTVRFiles,i_GTVR] = unique(GTVRfiles);
          Folder =(Folder(i_GTVR));
             
             %CONVERT CELL ARRAY TO CHARACTER
             if iscell(Folder)
                Folder = char(Folder);
             end    
             
             if iscell(GTVRFiles)
                GTVRFiles = char(GTVRFiles);
             end    
                    
      end %//~exist('GTVRfiles','var')                                        
           
     end  %//if checkSubfolders==1
     
end %//if exist('listFiles','var')

%=========================================END OF SearchGTVRGrid.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%                        END OF SearchGTVRgrid
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%              -----------------------------------------------------------
%(7)***********SUB-ROUTINE FOR CHECKING THE EXISTENCE OF TVGG-Tm GRID FILE
%              ------------------------------------------------------------
function [TVGG_grid,gridV_TVGG,folder] = SearchTVGGgrid()
%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *                                              
%            The "SearchTVGGgrid" subroutine searches Time-Varying Global..+
%            Gridded Mean Temperature(TVGG-Tm)model Grid file in a given...+
%            Directory/Folder as indicated in the sub-routine. Situation...+ 
%            where TVGG-Tm file is not found in the default directory/...  +   
%            folder,it searches recursively through all Sub-Directories ...+  
%            /Folders of the given Directory.TVGG-Tm file is extracted ... + 
%            by looping through all the listed files in the provided folder+
%            or sub-folders.Finally,if TVGG-Tm file is still not found in  + 
%            the default directory/folder and its sub-folders,the search...+
%            for TVGG-Tm file is extended to the current directory.        +           
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%%USAGE:                                                                   +
%       [TVGG_grid,gridV_TVGG,folder] = SearchTVGGgrid()                   +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%NOTE:                                                                     *
%     THE FF FOLDERS/DIRECTORIES ARE SEARCHED FOR TVGG GRID FILE           *
%1.   'TVGG-Tm file' folder; main folder for TVGG GRID                     *
%2.   'TropoGRIDS' folder; main folder for Tropospheric & Tm GRIDs         *
%3.   'data' folder; folder for goGPS data                                 *
%4.   'pwd'; Current folder/Directory. In goGPS, it is the 'goGPS' folder  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%--------------------------------------------------------------------------
%OUTPUTs                                                                   *
%1.   TVGG_grid : Flag to indicate presence of TVGG grid file.1 to indicate*
%                 there is file in the folder & 0 to mean absence of file  *
%2.  gridV_TVGG :  Extracted TVGG grid Values                              *
%3.     folder  : The found directory or folder(File Path) for TVGG grid...* 
%                 file.e.g.:                                               *
%                 'C:...data\TropoGRIDS\TVGG-Tm file\Coeff_TVGG_ERAI.mat   *       
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
% WRITTEN BY: 
%            OSAH SAMUEL, MSC GEOMATIC ENGINEERING (PhD STUDENT)           +
%            Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com                +
%            Phone:+233(0)246137410 / +233(0)509438484                     +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%**************************************************************************+

%1.SEARCH TVGG-Tm GRIDs from 'TVGG-Tm' folder
[TVGG_grid,TVGGgDir] = SearchTVGGGrid('../data/TropoGRIDS/TVGG-Tm file');

if isequal(TVGG_grid,0)%if TVGG grids are not in 'TVGG-Tm file' folder, try the 'TropoGRIDS'
    
   %2.SEARCH TVGG-Tm GRIDs from 'TropoGRIDS' folder 
   [TVGG_grid,TVGGgDir] = SearchTVGGGrid('../data/TropoGRIDS');

   if isequal(TVGG_grid,0)%if TVGG grids are not in 'TropoGRIDS' folder, try the 'data' folder
    
      [TVGG_grid,TVGGgDir] = SearchTVGGGrid('../data');%searching the 'data' folder
   
      if isequal(TVGG_grid,0) %if TVGG grids are not in the 'data' folder, try the 'current' directory
       
         [TVGG_grid,TVGGgDir] = SearchTVGGGrid(pwd);
      
         if isequal(TVGG_grid,0)%if TVGG grids are not in the 'current' directory, then issue a message
          
            beep %Give a beep sound
            errmsg2{1}=sprintf('No TVGG grid file found in goGPS directory.\n');
            errmsg2{2}=sprintf('Please Provide TVGG grid file & Try Again.\n');
            warndlg(errmsg2,'TVGG grid file Error','modal') 
            
            %ASSIGN EMPTY([]) MATRIX
            gridV_TVGG = [];
            folder     = [];
            
            return
            
         else  
             folder = TVGGgDir;%TVGG grid directory/full path [EX.C:\Users\...data\TropoGRIDS\TVGG-Tm file\Coeff_TVGG_ERAI.mat]
      
         end  
      
      else  
          folder = TVGGgDir; %TVGG grid directory/full path
       
      end  
   
   else   
       folder = TVGGgDir; %TVGG grid directory/full path
       
   end 
       
else 
     folder = TVGGgDir; %TVGG grid directory/full path  
    
end  

%            -------------------------------------------------
%************IF A DIRECTORY/FOLDER WITH TVGG GRIDs IS FOUND....
%            -------------------------------------------------
 if ~isempty(folder) %If Folder is not Empty
    
    %****LOAD TVGG-Tm MODEL COEFFICIENTS
    GridV_TVGG = load(folder);%512*256 ERA-INTERIM GRIDS
    gridV_TVGG = GridV_TVGG.Coeff_TVGG_ERAI;
         
 end
 
 %SAVE DIRECTORY & GRID VALUES 
 setappdata(0,'TVGGgDir',folder) 
 setappdata(0,'gridV_TVGG',gridV_TVGG) 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+ 
 
%            -----------------------------------------------------------
%***********MAIN SUB-ROUTINE FOR CHECKING THE EXISTENCE OF TVGG GRID FILE
%            ------------------------------------------------------------

function [TVGG_grid,Folder] = SearchTVGGGrid(directory)

%ASSIGNMENT
folder = directory;  %folder/Directory suppose to be containing TVGG FILE in goGPS

%****FIRST CHECK IF FOLDER IS EMPTY OR NOT
SizeListing=size(dir(folder)); %Size of folder content

%********CHECK IF FOLDER IS EMPTY OR NOT
if any([SizeListing(1,1)==2,SizeListing(1,1)==0])% size of 2 means folder is empty & 0 means folder doesn't exist
                  
   TVGG_grid = 0; %Flag to indicate the absence of TVGG grid file in directory
   Folder = []; %File directory is empty([]) if Coeff_TVGG_ERAI.mat is not found
   
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
          
       %*****GET INDEXEs FOR THE TVGG GRID FILE(i.e.Integer Val /[])
       TVGGIndex = regexp(fileList,'\w*Coeff_TVGG_ERAI.mat'); %matches any words with Coeff_TVGG_ERAI.mat.  
          
       %************GET THE GRID FILES
       TVGGFiles  = fileList(~cellfun(@isempty,TVGGIndex));%TVGG files(eg:'Coeff_TVGG_ERAI.mat')
       
    else
        TVGGFiles = [];
           
    end  %if ~isempty(fileList)
    
    %IF NO TVGG FILES EXIST IN THE CURRENT FOLDER, CHECK SUB-FOLDERS/              
    %SUB-DIRECTORIES FOR TVGG GRID FILEs
    %                              +                        ==========
     if isempty(fileList)
              
        checkSubfolders = 1; %Flag to search Sub-Directories/folders
           
     elseif isempty(TVGGFiles)
                   
            checkSubfolders = 1;%Flag to search Sub-Directories/folders 
                         
     else
          if ~isempty(TVGGFiles)
          
             checkSubfolders=0;%Flag not to search Sub-Directories/folders
             
             TVGG_grid = 1; %Flag to indicate the presence of TVGG grid file in directory
             
             %***CONSTRUCT THE FULL PATH [File directory [if Coeff_TVGG_ERAI.mat file is found)]       
             Folder = cellfun(@(x) fullfile(folder,x),TVGGFiles,'UniformOutput',false);%Prepend path to files [File directory [if GTVR-Tm.txt file is found)] 
             
             %SORT FILES & REMOVE DUPLICATE FILES
             [TVGGFiles,i_TVGG] = unique(TVGGFiles);%GET TVGG FILES with no repetitions & INDEX(i_TVGG)
             Folder =(Folder(i_TVGG));%FILES WITH PATH
             
             %CONVERT CELL ARRAY TO CHARACTER
             if iscell(Folder)
                Folder=char(Folder);
             end    
             
             if iscell(TVGGFiles)
                TVGGFiles = char(TVGGFiles);
             end    
             
          end
           
     end  %\\if isempty(fileList)
                
     
     if checkSubfolders==1 %if this condition is true,then either files in the current
                           %Directory were not TVGG grid file or there were no files at all.        
                           %System has to look through Sub-Folders or Sub-directories for TVGG Grid File                                                              
                                                                                                                                                                                                                                                                                  
        %***OPEN ALL SUB-FOLDERS & GET TVGG GRID FILE
                
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
            
                     j=1;%Loop index for the TVGG grid file type
                     
                      for i = 1:length(filelists)%****LOOP OVER THE FILES & CHECK FOR TVGG GRID FILE
                                                 %Identify the needed files and extract them separately
                                                 %into different file name(TVGGfiles)
                                                 %Note that dir also lists the directories,so you have to check for them. 
                                                         
                           %***GET NAME OF FILES FROM THE LIST OF FILES FOR EACH LOOP               
                           
                           FileName = filelists{i,1};
                           filepath = filepaths{i,1};

                           %********MATCH ANY WORD WITH Coeff_TVGG_ERAI.mat
                           if regexpi(FileName,'\w*Coeff_TVGG_ERAI.mat') 

                             %When Tested & the file is a TVGG file
                             TVGGfiles{j,1}=FileName; %#ok<*NASGU>
                             TVGGpath{j,1}=filepath;%TVGG file Path in cells
                             
                             j=j+1;%Update index
                                                                           
                           end  %//if regexpi(FileName,'\w*Coeff_TVGG_ERAI.mat') 
                           
                      end %//for i = 1:length(fileLists) 
                                                                                                               
                  end   %//if ~isempty(fileLists)
                      
            else                          
                 %WRITE ERROR MESSAGE
                 erMsg{q,1}=sprintf('The Folder %s and  It''s Sub-Folder, %s , are empty / Contain No TVGG Grid file .\n',folder,subDirs{iDir});
                         
                 q=q+1; %update index
                        
            end   %if Sizesublist(1,1)~=2
                   
        end    %for iDir = find(validIndex)
        
      if ~exist('TVGGfiles','var')
         
         TVGG_grid = 0;%Flag to indicate the absence of TVGG grid file in directory
         
         Folder = []; %File directory is empty([]) if Coeff_TVGG_ERAI.mat file is not found
         
         return
     
      else 
          TVGG_grid = 1; %Flag to indicate the presence of TVGG grid file in directory
          
          Folder = TVGGpath; %if Coeff_TVGG_ERAI.mat file is found
          
          %SORT FILES & REMOVE DUPLICATE FILES
          [TVGGFiles,i_TVGG] = unique(TVGGfiles);%GET TVGG FILES with no repetitions & INDEX(i_TVGG)
           Folder =(Folder(i_TVGG));%FILES WITH PATH
         
          %CONVERT CELL ARRAY TO CHARACTER
          if iscell(Folder)
             Folder=char(Folder);
          end 
          
          if iscell(TVGGFiles)
             TVGGFiles = char(TVGGFiles);
          end     
                    
      end %//~exist('TVGGfiles','var')                                        
           
     end  %//if checkSubfolders==1
     
end %//if exist('listFiles','var')

%=========================================END OF SearchTVGGGrid.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%                        END OF SearchTVGGgrid
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

%              -----------------------------------------------------------
%(8)***********SUB-ROUTINE FOR CHECKING THE EXISTENCE OF GTrop GRID FILE
%              ------------------------------------------------------------
function [GTrop_grid,gridV_GTrop,folder] =SearchGTropgrid()
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
 
%           -----------------------------------------------------------
%***********MAIN SUB-ROUTINE FOR CHECKING THE EXISTENCE OF GTrop GRID FILE
%           ------------------------------------------------------------

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


%              -----------------------------------------------------------
%(9)*****************SUB-ROUTINE TO SELECT PPP MAPPING FUNCTION
%              ------------------------------------------------------------

function [] = MF_ppp()

%CREATE A FIGURE[300 300 300 100]
S.fh = figure('units','pixels','position',[380 525 300 100],'menubar','none',...
              'name','Select Mapping Function','numbertitle','off','resize','off');

%**************ADD UICONTROLS
%1.POPUP MENU
S.pp = uicontrol('style','pop','unit','pix','position',[10 60 280 20],'fontsize',10,...
                 'fontweight','normal','string',{'NMF [Niell Mapping Function]';'GMF [Global Mapping Function]';...
                 'VMF grided / GPT series models (VMF1 or VMF3)'},...
                 'value',2,'callback',@MF);

%2.PUSH BUTTON(OK)
S.pb = uicontrol('style','push','unit','pix','position',[200 5 80 30],'fontsize',12,...
                 'fontweight','bold','string','Ok','callback',@PPPok);             
             
%SAVE FIG             
guidata(S.fh,S)

%MOVE GUI TO A SPECIFIED POSITION
% movegui('center')

%*****CALLBACK FUNCTION for POPUP MENU
function [] = MF(varargin)
%GET GUI/FIG DATA
S = guidata(gcbf);

%GET POPUP CONTENTs('string)
POPconts = get(S.pp,'string');%get popup menu contents as cell array
POPval = get(S.pp, 'Value');%get selected item value/Index (eg: 1 for first item )
mf = POPconts{POPval}; %Selected mapping funtion type

%SAVE MF INDEX FROM POPUP MENU
setappdata(0,'PPPmf_POPval',POPval)

if POPval == 3
    
   %GET SAVED goGPS GUI FIGURE HANDLES
   handles = getappdata(0,'handle'); 
   
   %CHECK IF VMF GRID FILEs ARE AVAILABLE IN goGPS
   %Call for the 'SearchVMFgrids.m' fxn
   [VMF_grid_found,VMFgrids] = SearchVMFgrids(handles); 
                        
   %*********SAVE STRUCT GRID FILES
   setappdata(0,'VMFgrids_ppp',VMFgrids)
         
   %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
   setappdata(0,'VMF_grid_found_ppp',VMF_grid_found)
   
   
   if VMF_grid_found == 1 %IF VMF grid FILES FOUND,LET USER CHOOSE GRID TYPE & RESOLUTION
                               %I.E.['VMF1(2°x2.5°)','VMF3(1°x1°)','VMF3(5°x5°)']
                               
      [VMF_type,grid_res] = VMFgrid_type(); 
      
      %SAVE GRID FILE VERSION & TYPE
      setappdata(0,'VMFgrid_type_ppp',VMF_type)
      setappdata(0,'VMFgrid_res_ppp',grid_res)
      
      grid_MF_ppp = [];%ASSIGN EMPTY([]) MATRIX TO GPT series grid values
      gridRES_MF_ppp = [];
      GPTmodel_ppp = []; %GPT MODEL WITH EMPTY([]) MATRIX
      
      %*********CHECK IF VMF_type,grid_res ARE BOTH EMPTY
      %NOTE:
      %     VMF_type,grid_res ARE ASSIGNED EMPTY WHEN USER CANCEL SELECTION
      %     OR CLOSE DIALOGUE BOX FOR THE SELECTION OF VMF GRID TYPE &
      %     RESOLUTION.IF IT HAPPENS SO, THE SELECTED TROPO DELAY MODEL,
      %     'VMF gridded ZTD' WOULD BE REPLACED WITH SAASTEMOIN MODEL
      %     OR THE POPUP MENU WILL BE SET TO SAASTAMOINEN AS DEFAULT
      
      if all([isempty(VMF_type),isempty(grid_res)])%IF USER CANCEL SELECTION/CLOSE DIALOGUE BOX
          
         beep %Give a beep sound
         warnmsg1{1}=sprintf('VMF grid type & resolution has not been selected .\n');
         warnmsg1{2}=sprintf('Try Again and select VMF grid type & resolution.\n');        
         warndlg(warnmsg1,'VMF grid file Error','modal')
         
         %READ FROM GPT3 GRID FILE(using GPT3 1° x 1° model)
         grid_MF_ppp=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDVALUES  
         gridRES_MF_ppp = 1;%GRID RESOLUTION
         
         MFh_ppp = 'VMF3(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
         MFw_ppp = 'VMF3(1° x 1°)';%WET MAPPING FUNCTION MODEL 
         GPTmodel_ppp = 'GPT3'; %GPT3 MODEL
         
         if isempty(grid_MF_ppp) %if GPT3 model grid is unavailable
            
            %READ FROM GPT2w GRID FILE(using GPT2w 1° x 1° model)
            grid_MF_ppp = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
            gridRES_MF_ppp = 1;%GRID RESOLUTION
         
            MFh_ppp = 'VMF1(1° x 1°)'; %HYDROSTATIC MAPPING FUNCTION MODEL 
            MFw_ppp = 'VMF1(1° x 1°)'; %WET MAPPING FUNCTION MODEL
            GPTmodel_ppp = 'GPT2w'; %GPT2w MODEL
        
         end        
            
      else %IF USER DOES NOT CANCEL SELECTION/CLOSE DIALOGUE BOX
           if any([strcmpi(VMF_type,'VMF1'),strfind(VMF_type,'VMF1')])
               
              MFh_ppp = 'VMF1(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
              MFw_ppp = 'VMF1(1° x 1°)';%WET MAPPING FUNCTION MODEL
              
           elseif any([strcmpi(VMF_type,'VMF3'),strfind(VMF_type,'VMF3')])
                  
                  if grid_res == 1
                      
                     MFh_ppp = 'VMF3(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
                     MFw_ppp = 'VMF3(1° x 1°)';%WET MAPPING FUNCTION MODEL
                     
                  elseif grid_res == 5
                      
                       MFh_ppp = 'VMF3(5° x 5°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
                       MFw_ppp = 'VMF3(5° x 5°)';%WET MAPPING FUNCTION MODEL
                       
                  end 
                      
           end %//if any([strcmpi(VMF_type,'VMF1'),strfind(VMF_type,'VMF1')])
           
      end %//if all([isempty(VMF_type),isempty(grid_res)])
      
   else %IF VMF grid FILES NOT FOUND
        
        %SAVE GRID FILE VERSION & TYPE WITH EMPTY MATRIX([])
        setappdata(0,'VMFgrid_type_ppp',[])
        setappdata(0,'VMFgrid_res_ppp',[])
           
        beep %Give a beep sound
        errmsg3{1}=sprintf('No VMF grid file(s) found in goGPS directory.\n');
        errmsg3{2}=sprintf('System will resort to GPT2w / GPT3 grid file(s) instead.\n');
        errmsg3{3}=sprintf('Mean while, VMF grid(s) can as well be downloaded from : http://vmf.geo.tuwien.ac.at/trop_products/GRID/\n');
        warndlg(errmsg3,'VMF grid file Error','modal')
        
        %READ FROM GPT3 GRID FILE(using GPT3 1° x 1° model)
        grid_MF_ppp=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE 
        gridRES_MF_ppp = 1;%GRID RESOLUTION
        
        MFh_ppp = 'VMF3(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
        MFw_ppp = 'VMF3(1° x 1°)';%WET MAPPING FUNCTION MODEL 
        GPTmodel_ppp = 'GPT3'; %GPT3 MODEL
        
        if isempty(grid_MF_ppp) %if GPT3 model grid is unavailable
            
           %READ FROM GPT2w GRID FILE(using GPT2w 1° x 1° model)
           grid_MF_ppp = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
         
           MFh_ppp = 'VMF1(1° x 1°)'; %HYDROSTATIC MAPPING FUNCTION MODEL 
           MFw_ppp = 'VMF1(1° x 1°)'; %WET MAPPING FUNCTION MODEL 
           GPTmodel_ppp = 'GPT2w'; %GPT2w MODEL
           
        end              
   
   end %//if VMF_grid_found == 1
   
    %SAVED GPT2w/GPT3 GRID VALUES & RESOLUTION
    %HYDROSTATIC
    setappdata(0,'grid_MFh_ppp',grid_MF_ppp)
    setappdata(0,'gridRES_MFh_ppp',gridRES_MF_ppp)
      
    %WET
    setappdata(0,'grid_MFw_ppp',grid_MF_ppp)
    setappdata(0,'gridRES_MFw_ppp',gridRES_MF_ppp)
    
    %GPT MODEL
    setappdata(0,'GPTmodel_ppp',GPTmodel_ppp)
     
elseif POPval == 2
       MFh_ppp = 'GMF';
       MFw_ppp = 'GMF';
       
    
elseif POPval == 1
       MFh_ppp = 'NMF';
       MFw_ppp = 'NMF';   
    
 end   
    %SAVE MAPPING FUNCTION MODELS
    setappdata(0,'MFh_model_ppp',MFh_ppp)
    setappdata(0,'MFw_model_ppp',MFw_ppp)
     
    
%*****CALLBACK FUNCTION for pushbutton
function [] = PPPok(varargin)

%GET GUI/FIG DATA
S = guidata(gcbf);

%GET POPUP CONTENTs('string)
POPconts = get(S.pp,'string');%get popup menu contents as cell array
POPval = get(S.pp, 'Value');%get selected item value/Index (eg: 1 for first item )
mf = POPconts{POPval}; %Selected mapping funtion type

%SAVE MAPPING FUNCTION MODELS
if all([~isempty(getappdata(0,'MFh_model_ppp')),~isempty(getappdata(0,'MFw_model_ppp'))])
    
   %COMPARE SELECTED MF MODEL FROM POPUP MENU UNDER  MF SUB-ROUTINE & THAT
   %RETRIEVE UNDER PPPok SUB-ROUTINE 
   if getappdata(0,'PPPmf_POPval') == POPval
    
      if all([POPval == 1,any([strcmpi(getappdata(0,'MFh_model_ppp'),'NMF'),strcmpi(getappdata(0,'MFw_model_ppp'),'NMF')])])  
 
         setappdata(0,'MFh_model_ppp',getappdata(0,'MFh_model_ppp'))
         setappdata(0,'MFw_model_ppp',getappdata(0,'MFw_model_ppp'))
         
      elseif all([POPval == 2,any([strcmpi(getappdata(0,'MFh_model_ppp'),'GMF'),strcmpi(getappdata(0,'MFw_model_ppp'),'GMF')])])  
 
             setappdata(0,'MFh_model_ppp',getappdata(0,'MFh_model_ppp'))
             setappdata(0,'MFw_model_ppp',getappdata(0,'MFw_model_ppp'))
             
      elseif  all([POPval == 3,any([any([strcmpi(getappdata(0,'MFh_model_ppp'),'VMF1(1° x 1°)'),strcmpi(getappdata(0,'MFw_model_ppp'),'VMF1(1° x 1°)')]),...
                                    any([strcmpi(getappdata(0,'MFh_model_ppp'),'VMF3(1° x 1°)'),strcmpi(getappdata(0,'MFw_model_ppp'),'VMF3(1° x 1°)')]),...
                                    any([strcmpi(getappdata(0,'MFh_model_ppp'),'VMF1(5° x 5°)'),strcmpi(getappdata(0,'MFw_model_ppp'),'VMF1(5° x 5°)')])])])  
 
              setappdata(0,'MFh_model_ppp',getappdata(0,'MFh_model_ppp'))
              setappdata(0,'MFw_model_ppp',getappdata(0,'MFw_model_ppp'))
              
              %GET & SAVE STRUCT GRID FILES
              setappdata(0,'VMFgrids_ppp',getappdata(0,'VMFgrids_ppp'))
         
              %GET & SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
              setappdata(0,'VMF_grid_found_ppp',getappdata(0,'VMF_grid_found_ppp')) 
        
              %SAVE GRID FILE VERSION & TYPE WITH EMPTY MATRIX([])
              setappdata(0,'VMFgrid_type_ppp',getappdata(0,'VMFgrid_type_ppp'))
              setappdata(0,'VMFgrid_res_ppp',getappdata(0,'VMFgrid_res_ppp')) 
       
              %GET & SAVED GPT2w/GPT3 GRID VALUES & RESOLUTION
              %HYDROSTATIC
              setappdata(0,'grid_MFh_ppp',getappdata(0,'grid_MFh_ppp'))
              setappdata(0,'gridRES_MFh_ppp',getappdata(0,'gridRES_MFh_ppp'))
      
              %WET
              setappdata(0,'grid_MFw_ppp',getappdata(0,'grid_MFw_ppp'))
              setappdata(0,'gridRES_MFw_ppp',getappdata(0,'gridRES_MFw_ppp'))
              
              %GPT MODEL
              setappdata(0,'GPTmodel_ppp',getappdata(0,'GPTmodel_ppp'))
             
      else
          
          if POPval == 3
        
             %GET SAVED goGPS GUI FIGURE HANDLES
             handles = getappdata(0,'handle');  
         
             %CHECK IF VMF GRID FILEs ARE AVAILABLE IN goGPS
             %Call for the 'SearchVMFgrids.m' fxn
             [VMF_grid_found,VMFgrids] = SearchVMFgrids(handles); 
                        
             %*********SAVE STRUCT GRID FILES
             setappdata(0,'VMFgrids_ppp',VMFgrids)
         
             %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
             setappdata(0,'VMF_grid_found_ppp',VMF_grid_found)
   
   
             if VMF_grid_found == 1 %IF VMF grid FILES FOUND,LET USER CHOOSE GRID TYPE & RESOLUTION
                               %I.E.['VMF1(2°x2.5°)','VMF3(1°x1°)','VMF3(5°x5°)']
                               
                [VMF_type,grid_res] = VMFgrid_type(); 
      
                %SAVE GRID FILE VERSION & TYPE
                setappdata(0,'VMFgrid_type_ppp',VMF_type)
                setappdata(0,'VMFgrid_res_ppp',grid_res)
      
                grid_MF_ppp = [];%ASSIGN EMPTY([]) MATRIX TO GPT series grid values
                gridRES_MF_ppp = [];
                GPTmodel_ppp = []; %GPT MODEL WITH EMPTY([]) MATRIX 
      
                %*********CHECK IF VMF_type,grid_res ARE BOTH EMPTY
                %NOTE:
                %     VMF_type,grid_res ARE ASSIGNED EMPTY WHEN USER CANCEL SELECTION
                %     OR CLOSE DIALOGUE BOX FOR THE SELECTION OF VMF GRID TYPE &
                %     RESOLUTION.IF IT HAPPENS SO, THE SELECTED TROPO DELAY MODEL,
                %     'VMF gridded ZTD' WOULD BE REPLACED WITH SAASTEMOIN MODEL
                %     OR THE POPUP MENU WILL BE SET TO SAASTAMOINEN AS DEFAULT
      
                if all([isempty(VMF_type),isempty(grid_res)])%IF USER CANCEL SELECTION/CLOSE DIALOGUE BOX
          
                   beep %Give a beep sound
                   warnmsg1{1}=sprintf('VMF grid type & resolution has not been selected .\n');
                   warnmsg1{2}=sprintf('Try Again and select VMF grid type & resolution.\n');        
                   warndlg(warnmsg1,'VMF grid file Error','modal')
         
                  %READ FROM GPT3 GRID FILE(using GPT3 1° x 1° model)
                  grid_MF_ppp=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDVALUES  
                  gridRES_MF_ppp = 1;%GRID RESOLUTION
         
                  MFh_ppp = 'VMF3(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
                  MFw_ppp = 'VMF3(1° x 1°)';%WET MAPPING FUNCTION MODEL 
                  GPTmodel_ppp = 'GPT3'; %GPT3 MODEL 
                  
                  if isempty(grid_MF_ppp) %if GPT3 model grid is unavailable
            
                     %READ FROM GPT2w GRID FILE(using GPT2w 1° x 1° model)
                     grid_MF_ppp = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
                     gridRES_MF_ppp = 1;%GRID RESOLUTION
         
                     MFh_ppp = 'VMF1(1° x 1°)'; %HYDROSTATIC MAPPING FUNCTION MODEL 
                     MFw_ppp = 'VMF1(1° x 1°)'; %WET MAPPING FUNCTION MODEL 
                     GPTmodel_ppp = 'GPT2w'; %GPT2w MODEL 
                     
                  end          
            
                else %IF USER DOES NOT CANCEL SELECTION/CLOSE DIALOGUE BOX
                    if any([strcmpi(VMF_type,'VMF1'),strfind(VMF_type,'VMF1')])
               
                       MFh_ppp = 'VMF1(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
                       MFw_ppp = 'VMF1(1° x 1°)';%WET MAPPING FUNCTION MODEL
              
                    elseif any([strcmpi(VMF_type,'VMF3'),strfind(VMF_type,'VMF3')])
                  
                           if grid_res == 1
                      
                              MFh_ppp = 'VMF3(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
                              MFw_ppp = 'VMF3(1° x 1°)';%WET MAPPING FUNCTION MODEL
                     
                           elseif grid_res == 5
                      
                                  MFh_ppp = 'VMF3(5° x 5°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
                                  MFw_ppp = 'VMF3(5° x 5°)';%WET MAPPING FUNCTION MODEL
                       
                           end   
                      
                    end %//if any([strcmpi(VMF_type,'VMF1'),strfind(VMF_type,'VMF1')])
           
                end  %//if all([isempty(VMF_type),isempty(grid_res)])
      
             else %IF VMF grid FILES NOT FOUND
            
                  %SAVE GRID FILE VERSION & TYPE WITH EMPTY MATRIX([])
                  setappdata(0,'VMFgrid_type_ppp',[])
                  setappdata(0,'VMFgrid_res_ppp',[])
        
                  beep %Give a beep sound
                  errmsg3{1}=sprintf('No VMF grid file(s) found in goGPS directory.\n');
                  errmsg3{2}=sprintf('System will resort to GPT2w / GPT3 grid file(s) instead.\n');
                  errmsg3{3}=sprintf('Mean while, VMF grid(s) can as well be downloaded from : http://vmf.geo.tuwien.ac.at/trop_products/GRID/\n');
                  warndlg(errmsg3,'VMF grid file Error','modal')
        
                  %READ FROM GPT3 GRID FILE(using GPT3 1° x 1° model)
                  grid_MF_ppp=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE 
                  gridRES_MF_ppp = 1;%GRID RESOLUTION
        
                  MFh_ppp = 'VMF3(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
                  MFw_ppp = 'VMF3(1° x 1°)';%WET MAPPING FUNCTION MODEL 
                  GPTmodel_ppp = 'GPT3'; %GPT3 MODEL 
                   
                  if isempty(grid_MF_ppp) %if GPT3 model grid is unavailable
            
                     %READ FROM GPT2w GRID FILE(using GPT2w 1° x 1° model)
                     grid_MF_ppp = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
         
                     MFh_ppp = 'VMF1(1° x 1°)'; %HYDROSTATIC MAPPING FUNCTION MODEL 
                     MFw_ppp = 'VMF1(1° x 1°)'; %WET MAPPING FUNCTION MODEL 
                     GPTmodel_ppp = 'GPT2w'; %GPT2w MODEL 
                     
                  end                
   
             end  %//if VMF_grid_found == 1
   
            %SAVED GPT2w/GPT3 GRID VALUES & RESOLUTION
            %HYDROSTATIC
            setappdata(0,'grid_MFh_ppp',grid_MF_ppp)
            setappdata(0,'gridRES_MFh_ppp',gridRES_MF_ppp)
      
            %WET
            setappdata(0,'grid_MFw_ppp',grid_MF_ppp)
            setappdata(0,'gridRES_MFw_ppp',gridRES_MF_ppp)
     
            %GPT MODEL
            setappdata(0,'GPTmodel_ppp',GPTmodel_ppp)
            
          elseif POPval == 2
                 MFh_ppp = 'GMF';
                 MFw_ppp = 'GMF';
       
          elseif POPval == 1
                 MFh_ppp = 'NMF';
                 MFw_ppp = 'NMF';   
    
          end %//if POPval == 3  
     
          %SAVE MAPPING FUNCTION MODELS
          setappdata(0,'MFh_model_ppp',MFh_ppp)
          setappdata(0,'MFw_model_ppp',MFw_ppp)
              
      end %//if all([POPval == 1,any([strcmpi(getappdata(0,'MFh_model_ppp'),'NMF'),strcmpi(getappdata(0,'MFw_model_ppp'),'NMF')])])     
   
   else %if getappdata(0,'PPPmf_POPval') ~= POPval  
    
     if POPval == 3
        
        %GET SAVED goGPS GUI FIGURE HANDLES
        handles = getappdata(0,'handle');  
         
        %CHECK IF VMF GRID FILEs ARE AVAILABLE IN goGPS
        %Call for the 'SearchVMFgrids.m' fxn
        [VMF_grid_found,VMFgrids] = SearchVMFgrids(handles); 
                        
        %*********SAVE STRUCT GRID FILES
        setappdata(0,'VMFgrids_ppp',VMFgrids)
         
        %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
        setappdata(0,'VMF_grid_found_ppp',VMF_grid_found)
   
   
        if VMF_grid_found == 1 %IF VMF grid FILES FOUND,LET USER CHOOSE GRID TYPE & RESOLUTION
                               %I.E.['VMF1(2°x2.5°)','VMF3(1°x1°)','VMF3(5°x5°)']
                               
           [VMF_type,grid_res] = VMFgrid_type(); 
      
           %SAVE GRID FILE VERSION & TYPE
           setappdata(0,'VMFgrid_type_ppp',VMF_type)
           setappdata(0,'VMFgrid_res_ppp',grid_res)
      
           grid_MF_ppp = [];%ASSIGN EMPTY([]) MATRIX TO GPT series grid values
           gridRES_MF_ppp = [];
           GPTmodel_ppp = []; %GPT MODEL (EMPTY([])
            
           %*********CHECK IF VMF_type,grid_res ARE BOTH EMPTY
           %NOTE:
           %     VMF_type,grid_res ARE ASSIGNED EMPTY WHEN USER CANCEL SELECTION
           %     OR CLOSE DIALOGUE BOX FOR THE SELECTION OF VMF GRID TYPE &
           %     RESOLUTION.IF IT HAPPENS SO, THE SELECTED TROPO DELAY MODEL,
           %     'VMF gridded ZTD' WOULD BE REPLACED WITH SAASTEMOIN MODEL
           %     OR THE POPUP MENU WILL BE SET TO SAASTAMOINEN AS DEFAULT
      
           if all([isempty(VMF_type),isempty(grid_res)])%IF USER CANCEL SELECTION/CLOSE DIALOGUE BOX
          
              beep %Give a beep sound
              warnmsg1{1}=sprintf('VMF grid type & resolution has not been selected .\n');
              warnmsg1{2}=sprintf('Try Again and select VMF grid type & resolution.\n');        
              warndlg(warnmsg1,'VMF grid file Error','modal')
         
              %READ FROM GPT3 GRID FILE(using GPT3 1° x 1° model)
              grid_MF_ppp=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDVALUES  
              gridRES_MF_ppp = 1;%GRID RESOLUTION
         
              MFh_ppp = 'VMF3(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
              MFw_ppp = 'VMF3(1° x 1°)';%WET MAPPING FUNCTION MODEL 
              GPTmodel_ppp = 'GPT3'; %GPT3 MODEL 
               
              if isempty(grid_MF_ppp) %if GPT3 model grid is unavailable
            
                 %READ FROM GPT2w GRID FILE(using GPT2w 1° x 1° model)
                 grid_MF_ppp = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
                 gridRES_MF_ppp = 1;%GRID RESOLUTION
         
                 MFh_ppp = 'VMF1(1° x 1°)'; %HYDROSTATIC MAPPING FUNCTION MODEL 
                 MFw_ppp = 'VMF1(1° x 1°)'; %WET MAPPING FUNCTION MODEL 
                 GPTmodel_ppp = 'GPT2w'; %GPT2w MODEL 
                  
              end         
            
           else %IF USER DOES NOT CANCEL SELECTION/CLOSE DIALOGUE BOX
                if any([strcmpi(VMF_type,'VMF1'),strfind(VMF_type,'VMF1')])
               
                   MFh_ppp = 'VMF1(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
                   MFw_ppp = 'VMF1(1° x 1°)';%WET MAPPING FUNCTION MODEL
              
                elseif any([strcmpi(VMF_type,'VMF3'),strfind(VMF_type,'VMF3')])
                  
                       if grid_res == 1
                      
                          MFh_ppp = 'VMF3(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
                          MFw_ppp = 'VMF3(1° x 1°)';%WET MAPPING FUNCTION MODEL
                     
                       elseif grid_res == 5
                      
                              MFh_ppp = 'VMF3(5° x 5°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
                              MFw_ppp = 'VMF3(5° x 5°)';%WET MAPPING FUNCTION MODEL
                       
                       end  
                      
                end %//if any([strcmpi(VMF_type,'VMF1'),strfind(VMF_type,'VMF1')])
           
           end  %//if all([isempty(VMF_type),isempty(grid_res)])
      
        else  %IF VMF grid FILES NOT FOUND
            
             %SAVE GRID FILE VERSION & TYPE WITH EMPTY MATRIX([])
             setappdata(0,'VMFgrid_type_ppp',[])
             setappdata(0,'VMFgrid_res_ppp',[])
        
             beep %Give a beep sound
             errmsg3{1}=sprintf('No VMF grid file(s) found in goGPS directory.\n');
             errmsg3{2}=sprintf('System will resort to GPT2w / GPT3 grid file(s) instead.\n');
             errmsg3{3}=sprintf('Mean while, VMF grid(s) can as well be downloaded from : http://vmf.geo.tuwien.ac.at/trop_products/GRID/\n');
             warndlg(errmsg3,'VMF grid file Error','modal')
        
            %READ FROM GPT3 GRID FILE(using GPT3 1° x 1° model)
            grid_MF_ppp=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE 
            gridRES_MF_ppp = 1;%GRID RESOLUTION
        
            MFh_ppp = 'VMF3(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
            MFw_ppp = 'VMF3(1° x 1°)';%WET MAPPING FUNCTION MODEL 
            GPTmodel_ppp = 'GPT3'; %GPT3 MODEL
            
            if isempty(grid_MF_ppp) %if GPT3 model grid is unavailable
            
               %READ FROM GPT2w GRID FILE(using GPT2w 1° x 1° model)
               grid_MF_ppp = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
         
               MFh_ppp = 'VMF1(1° x 1°)'; %HYDROSTATIC MAPPING FUNCTION MODEL 
               MFw_ppp = 'VMF1(1° x 1°)'; %WET MAPPING FUNCTION MODEL 
               GPTmodel_ppp = 'GPT2w'; %GPT2w MODEL
               
            end               
   
        end  %//if VMF_grid_found == 1
   
        %SAVED GPT2w/GPT3 GRID VALUES & RESOLUTION
        %HYDROSTATIC
        setappdata(0,'grid_MFh_ppp',grid_MF_ppp)
        setappdata(0,'gridRES_MFh_ppp',gridRES_MF_ppp)
      
        %WET
        setappdata(0,'grid_MFw_ppp',grid_MF_ppp)
        setappdata(0,'gridRES_MFw_ppp',gridRES_MF_ppp)
     
        %GPT MODEL
        setappdata(0,'GPTmodel_ppp',GPTmodel_ppp)
        
     elseif POPval == 2
            MFh_ppp = 'GMF';
            MFw_ppp = 'GMF';
       
     elseif POPval == 1
            MFh_ppp = 'NMF';
            MFw_ppp = 'NMF';   
    
     end  
     
     %SAVE MAPPING FUNCTION MODELS
     setappdata(0,'MFh_model_ppp',MFh_ppp)
     setappdata(0,'MFw_model_ppp',MFw_ppp)
     
     
   end %//if getappdata(0,'PPPmf_POPval') == POPval
       
else %IF MAPPING FUNCTION MODELS WERE NOT SELECTED 
    
     if POPval == 3
        
        %GET SAVED goGPS GUI FIGURE HANDLES
        handles = getappdata(0,'handle');  
         
        %CHECK IF VMF GRID FILEs ARE AVAILABLE IN goGPS
        %Call for the 'SearchVMFgrids.m' fxn
        [VMF_grid_found,VMFgrids] = SearchVMFgrids(handles); 
                        
        %*********SAVE STRUCT GRID FILES
        setappdata(0,'VMFgrids_ppp',VMFgrids)
         
        %SAVE VMF grid FILE FOUND/NOT FOUND INDICATOR
        setappdata(0,'VMF_grid_found_ppp',VMF_grid_found)
   
   
        if VMF_grid_found == 1 %IF VMF grid FILES FOUND,LET USER CHOOSE GRID TYPE & RESOLUTION
                               %I.E.['VMF1(2°x2.5°)','VMF3(1°x1°)','VMF3(5°x5°)']
                               
           [VMF_type,grid_res] = VMFgrid_type(); 
      
           %SAVE GRID FILE VERSION & TYPE
           setappdata(0,'VMFgrid_type_ppp',VMF_type)
           setappdata(0,'VMFgrid_res_ppp',grid_res)
      
           grid_MF_ppp = [];%ASSIGN EMPTY([]) MATRIX TO GPT series grid values
           gridRES_MF_ppp = [];
           GPTmodel_ppp = []; %GPT3 MODEL (EMPTY([])
            
           %*********CHECK IF VMF_type,grid_res ARE BOTH EMPTY
           %NOTE:
           %     VMF_type,grid_res ARE ASSIGNED EMPTY WHEN USER CANCEL SELECTION
           %     OR CLOSE DIALOGUE BOX FOR THE SELECTION OF VMF GRID TYPE &
           %     RESOLUTION.IF IT HAPPENS SO, THE SELECTED TROPO DELAY MODEL,
           %     'VMF gridded ZTD' WOULD BE REPLACED WITH SAASTEMOIN MODEL
           %     OR THE POPUP MENU WILL BE SET TO SAASTAMOINEN AS DEFAULT
      
           if all([isempty(VMF_type),isempty(grid_res)])%IF USER CANCEL SELECTION/CLOSE DIALOGUE BOX
          
              beep %Give a beep sound
              warnmsg1{1}=sprintf('VMF grid type & resolution has not been selected .\n');
              warnmsg1{2}=sprintf('Try Again and select VMF grid type & resolution.\n');        
              warndlg(warnmsg1,'VMF grid file Error','modal')
         
              %READ FROM GPT3 GRID FILE(using GPT3 1° x 1° model)
              grid_MF_ppp=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDVALUES  
              gridRES_MF_ppp = 1;%GRID RESOLUTION
         
              MFh_ppp = 'VMF3(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
              MFw_ppp = 'VMF3(1° x 1°)';%WET MAPPING FUNCTION MODEL 
              GPTmodel_ppp = 'GPT3'; %GPT3 MODEL 
              
              if isempty(grid_MF_ppp) %if GPT3 model grid is unavailable
            
                 %READ FROM GPT2w GRID FILE(using GPT2w 1° x 1° model)
                 grid_MF_ppp = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
                 gridRES_MF_ppp = 1;%GRID RESOLUTION
         
                 MFh_ppp = 'VMF1(1° x 1°)'; %HYDROSTATIC MAPPING FUNCTION MODEL 
                 MFw_ppp = 'VMF1(1° x 1°)'; %WET MAPPING FUNCTION MODEL 
                 GPTmodel_ppp = 'GPT2w'; %GPT2w MODEL 
                  
              end         
            
           else %IF USER DOES NOT CANCEL SELECTION/CLOSE DIALOGUE BOX
                if any([strcmpi(VMF_type,'VMF1'),strfind(VMF_type,'VMF1')])
               
                   MFh_ppp = 'VMF1(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
                   MFw_ppp = 'VMF1(1° x 1°)';%WET MAPPING FUNCTION MODEL
              
                elseif any([strcmpi(VMF_type,'VMF3'),strfind(VMF_type,'VMF3')])
                  
                       if grid_res == 1
                      
                          MFh_ppp = 'VMF3(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
                          MFw_ppp = 'VMF3(1° x 1°)';%WET MAPPING FUNCTION MODEL
                     
                       elseif grid_res == 5
                      
                              MFh_ppp = 'VMF3(5° x 5°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
                              MFw_ppp = 'VMF3(5° x 5°)';%WET MAPPING FUNCTION MODEL
                       
                       end  
                      
                end %//if any([strcmpi(VMF_type,'VMF1'),strfind(VMF_type,'VMF1')])
           
           end  %//if all([isempty(VMF_type),isempty(grid_res)])
      
        else  %IF VMF grid FILES NOT FOUND
            
             %SAVE GRID FILE VERSION & TYPE WITH EMPTY MATRIX([])
             setappdata(0,'VMFgrid_type_ppp',[])
             setappdata(0,'VMFgrid_res_ppp',[])
        
             beep %Give a beep sound
             errmsg3{1}=sprintf('No VMF grid file(s) found in goGPS directory.\n');
             errmsg3{2}=sprintf('System will resort to GPT2w / GPT3 grid file(s) instead.\n');
             errmsg3{3}=sprintf('Mean while, VMF grid(s) can as well be downloaded from : http://vmf.geo.tuwien.ac.at/trop_products/GRID/\n');
             warndlg(errmsg3,'VMF grid file Error','modal')
        
            %READ FROM GPT3 GRID FILE(using GPT3 1° x 1° model)
            grid_MF_ppp=readGPTgrid('gpt3_1.mat','GPT3',1);%GRIDFILE 
            gridRES_MF_ppp = 1;%GRID RESOLUTION
        
            MFh_ppp = 'VMF3(1° x 1°)';%HYDROSTATIC MAPPING FUNCTION MODEL 
            MFw_ppp = 'VMF3(1° x 1°)';%WET MAPPING FUNCTION MODEL 
            GPTmodel_ppp = 'GPT3'; %GPT3 MODEL 
            
            if isempty(grid_MF_ppp) %if GPT3 model grid is unavailable
            
               %READ FROM GPT2w GRID FILE(using GPT2w 1° x 1° model)
               grid_MF_ppp = readGPTgrid('gpt2_1w.mat','GPT2w',1);%GRIDFILE
         
               MFh_ppp = 'VMF1(1° x 1°)'; %HYDROSTATIC MAPPING FUNCTION MODEL 
               MFw_ppp = 'VMF1(1° x 1°)'; %WET MAPPING FUNCTION MODEL 
               GPTmodel_ppp = 'GPT2w'; %GPT2w MODEL 
               
            end               
   
        end  %//if VMF_grid_found == 1
   
        %SAVED GPT2w/GPT3 GRID VALUES & RESOLUTION
        %HYDROSTATIC
        setappdata(0,'grid_MFh_ppp',grid_MF_ppp)
        setappdata(0,'gridRES_MFh_ppp',gridRES_MF_ppp)
      
        %WET
        setappdata(0,'grid_MFw_ppp',grid_MF_ppp)
        setappdata(0,'gridRES_MFw_ppp',gridRES_MF_ppp)
        
        %GPT MODEL
        setappdata(0,'GPTmodel_ppp',GPTmodel_ppp)
        
     elseif POPval == 2
            MFh_ppp = 'GMF';
            MFw_ppp = 'GMF';
       
     elseif POPval == 1
            MFh_ppp = 'NMF';
            MFw_ppp = 'NMF';   
    
     end  
     
     %SAVE MAPPING FUNCTION MODELS
     setappdata(0,'MFh_model_ppp',MFh_ppp)
     setappdata(0,'MFw_model_ppp',MFw_ppp)
     
     
end %//if all([~isempty(getappdata(0,'MFh_model_ppp')),~isempty(getappdata(0,'MFw_model_ppp'))]) 
 
%SAVE MF INDEX FROM POPUP MENU
setappdata(0,'PPPmf_POPval',POPval)

close(gcbf)%close figure

%=========================================END OF MF_ppp.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 



%              -----------------------------------------------------------
%(10)***********SUB-ROUTINE TO SAVE THE VARIOUS ATMOSPHERIC DELAY  MODELLING
%              GUI COMPONENTS 
%              ------------------------------------------------------------
function saveGUIState(handles)
%**************************************************************************
%DESCRIPTION:                                                              * 
%           "saveGUIState" is a Sub-routine that saves the  various        *
%            Atmospheric Delay Modelling GUI  components in goGPS.         *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================
% WRITTEN BY: 
%            OSAH SAMUEL, MSC GEOMATIC ENGINEERING (PhD STUDENT)           +
%            Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com                +
%            Phone:+233(0)246137410 / +233(0)509438484                     +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%**************************************************************************+
%********GUI NAME
%******GET UICONTROLS
%--------------------------------------------------------------------------
%*****GET MASTER STATION INFO
state.master_pos_ena=get(handles.master_pos,'enable');%GET STATE OF MASTER STATION POSITION INPUT(ENABLE OR NOT ENABLE) 
state.master_pos_val=get(handles.master_pos,'value');%GET MASTER STATION POSITION INPUT(VALUE:1 OR 0) 

%GET POSITION TYPE/FORMAT(I.E.ECEF(XYZ) OR GEOGRAPHIC(LAT,LON,h)
state.crs_ena=get(handles.crs,'enable');%MASTER POSITION FORMAT/TYPE POPUP MENU CONTENTS
state.crs_str=get(handles.crs,'string');%MASTER POSITION FORMAT/TYPE POPUP MENU CONTENTS
state.crs_val=get(handles.crs,'value');%POSITION(VALUE) OF SELECTED MASTER POSITION FORMAT/TYPE

%GET STATE OF POSITION INPUT STATIC & EDIT TEXT
%A.ECEF(XYZ) INPUT
%1.STATIC TEXT
state.text_master_X_ena=get(handles.text_master_X,'enable');%GET TEXT STATE(ENABLE OR NOT ENABLE)
state.text_master_Y_ena=get(handles.text_master_Y,'enable');%GET TEXT STATE(ENABLE OR NOT ENABLE)
state.text_master_Z_ena=get(handles.text_master_Z,'enable');%GET TEXT STATE(ENABLE OR NOT ENABLE)

%UNIT
state.text_master_X_unit_ena=get(handles.text_master_X_unit,'enable');%GET TEXT UNIT STATE(ENABLE OR NOT ENABLE)
state.text_text_master_Y_unit_ena=get(handles.text_master_Y_unit,'enable');%GET TEXT UNIT STATE(ENABLE OR NOT ENABLE)
state.text_master_Z_unit_ena=get(handles.text_master_Z_unit,'enable');%GET TEXT UNIT STATE(ENABLE OR NOT ENABLE)

%1.EDIT TEXT
state.master_X_ena=get(handles.master_X,'enable');%GET EDIT TEXT STATE(ENABLE OR NOT ENABLE)
state.master_Y_ena=get(handles.master_Y,'enable');%GET EDIT TEXT STATE(ENABLE OR NOT ENABLE)
state.master_Z_ena=get(handles.master_Z,'enable');%GET EDIT TEXT STATE(ENABLE OR NOT ENABLE)
state.master_X_str=get(handles.master_X,'string');%GET EDIT TEXT STATE VALUE (I.E. STRING)
state.master_Y_str=get(handles.master_Y,'string');%GET EDIT TEXT STATE VALUE (I.E. STRING)
state.master_Z_str=get(handles.master_Z,'string');%GET EDIT TEXT STATE VALUE (I.E. STRING)


%B.GEODETIC(LAT,LON h) INPUT
%1.STATIC TEXT
state.text_master_lat_ena=get(handles.text_master_lat,'enable');%GET TEXT STATE(ENABLE OR NOT ENABLE)
state.text_master_lon_ena=get(handles.text_master_lon,'enable');%GET TEXT STATE(ENABLE OR NOT ENABLE)
state.text_master_h_ena=get(handles.text_master_h,'enable');%GET TEXT STATE(ENABLE OR NOT ENABLE)

%UNIT
state.text_master_lat_unit_ena=get(handles.text_master_lat_unit,'enable');%GET TEXT UNIT STATE(ENABLE OR NOT ENABLE)
state.text_text_master_lon_unit_ena=get(handles.text_master_lon_unit,'enable');%GET TEXT UNIT STATE(ENABLE OR NOT ENABLE)
state.text_master_h_unit_ena=get(handles.text_master_h_unit,'enable');%GET TEXT UNIT STATE(ENABLE OR NOT ENABLE)

%1.EDIT TEXT
state.master_lat_ena=get(handles.master_lat,'enable');%GET EDIT TEXT STATE(ENABLE OR NOT ENABLE)
state.master_lon_ena=get(handles.master_lon,'enable');%GET EDIT TEXT STATE(ENABLE OR NOT ENABLE)
state.master_h_ena=get(handles.master_h,'enable');%GET EDIT TEXT STATE(ENABLE OR NOT ENABLE)
state.master_lat_str=get(handles.master_lat,'string');%GET EDIT TEXT STATE VALUE (I.E. STRING)
state.master_lon_str=get(handles.master_lon,'string');%GET EDIT TEXT STATE VALUE (I.E. STRING)
state.master_h_str=get(handles.master_h,'string');%GET EDIT TEXT STATE VALUE (I.E. STRING)


%%1.IONOSPHERIC DELAY MODELLING UICONTROLS
state.pop_lIono_correction_str=get(handles.pop_lIono_correction,'string');%IONOSPHERIC POPUP MENU CONTENTS
state.pop_lIono_correction_val=get(handles.pop_lIono_correction,'value');%POSITION(VALUE) OF SELECTED IONOSPHERIC MODEL
%--------------------------------------------------------------------------
%2.TROPOSPHERIC DELAY MODELLING UICONTROLS
%I.TROPOSPHERIC RADIO BUTTON       
state.rb_Tropo_combine=get(handles.rb_Tropo_combine,'value');%COMBINE(DRY+WET)
state.rb_Tropo_separate=get(handles.rb_Tropo_separate,'value');%SEPARATE(DRY OR WET)

%II.POPUP MENUs
%A.COMBINE(DRY+WET) MODELS
state.pop_lTropo_combine_str=get(handles.pop_lTropo_combine,'string');% POPUP MENU CONTENTS
state.pop_lTropo_combine_val=get(handles.pop_lTropo_combine,'value');%POSITION(VALUE) OF SELECTED COMBINE TROPO MODEL
state.pop_lTropo_combine_enable=get(handles.pop_lTropo_combine,'Enable');%POPUP MENUS STATE(ENABLE OR NOT ENABLE)

%TEXT
state.text_Tropo_combine=get(handles.text_Tropo_combine,'Enable');% TEXT STATE(ENABLE OR NOT ENABLE) FOR  COMBINE(DRY+WET) MODELS 

%B.SEPARATE(DRY OR WET) MODELS
%I.HYDROSTATIC MODELs
state.pop_ITropo_hydrostatic_str=get(handles.pop_ITropo_hydrostatic,'string');% POPUP MENU CONTENTS
state.pop_ITropo_hydrostatic_val=get(handles.pop_ITropo_hydrostatic,'value');%POSITION(VALUE) OF SELECTED HYDROSTATIC MODEL
state.pop_ITropo_hydrostatic_enable=get(handles.pop_ITropo_hydrostatic,'Enable');%POPUP MENUS STATE(ENABLE OR NOT ENABLE)

%TEXT
state.text_Tropo_hydrostatic=get(handles.text_Tropo_hydrostatic,'Enable');% TEXT STATE(ENABLE OR NOT ENABLE) FOR  HYDROSTATIC MODELS 

%I.WET MODELs
state.pop_ITropo_wet_str=get(handles.pop_ITropo_wet,'string');% POPUP MENU CONTENTS
state.pop_ITropo_wet_val=get(handles.pop_ITropo_wet,'value');%POSITION(VALUE) OF SELECTED WET MODEL
state.pop_ITropo_wet_enable=get(handles.pop_ITropo_wet,'Enable');%POPUP MENUS STATE(ENABLE OR NOT ENABLE)

%TEXT
state.text_Tropo_wet=get(handles.text_Tropo_wet,'Enable');% TEXT STATE(ENABLE OR NOT ENABLE) FOR  WET MODELS 
%--------------------------------------------------------------------------

%3.********MAPPING FUNCTIONS(MF)
%I.RADIO BUTTONS(rb)
state.rb_model_MF=get(handles.rb_model_MF,'value');%Model Mapping Function rb
state.rb_different_MF=get(handles.rb_different_MF,'value');% other Mapping Models rb
state.rb_different_MF_enable=get(handles.rb_different_MF,'Enable');% RADIO BUTTON STATE(ENABLE OR NOT ENABLE)  
state.rb_model_MF_enable=get(handles.rb_model_MF,'Enable');% RADIO BUTTON STATE(ENABLE OR NOT ENABLE)  

%TEXT
state.text_mapping_function=get(handles.text_mapping_function,'Enable');% TEXT STATE(ENABLE OR NOT ENABLE) FOR  MF RB

%II. POPUP MENUS
%I.HYDROSTATIC MF MODELs
state.pop_IHydrostic_MF_str=get(handles.pop_IHydrostic_MF,'string');% POPUP MENU CONTENTS
state.pop_IHydrostic_MF_val=get(handles.pop_IHydrostic_MF,'value');%POSITION(VALUE) OF SELECTED HYDROSTATIC MF MODEL
state.pop_IHydrostic_MF_enable=get(handles.pop_IHydrostic_MF,'Enable');%POPUP MENUS STATE(ENABLE OR NOT ENABLE)

%TEXT
state.text__Hydrostatic_MF=get(handles.text__Hydrostatic_MF,'Enable');% TEXT STATE(ENABLE OR NOT ENABLE) FOR  HYDROSTATIC MODELS 

%I.WET MF MODELs
state.pop_IWet_MF_str=get(handles.pop_IWet_MF,'string');% POPUP MENU CONTENTS
state.pop_IWet_MF_val=get(handles.pop_IWet_MF,'value');%POSITION(VALUE) OF SELECTED WET MF MODEL
state.pop_IWet_MF_enable=get(handles.pop_IWet_MF,'Enable');%POPUP MENUS STATE(ENABLE OR NOT ENABLE)

%TEXT
state.text__Wet_MF=get(handles.text__Wet_MF,'Enable');% TEXT STATE(ENABLE OR NOT ENABLE) FOR  WET MODELS 
%--------------------------------------------------------------------------
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%3.SOURCE OF METEOROLOGICAL PARAMETERS
%I. POPUP MENUS
state.pop_source_metpara_str=get(handles.pop_source_metpara,'string');% POPUP MENU CONTENTS
state.pop_source_metpara_val=get(handles.pop_source_metpara,'value');%POPUP MENU CONTENTS [SOURCE OF MET PARA]
state.pop_source_metpara_vis=get(handles.pop_source_metpara,'visible');%POPUP MENU VISIBILITY(On/Off)
state.pop_source_metpara_ena=get(handles.pop_source_metpara,'enable');%POPUP MENU ENABILITY(On/Off)

state.pop_met_manual_str=get(handles.pop_met_manual,'string');%POSITION(VALUE) OF SELECTED SOURCE OF USER DEFINED MET PARA
state.pop_met_manual_val=get(handles.pop_met_manual,'value');%POSITION(VALUE) OF SELECTED SOURCE OF USER DEFINED MET PARA
state.pop_met_manual_vis=get(handles.pop_met_manual,'visible');%SOURCE OF USER DEFINED MET PARA POPUP MENU visibility(on or off)
state.pop_met_manual_ena=get(handles.pop_met_manual,'enable');%SOURCE OF USER DEFINED MET PARA POPUP MENU enability(on or off)

%II.GRID RESOLUTION RADIO BUTTONS
state.rb_grid_resolution_1_val=get(handles.rb_grid_resolution_1,'value');%1° grid resolution (VALUE) OF SELECTION (0 OR 1)
state.rb_grid_resolution_5_val=get(handles.rb_grid_resolution_5,'value');% 5° grid resolution (VALUE) OF SELECTION(0 OR 1)
state.rb_grid_resolution_1_vis=get(handles.rb_grid_resolution_1,'visible');%1° grid resolution visibility(on or off)
state.rb_grid_resolution_5_vis=get(handles.rb_grid_resolution_5,'visible');% 5° grid resolution visibility(on or off)
state.rb_grid_resolution_1_fcolor=get(handles.rb_grid_resolution_1,'ForegroundColor');% 1° grid resolution ForegroundColor
state.rb_grid_resolution_5_fcolor=get(handles.rb_grid_resolution_5,'ForegroundColor');% 5° grid resolution ForegroundColor

%TEXT
state.text_grid_resolution=get(handles.text_grid_resolution,'Visible');% TEXT STATE(VISIBLE OR NOT VISIBLE)
state.text_metmanual_source_vis=get(handles.text_metmanual_source,'Visible');% TEXT STATE(VISIBLE OR NOT VISIBLE)
state.text_metmanual_source_ena=get(handles.text_metmanual_source,'enable');%TEXT STATE [enability(on or off)]


%4. EXTRACT VMF GRIDDED ZENITH DELAYS CHECKBOX
state.cb_extract_VMF_ZTDs_val=get(handles.cb_extract_VMF_ZTDs,'Value');%CHECKBOX (VALUE) OF SELECTION(0 OR 1)

%5.RETRIEVAL OF PRECIPITABLE WATER VAPOUR(PWV)
%A.retrievePWV CHECKBOX
state.cb_retrievePWV_val=get(handles.cb_retrievePWV,'Value');%CHECKBOX (VALUE) OF SELECTION(0 OR 1)

%B.WEIGHTED MEAN TEMPERATURE(Tm) POPUP MENU
state.pop_Tm_str=get(handles.pop_Tm,'string');% POPUP MENU CONTENTS
state.pop_Tm_val=get(handles.pop_Tm,'value');%POSITION(VALUE) OF SELECTED Tm MODEL
state.pop_Tm_vis=get(handles.pop_Tm,'visible');%Get POPUP MENU visibility(on or off)
state.pop_Tm_ena=get(handles.pop_Tm,'enable');%Get POPUP MENU enability(on or off)

%C.GRID RESOLUTION RADIO BUTTONS
state.rb_grid_resolution_1_pwv_val=get(handles.rb_grid_resolution_1_pwv,'value');%1° grid resolution (VALUE) OF SELECTION (0 OR 1)
state.rb_grid_resolution_5_pwv_val=get(handles.rb_grid_resolution_5_pwv,'value');% 5° grid resolution (VALUE) OF SELECTION(0 OR 1)
state.rb_grid_resolution_1_pwv_vis=get(handles.rb_grid_resolution_1_pwv,'visible');%1° grid resolution visibility(on or off)
state.rb_grid_resolution_5_pwv_vis=get(handles.rb_grid_resolution_5_pwv,'visible');% 5° grid resolution visibility(on or off)
state.rb_grid_resolution_1_pwv_fcolor=get(handles.rb_grid_resolution_1_pwv,'ForegroundColor');% 1° grid resolution ForegroundColor
state.rb_grid_resolution_5_pwv_fcolor=get(handles.rb_grid_resolution_5_pwv,'ForegroundColor');% 5° grid resolution ForegroundColor

%D.TEXT
state.text_Tm_vis=get(handles.text_Tm,'visible');% TEXT visibility(VISIBLE OR NOT VISIBLE)
state.text_Tm_ena=get(handles.text_Tm,'enable');% TEXT enability(on or off)
state.text_grid_resolution_pwv=get(handles.text_grid_resolution_pwv,'Visible');% TEXT STATE(VISIBLE OR NOT VISIBLE)
%--------------------------------------------------------------------------
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

save('atmosphere.mat','state');
%=========================================END OF saveGUIState.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 


%              -----------------------------------------------------------
%(11)***********SUB-ROUTINE TO RETRIEVE THE SAVED VARIOUS ATMOSPHERIC DELAY  
%               MODELLING GUI COMPONENTS 
%              ------------------------------------------------------------
function loadGUIState(handles)
%**************************************************************************
%DESCRIPTION:                                                              * 
%           "loadGUIState" is a Sub-routine that Loads the  various        *
%            Atmospheric Delay Modelling GUI  components in goGPS Set them *
%            as previously selected by user.
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%==========================================================================
% WRITTEN BY: 
%            OSAH SAMUEL, MSC GEOMATIC ENGINEERING (PhD STUDENT)           +
%            Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com                +
%            Phone:+233(0)246137410 / +233(0)509438484                     +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%**************************************************************************+
%=========================================END OF loadGUIState.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
%Load save file-'atmosphere.mat')

if exist('atmosphere.mat','file')

  %LOAD fileName
  load('atmosphere.mat')
  
%********GUI NAME
%******SET UICONTROLS
%--------------------------------------------------------------------------

%*****GET MASTER STATION INFO
set(handles.master_pos,'enable',state.master_pos_ena);%STATE OF MASTER STATION POSITION INPUT(ENABLE OR NOT ENABLE) 
set(handles.master_pos,'value',state.master_pos_val);%VALUE OF  MASTER STATION POSITION INPUT(VALUE:1 OR 0) 

%GET POSITION TYPE/FORMAT(I.E.ECEF(XYZ) OR GEOGRAPHIC(LAT,LON,h)
set(handles.crs,'enable',state.crs_ena);%STATE OF MASTER POSITION FORMAT/TYPE POPUP MENU CONTENTS(ENABLE OR NOT ENABLE) 
set(handles.crs,'string',state.crs_str);%MASTER POSITION FORMAT/TYPE POPUP MENU CONTENTS
set(handles.crs,'value',state.crs_val);%POSITION(VALUE) OF SELECTED MASTER POSITION FORMAT/TYPE

%GET STATE OF POSITION INPUT STATIC & EDIT TEXT
%A.ECEF(XYZ) INPUT
%1.STATIC TEXT
set(handles.text_master_X,'enable',state.text_master_X_ena);%TEXT STATE(ENABLE OR NOT ENABLE)
set(handles.text_master_Y,'enable',state.text_master_Y_ena);%TEXT STATE(ENABLE OR NOT ENABLE)
set(handles.text_master_Z,'enable',state.text_master_Z_ena);%TEXT STATE(ENABLE OR NOT ENABLE)

%UNIT
set(handles.text_master_X_unit,'enable',state.text_master_X_unit_ena);%TEXT UNIT STATE(ENABLE OR NOT ENABLE)
set(handles.text_master_Y_unit,'enable',state.text_text_master_Y_unit_ena);%TEXT UNIT STATE(ENABLE OR NOT ENABLE)
set(handles.text_master_Z_unit,'enable',state.text_master_Z_unit_ena);%TEXT UNIT STATE(ENABLE OR NOT ENABLE)

%1.EDIT TEXT
set(handles.master_X,'enable',state.master_X_ena);%EDIT TEXT STATE(ENABLE OR NOT ENABLE)
set(handles.master_Y,'enable',state.master_Y_ena);%EDIT TEXT STATE(ENABLE OR NOT ENABLE)
set(handles.master_Z,'enable',state.master_Z_ena);%EDIT TEXT STATE(ENABLE OR NOT ENABLE)
set(handles.master_X,'string',state.master_X_str);%GEDIT TEXT STATE VALUE (I.E. STRING)
set(handles.master_Y,'string',state.master_Y_str);%EDIT TEXT STATE VALUE (I.E. STRING)
set(handles.master_Z,'string',state.master_Z_str);%EDIT TEXT STATE VALUE (I.E. STRING)


%B.GEODETIC(LAT,LON h) INPUT
%1.STATIC TEXT
set(handles.text_master_lat,'enable',state.text_master_lat_ena);%TEXT STATE(ENABLE OR NOT ENABLE)
set(handles.text_master_lon,'enable',state.text_master_lon_ena);%TEXT STATE(ENABLE OR NOT ENABLE)
set(handles.text_master_h,'enable',state.text_master_h_ena);%TEXT STATE(ENABLE OR NOT ENABLE)

%UNIT
set(handles.text_master_lat_unit,'enable',state.text_master_lat_unit_ena);%TEXT UNIT STATE(ENABLE OR NOT ENABLE)
set(handles.text_master_lon_unit,'enable',state.text_text_master_lon_unit_ena);%TEXT UNIT STATE(ENABLE OR NOT ENABLE)
set(handles.text_master_h_unit,'enable',state.text_master_h_unit_ena);%TEXT UNIT STATE(ENABLE OR NOT ENABLE)

%1.EDIT TEXT
set(handles.master_lat,'enable',state.master_lat_ena);%EDIT TEXT STATE(ENABLE OR NOT ENABLE)
set(handles.master_lon,'enable',state.master_lon_ena);%EDIT TEXT STATE(ENABLE OR NOT ENABLE)
set(handles.master_h,'enable',state.master_h_ena);%EDIT TEXT STATE(ENABLE OR NOT ENABLE)
set(handles.master_lat,'string',state.master_lat_str);%EDIT TEXT STATE VALUE (I.E. STRING)
set(handles.master_lon,'string',state.master_lon_str);%EDIT TEXT STATE VALUE (I.E. STRING)
set(handles.master_h,'string',state.master_h_str);%EDIT TEXT STATE VALUE (I.E. STRING)

%%1.IONOSPHERIC DELAY MODELLING UICONTROLS
set(handles.pop_lIono_correction,'string',state.pop_lIono_correction_str);%IONOSPHERIC POPUP MENU CONTENTS
% set(handles.pop_lIono_correction,'value',pop_lIono_correction_val);%POSITION(VALUE) OF SELECTED IONOSPHERIC MODEL
%--------------------------------------------------------------------------

%2.TROPOSPHERIC DELAY MODELLING UICONTROLS
%I.TROPOSPHERIC RADIO BUTTON       
set(handles.rb_Tropo_combine,'value',state.rb_Tropo_combine);%COMBINE(DRY+WET)
set(handles.rb_Tropo_separate,'value',state.rb_Tropo_separate);%SEPARATE(DRY OR WET)

%II.POPUP MENUs
%A.COMBINE(DRY+WET) MODELS
set(handles.pop_lTropo_combine,'string',state.pop_lTropo_combine_str);% POPUP MENU CONTENTS
set(handles.pop_lTropo_combine,'value',state.pop_lTropo_combine_val);%POSITION(VALUE) OF SELECTED COMBINE TROPO MODEL
set(handles.pop_lTropo_combine,'Enable',state.pop_lTropo_combine_enable);%POPUP MENUS STATE(ENABLE OR NOT ENABLE)

%TEXT
set(handles.text_Tropo_combine,'Enable',state.text_Tropo_combine);% TEXT STATE(ENABLE OR NOT ENABLE) FOR  COMBINE(DRY+WET) MODELS 

%B.SEPARATE(DRY OR WET) MODELS
%I.HYDROSTATIC MODELs
set(handles.pop_ITropo_hydrostatic,'string',state.pop_ITropo_hydrostatic_str);% POPUP MENU CONTENTS
set(handles.pop_ITropo_hydrostatic,'value',state.pop_ITropo_hydrostatic_val);%POSITION(VALUE) OF SELECTED HYDROSTATIC MODEL
set(handles.pop_ITropo_hydrostatic,'Enable',state.pop_ITropo_hydrostatic_enable);%POPUP MENUS STATE(ENABLE OR NOT ENABLE)

%TEXT
set(handles.text_Tropo_hydrostatic,'Enable',state.text_Tropo_hydrostatic);% TEXT STATE(ENABLE OR NOT ENABLE) FOR  HYDROSTATIC MODELS 

%I.WET MODELs
set(handles.pop_ITropo_wet,'string',state.pop_ITropo_wet_str);% POPUP MENU CONTENTS
set(handles.pop_ITropo_wet,'value',state.pop_ITropo_wet_val);%POSITION(VALUE) OF SELECTED WET MODEL
set(handles.pop_ITropo_wet,'Enable',state.pop_ITropo_wet_enable);%POPUP MENUS STATE(ENABLE OR NOT ENABLE)

%TEXT
set(handles.text_Tropo_wet,'Enable',state.text_Tropo_wet);% TEXT STATE(ENABLE OR NOT ENABLE) FOR  WET MODELS 
%--------------------------------------------------------------------------


%3.********MAPPING FUNCTIONS(MF)
%I.RADIO BUTTONS(rb)
set(handles.rb_model_MF,'value',state.rb_model_MF);%Model Mapping Function rb
set(handles.rb_different_MF,'value',state.rb_different_MF);% other Mapping Models rb
set(handles.rb_different_MF,'Enable',state.rb_different_MF_enable);% RADIO BUTTON STATE(ENABLE OR NOT ENABLE)  
set(handles.rb_model_MF,'Enable',state.rb_model_MF_enable);% RADIO BUTTON STATE(ENABLE OR NOT ENABLE)  

%TEXT
set(handles.text_mapping_function,'Enable',state.text_mapping_function);% TEXT STATE(ENABLE OR NOT ENABLE) FOR  MF RB

%II. POPUP MENUS
%I.HYDROSTATIC MF MODELs
set(handles.pop_IHydrostic_MF,'string',state.pop_IHydrostic_MF_str);% POPUP MENU CONTENTS
set(handles.pop_IHydrostic_MF,'value',state.pop_IHydrostic_MF_val);%POSITION(VALUE) OF SELECTED HYDROSTATIC MF MODEL
set(handles.pop_IHydrostic_MF,'Enable',state.pop_IHydrostic_MF_enable);%POPUP MENUS STATE(ENABLE OR NOT ENABLE)

%TEXT
set(handles.text__Hydrostatic_MF,'Enable',state.text__Hydrostatic_MF);% TEXT STATE(ENABLE OR NOT ENABLE) FOR  HYDROSTATIC MODELS 

%I.WET MF MODELs
set(handles.pop_IWet_MF,'string',state.pop_IWet_MF_str);% POPUP MENU CONTENTS
set(handles.pop_IWet_MF,'value',state.pop_IWet_MF_val);%POSITION(VALUE) OF SELECTED WET MF MODEL
set(handles.pop_IWet_MF,'Enable',state.pop_IWet_MF_enable);%POPUP MENUS STATE(ENABLE OR NOT ENABLE)

%TEXT
set(handles.text__Wet_MF,'Enable',state.text__Wet_MF);% TEXT STATE(ENABLE OR NOT ENABLE) FOR  WET MODELS 
%--------------------------------------------------------------------------
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%3.SOURCE OF METEOROLOGICAL PARAMETERS
%I. POPUP MENUS
set(handles.pop_source_metpara,'string',state.pop_source_metpara_str);% POPUP MENU CONTENTS
set(handles.pop_source_metpara,'value',state.pop_source_metpara_val);%POSITION(VALUE) OF SELECTED HYDROSTATIC MF MODEL
set(handles.pop_source_metpara,'visible',state.pop_source_metpara_vis);%POPUP MENU VISIBILITY(On/Off)
set(handles.pop_source_metpara,'enable',state.pop_source_metpara_ena);%POPUP MENU ENABILITY(On/Off)

%USER ENTERED/IMPORTED MET PARA SOURCE
set(handles.pop_met_manual,'string',state.pop_met_manual_str);%POPUP MENU CONTENTS[SOURCE OF USER DEFINED MET PARA]
set(handles.pop_met_manual,'value',state.pop_met_manual_val);%POSITION(VALUE) OF SELECTED SOURCE OF USER DEFINED MET PARA
set(handles.pop_met_manual,'visible',state.pop_met_manual_vis);%SOURCE OF USER DEFINED MET PARA POPUP MENU enability(on or off)
set(handles.pop_met_manual,'enable',state.pop_met_manual_ena)

%II.GRID RESOLUTION RADIO BUTTONS
set(handles.rb_grid_resolution_1,'value',state.rb_grid_resolution_1_val);%1° grid resolution (VALUE) OF SELECTION (0 OR 1)
set(handles.rb_grid_resolution_5,'value',state.rb_grid_resolution_5_val);% 5° grid resolution (VALUE) OF SELECTION(0 OR 1)
set(handles.rb_grid_resolution_1,'visible',state.rb_grid_resolution_1_vis);%1° grid resolution visibility(on or off)
set(handles.rb_grid_resolution_5,'visible',state.rb_grid_resolution_5_vis);% 5° grid resolution visibility(on or off)
set(handles.rb_grid_resolution_1,'ForegroundColor',state.rb_grid_resolution_1_fcolor);% 1° grid resolution ForegroundColor
set(handles.rb_grid_resolution_5,'ForegroundColor',state.rb_grid_resolution_5_fcolor);% 5° grid resolution ForegroundColor


%TEXT
set(handles.text_grid_resolution,'Visible',state.text_grid_resolution);% TEXT STATE(VISIBLE OR NOT VISIBLE)
set(handles.text_metmanual_source,'Visible',state.text_metmanual_source_vis);% TEXT STATE(VISIBLE OR NOT VISIBLE)
set(handles.text_metmanual_source,'Enable',state.text_metmanual_source_ena);% TEXT STATE(VISIBLE OR NOT VISIBLE)


%4. EXTRACT VMF GRIDDED ZENITH DELAYS CHECKBOX
set(handles.cb_extract_VMF_ZTDs,'Value',state.cb_extract_VMF_ZTDs_val);%CHECKBOX (VALUE) OF SELECTION(0 OR 1)

%5.RETRIEVAL OF PRECIPITABLE WATER VAPOUR(PWV)
%A.retrievePWV CHECKBOX
set(handles.cb_retrievePWV,'Value',state.cb_retrievePWV_val);%CHECKBOX (VALUE) OF SELECTION(0 OR 1)

%B.WEIGHTED MEAN TEMPERATURE(Tm) POPUP MENU
set(handles.pop_Tm,'string',state.pop_Tm_str);% POPUP MENU CONTENTS
set(handles.pop_Tm,'value',state.pop_Tm_val);%POSITION(VALUE) OF SELECTED Tm MODEL
set(handles.pop_Tm,'visible',state.pop_Tm_vis);%set visibility(on or off)
set(handles.pop_Tm,'enable',state.pop_Tm_ena);%set enability(on or off)

%C.GRID RESOLUTION RADIO BUTTONS
set(handles.rb_grid_resolution_1_pwv,'value',state.rb_grid_resolution_1_pwv_val);%1° grid resolution (VALUE) OF SELECTION (0 OR 1)
set(handles.rb_grid_resolution_5_pwv,'value',state.rb_grid_resolution_5_pwv_val);% 5° grid resolution (VALUE) OF SELECTION(0 OR 1)
set(handles.rb_grid_resolution_1_pwv,'visible',state.rb_grid_resolution_1_pwv_vis);%1° grid resolution visibility(on or off)
set(handles.rb_grid_resolution_5_pwv,'visible',state.rb_grid_resolution_5_pwv_vis);% 5° grid resolution visibility(on or off)
set(handles.rb_grid_resolution_1_pwv,'ForegroundColor',state.rb_grid_resolution_1_pwv_fcolor);% 1° grid resolution ForegroundColor
set(handles.rb_grid_resolution_5_pwv,'ForegroundColor',state.rb_grid_resolution_5_pwv_fcolor);% 5° grid resolution ForegroundColor

%D.TEXT
set(handles.text_Tm,'Visible',state.text_Tm_vis);% TEXT STATE(VISIBLE OR NOT VISIBLE)
set(handles.text_Tm,'enable',state.text_Tm_ena);% set enability(on or off)
set(handles.text_grid_resolution_pwv,'Visible',state.text_grid_resolution_pwv);% TEXT STATE(VISIBLE OR NOT VISIBLE)

%--------------------------------------------------------------------------
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

else
    
    return

end