%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

% clear all variables
% NOTE: using only 'clearvars' does not clear global variables, while using
% 'clear all' removes breakpoints
clearvars -global -except ini_settings_file use_gui; %  exceptions for goGPSgo
clearvars -except ini_settings_file use_gui; % exceptions for goGPSgo

global report

% if the plotting gets slower than usual, there might be problems with the
% Java garbage collector. In case, you can try to use the following
% command:
%
% java.lang.System.gc() %clear the Java garbage collector
%
% or:
%
% clear java

%cd(fileparts(which('goGPS')));
%pwd

% close all windows
% close all
fclose('all');
flag_init_out = false;

% clear the command prompt
%clc

% disable warnings
warning off;

% include all subdirectories
if (~isdeployed)
    addpath(genpath(pwd));
end

core = Core.getInstance();
core.showTextHeader();

log = Logger.getInstance();

% Pointer to the global settings:
gs = Go_State.getInstance();
state = gs.getCurrentSettings();
%settings_file = checkPath('..\data\project\default_PPP\config\settings.ini');
if exist('ini_settings_file', 'var')
    state.importIniFile(ini_settings_file);
end

%----------------------------------------------------------------------------------------------
% INTERFACE TYPE DEFINITION
%----------------------------------------------------------------------------------------------

if exist('use_gui', 'var')
    mode_user = use_gui;
else
    mode_user =   1; % user interface type
    % mode_user = 0 --> use text interface
    % mode_user = 1 --> use GUI
end

% Init output interfaces (singletons)
w_bar = Go_Wait_Bar.getInstance(100,'Welcome to goGPS', 0);  % 0 means text, 1 means GUI, 5 both

%if mode_user == 1
%    w_bar.setOutputType(1); % 0 means text, 1 means GUI, 5 both
%else
    w_bar.setOutputType(0); % 0 means text, 1 means GUI, 5 both
%end

% Kalman filter cannot be initialized when goGPS starts
kalman_initialized = false;

%----------------------------------------------------------------------------------------------
% INTERFACE STARTUP
%----------------------------------------------------------------------------------------------

global order o1 o2 o3 cutoff weights t nC
global cs_threshold
global iono_model tropo_model
global flag_outlier flag_outlier_OLOO SPP_threshold
global apriori_ZHD


% Set global variable for goGPS obj mode
clearvars -global goObj;
global goObj;
% For future development the flag goObs will guide a possible migration to the
% use of generic objects (such as goObservations) able to automatically manage
% multiple modes
goObj = 0;  % this variable is set in the interface.

if (mode_user == 1)
    
    % Now there's a unique interface for goGPS
    % to be compatible among various OSs the property "unit" of all the
    % elements must be set to "pixels"
    % (default unit is "character", but the size of a character is OS dependent)
    [ok_go] = gui_goGPS;
    if (~ok_go)
        return
    end
    
    [mode, mode_vinc, mode_data, mode_ref, flag_ms_pos, flag_ms, flag_ge, flag_cov, flag_NTRIP, flag_amb, ...
        flag_skyplot, flag_plotproc, flag_var_dyn_model, flag_stopGOstop, flag_SBAS, flag_IAR, ...
        filerootIN, ~, ~ , ~, ~, filename_ref, filename_pco, filename_blq, pos_M_man, protocol_idx, multi_antenna_rf, iono_model, tropo_model, fsep_char, ...
        flag_ocean, flag_outlier, flag_outlier_OLOO, flag_tropo, flag_tropo_gradient, frequencies, flag_SEID, processing_interval, obs_comb, flag_full_prepro, filename_sta, filename_met] = gs.settingsToGo(state);
    
else
    %----------------------------------------------------------------------------------------------
    % USER-DEFINED SETTINGS
    %----------------------------------------------------------------------------------------------
    
    [mode, mode_vinc, mode_data, mode_ref, flag_ms_pos, flag_ms, flag_ge, flag_cov, flag_NTRIP, flag_amb, ...
        flag_skyplot, flag_plotproc, flag_var_dyn_model, flag_stopGOstop, flag_SBAS, flag_IAR, ...
        filerootIN, ~, ~, ~, ~, filename_ref, filename_pco, filename_blq, pos_M_man, protocol_idx, multi_antenna_rf, iono_model, tropo_model, fsep_char, ...
        flag_ocean, flag_outlier, flag_outlier_OLOO, flag_tropo, flag_tropo_gradient, frequencies, flag_SEID, processing_interval, obs_comb, flag_full_prepro, filename_sta, filename_met] = gs.settingsToGo();
end

%-------------------------------------------------------------------------------------------
%% GO goGPS - here the computations start
%-------------------------------------------------------------------------------------------

log.newLine();
state.showTextMode();

gs.initProcessing(); % Set up / download observations and navigational files

cc = state.getConstellationCollector();

% start evaluating computation time
tic;

%-------------------------------------------------------------------------------------------
%% DETECT IF NAVIGATION FILE IS IN SP3 FORMAT
%-------------------------------------------------------------------------------------------
filename_nav = state.getFullNavEphPath(1);
if goGNSS.isPP(mode)
    fid = fopen(filename_nav);
    line1 = fgetl(fid); line2 = fgetl(fid);
    if (strcmp(line1(1),'#') && strcmp(line2(1:2),'##'))
        flag_SP3 = 1;
    else
        flag_SP3 = 0;
    end
    fclose(fid);
else
    flag_SP3 = 0;
end

%-------------------------------------------------------------------------------------------
%% STARTING BATCH
%-------------------------------------------------------------------------------------------

% Starting batch!!!
f_ref_rec = state.getReferencePath();
num_ref_rec = numel(f_ref_rec);
f_trg_rec = state.getTargetPath();
num_trg_rec = numel(f_trg_rec);
num_session = numel(f_trg_rec{1});
f_mst_rec = state.getMasterPath();
num_mst_rec = numel(f_mst_rec);

% get short name for File_Name_Processor
fnp = File_Name_Processor();

initial_mode = state.getMode();
if num_session > 1
    is_batch = true;
    w_bar.setOutputType(0);
    %log.setColorMode(0);
else
    is_batch = false;
end


% New implementation of goGPS
goGPS_new = false; % set to true to test the new goGPS implementation (VERY VERY PRELIMINAR))
if goGPS_new
    %% NEW goGPS
    state.showTextMode();
    
    sky = Core_Sky.getInstance();
    for s = 1 : num_session
        %-------------------------------------------------------------------------------------------
        % SESSION START
        %-------------------------------------------------------------------------------------------
        
        fprintf('\n--------------------------------------------------------------------------\n');
        log.addMessage(sprintf('Starting session %d of %d', s, num_session));
        fprintf('--------------------------------------------------------------------------\n');

        % Init sky
        fr = File_Rinex(f_trg_rec{1}{s},100);
        cur_date_start = fr.first_epoch.last();
        cur_date_stop = fr.last_epoch.first();
        sky.initSession(cur_date_start, cur_date_stop);

        % prepare reference time
        % processing time will start with the receiver with the last first epoch
        %          and it will stop  with the receiver with the first last epoch
        clear p_time;
        p_time_zero = round(cur_date_start.getMatlabTime() * 24)/24; % get the reference time
        p_time_start = -inf;
        p_time_stop = inf;
        p_rate = 1e-6;

        clear rec;  % handle to all the receivers
        clear mst;
        r = 0;
        for i = 1 : num_mst_rec
            log.newLine();
            log.addMessage(sprintf('Reading master %d of %d', i, num_mst_rec));
            fprintf('--------------------------------------------------------------------------\n\n');
            
            r = r + 1;
            mst(i) = Receiver(cc); %#ok<SAGROW>
            mst(i).loadRinex(f_mst_rec{i}{s});
            mst(i).correctTimeDesync();
            rec(r) = mst(i);
            
            % recompute the parameters for the ref_time estimation
            p_time_start = max(p_time_start,  round(rec(r).time.first.getRefTime(p_time_zero) * rec(r).time.getRate) / rec(r).time.getRate);
            p_time_stop = min(p_time_stop,  round(rec(r).time.last.getRefTime(p_time_zero) * rec(r).time.getRate) / rec(r).time.getRate);
            p_rate = lcm(round(p_rate * 1e6), round(rec(r).time.getRate * 1e6)) * 1e-6;
        end
        
        clear ref;
        for i = 1 : num_ref_rec
            log.newLine();
            log.addMessage(sprintf('Reading reference %d of %d', i, num_ref_rec));
            fprintf('--------------------------------------------------------------------------\n\n');
            
            r = r + 1;
            ref(i) = Receiver(cc); %#ok<SAGROW>
            ref(i).loadRinex(f_ref_rec{i}{s});
            ref(i).correctTimeDesync();
            rec(r) = ref(i);

            % recompute the parameters for the ref_time estimation
            p_time_start = max(p_time_start,  round(rec(r).time.first.getRefTime(p_time_zero) * rec(r).time.getRate) / rec(r).time.getRate);
            p_time_stop = min(p_time_stop,  round(rec(r).time.last.getRefTime(p_time_zero) * rec(r).time.getRate) / rec(r).time.getRate);
            p_rate = lcm(round(p_rate * 1e6), round(rec(r).time.getRate * 1e6)) * 1e-6;
        end
        
        clear trg;
        for i = 1 : num_trg_rec
            log.newLine();
            log.addMessage(sprintf('Reading target %d of %d', i, num_trg_rec));
            fprintf('--------------------------------------------------------------------------\n\n');
                        
            r = r + 1;
            trg(i) = Receiver(cc); %#ok<SAGROW>
            trg(i).loadRinex(f_trg_rec{i}{s});
            trg(i).static = state.kf_mode == 0; %#ok<SAGROW>
            trg(i).correctTimeDesync();
            rec(r) = trg(i);

            % recompute the parameters for the ref_time estimation
            % not that in principle I can have up to num_trg_rec ref_time
            % in case of multiple targets the reference times should be independent
            % so here I keep the temporary rt0 rt1 r_rate var
            % instead of ref_time_start, ref_time_stop, ref_rate
            pt0 = max(p_time_start, round(rec(r).time.first.getRefTime(p_time_zero) * rec(r).time.getRate) / rec(r).time.getRate);
            pt1 = min(p_time_stop, round(rec(r).time.last.getRefTime(p_time_zero) * rec(r).time.getRate) / rec(r).time.getRate);
            pr = lcm(lcm(round(p_rate * 1e6), round(rec(r).time.getRate * 1e6)), state.getProcessingRate() * 1e6) * 1e-6;
            
            p_time(i) = GPS_Time(p_time_zero, (pt0 : pr : pt1)); %#ok<SAGROW>
            p_time(i).toUnixTime();
            clear pt0 pt1 pr;
        end
        clear p_time_start p_time_stop p_rate r i;
        log.addMessage('Syncing times, computing reference time');

        keyboard;
    end
else

if is_batch
    %-------------------------------------------------------------------------------------------
    %% DISABLE FUNCTIONS NOT USED FOR BATCH PROCESSING
    %-------------------------------------------------------------------------------------------
    
    mode_vinc = 0;          % navigation mode
    mode_ref = 0;           % reference path mode
    flag_ms = 0;            % plot master station position --> no=0, yes=1
    flag_ge = 0;            % use google earth --> no=0, yes=1
    flag_cov = 0;           % plot error ellipse --> no=0, yes=1
    flag_NTRIP = 1;         % use NTRIP --> no=0, yes=1
    flag_amb = 0;           % plot ambiguities (only in post-processing)
    flag_skyplot = 0;       % draw skyplot and SNR graph (save CPU) --> no=0, yes=1
    flag_plotproc = 0;      % plot while processing
    flag_stopGOstop = 0;    % use a stop-go-stop procedure for direction estimation --> no=0, yes=1
    flag_var_dyn_model = 0; % variable dynamic model --> no=0, yes=1
    fsep_char = 'default';
end

active_sys = cc.getActive();
[constellations] = goGNSS.initConstellation(active_sys(1), active_sys(2), active_sys(3), active_sys(4), active_sys(5), active_sys(6));

nSatTot = cc.getNumSat();

%number of enabled constellations
n_sys = sum(cc.getActive);

for session = 1 : num_session
    %-------------------------------------------------------------------------------------------
    %% SESSION START
    %-------------------------------------------------------------------------------------------
    
    fprintf('\n--------------------------------------------------------------------------\n');
    log.addMessage(sprintf('Starting session %d of %d', session, num_session));
    fprintf('--------------------------------------------------------------------------\n\n');
    state.showTextMode();
    
    if is_batch
        clear X_KAL Xhat_t_t
    end
    
    %initialization of global variables/constants
    global_init;
    
    filename_R_obs = {};
    filename_M_obs = {};
    
    % Internally the current versionn of goGPS uses filename_R_obs for
    % multiple receivers and filename_M_obs for single, in SEID it is
    % necessary to assign as M the target and as R the reference receivers
    if state.isModeSEID()
        for r = 1 : num_ref_rec
            filename_R_obs{r} = f_ref_rec{r}{session}; %#ok<SAGROW>
        end
        filename_M_obs = f_trg_rec{1}{session};
        fr = File_Rinex([filename_R_obs(:); {filename_M_obs}],100);
    else
        for r = 1 : num_trg_rec
            filename_R_obs{r} = f_trg_rec{r}{session}; %#ok<SAGROW>
        end
        if state.isModeDD()
            for r = 1 : num_mst_rec
                filename_M_obs{r} = f_mst_rec{r}{min(session, numel(f_mst_rec{r}))}; %#ok<SAGROW>
            end
            fr = File_Rinex([filename_R_obs(:); filename_M_obs(:)],100);
        else
            fr = File_Rinex(filename_R_obs(:),100);
        end
    end
    
    if ~is_batch
        % close all the opened files
        fclose('all');
    end
    
    if ~fr.isValid()
        log.newLine();
        log.addWarning(sprintf('Skipping processing of "%s"', fnp.getRelDirPath(filename_R_obs{1})));
        log.newLine();
    else
        cur_date_start = fr.first_epoch.last();
        cur_date_stop = fr.last_epoch.first();
        
        % updating the file path of the output -> special key are now supported
        state.updateOutPath(cur_date_start, session);
        filerootOUT = state.getFullOutPath();
        
        % create a new directory when required
        dir_name = fileparts(filerootOUT);
        if ~(exist(dir_name, 'dir') == 7)
            mkdir(dir_name);
        end
        
        %-------------------------------------------------------------------------------------------
        %% REPORT INITIALIZATION
        %-------------------------------------------------------------------------------------------
        report = [];
        report.opt.write = 0;
        if state.isModePP() % report only if postprocessing
            report.opt.write = 1;
            
            if exist('is_batch','var') && is_batch
                report.opt.is_batch = 'YES';
            else
                report.opt.is_batch = 'NO';
            end
            
            report.opt.constellations = '';
            report.opt.sat_index=[];
            report.opt.sat_id={};
            flag_const='';
            if constellations.GPS.enabled == 1
                report.opt.constellations = [report.opt.constellations,flag_const,'GPS (G)'];
                flag_const=', ';
                report.opt.sat_index=[report.opt.sat_index,constellations.GPS.indexes];
                for i = 1:length(constellations.GPS.indexes)
                    report.opt.sat_id(constellations.GPS.indexes(i))=cellstr(sprintf(' G%02d',constellations.GPS.PRN(i)));
                end
            end
            if constellations.Galileo.enabled == 1
                report.opt.constellations = [report.opt.constellations,flag_const,'Galileo (E)'];
                flag_const=', ';
                report.opt.sat_index=[report.opt.sat_index,constellations.Galileo.indexes];
                for i = 1:length(constellations.Galileo.indexes)
                    report.opt.sat_id(constellations.Galileo.indexes(i))=cellstr(sprintf(' E%02d',constellations.Galileo.PRN(i)));
                end
            end
            if constellations.GLONASS.enabled == 1
                report.opt.constellations = [report.opt.constellations,flag_const,'GLONASS (R)'];
                flag_const=', ';
                report.opt.sat_index=[report.opt.sat_index,constellations.GLONASS.indexes];
                for i = 1:length(constellations.GLONASS.indexes)
                    report.opt.sat_id(constellations.GLONASS.indexes(i))=cellstr(sprintf(' R%02d',constellations.GLONASS.PRN(i)));
                end
            end
            if constellations.BeiDou.enabled == 1
                report.opt.constellations = [report.opt.constellations,flag_const,'BeiDou (C)'];
                flag_const=', ';
                report.opt.sat_index=[report.opt.sat_index,constellations.BeiDou.indexes];
                for i = 1:length(constellations.BeiDou.indexes)
                    report.opt.sat_id(constellations.BeiDou.indexes(i))=cellstr(sprintf(' C%02d',constellations.BeiDou.PRN(i)));
                end
            end
            if constellations.QZSS.enabled == 1
                report.opt.constellations = [report.opt.constellations,flag_const,'QZSS (J)'];
                flag_const=', ';
                report.opt.sat_index=[report.opt.sat_index,constellations.QZSS.indexes];
                for i = 1:length(constellations.QZSS.indexes)
                    report.opt.sat_id(constellations.QZSS.indexes(i))=cellstr(sprintf('J%03d',constellations.QZSS.PRN(i)));
                end
            end
            if constellations.SBAS.enabled == 1
                report.opt.constellations = [report.opt.constellations,flag_const,'SBAS (S)'];
                report.opt.sat_index=[report.opt.sat_index,constellations.SBAS.indexes];
                for i = 1:length(constellations.SBAS.indexes)
                    report.opt.sat_id(constellations.SBAS.indexes(i))=cellstr(sprintf(' S%02d',constellations.SBAS.PRN(i)));
                end
            end
            
            report.opt.outfolder = state.getFullOutPath();
            report.opt.mode = mode;
            report.opt.flag_SP3 = flag_SP3;
            report.opt.flag_ms_pos = flag_ms_pos;
            report.opt.flag_SBAS = flag_SBAS;
            report.opt.flag_IAR = flag_IAR;
            report.opt.cutoff = cutoff;
            report.opt.weights = weights;
            if exist('iniFile','var')
                report.inp.iniFile = fnp.getRelDirPath(iniFile, state.getHomeDir);
            end
            if goGNSS.isDD(mode)
                report.inp.filename_M_obs = fnp.getRelDirPath(filename_M_obs, state.getHomeDir);
            end
            report.inp.filename_R_obs = fnp.getRelDirPath(filename_R_obs, state.getHomeDir);
            report.inp.filename_nav = fnp.getRelDirPath(filename_nav, state.getHomeDir);
            report.inp.filename_pco = fnp.getRelDirPath(filename_pco, state.getHomeDir);
            report.inp.filename_blq = fnp.getRelDirPath(filename_blq, state.getHomeDir);
            if exist('filename_sta','var')
                report.inp.filename_sta = fnp.getRelDirPath(filename_sta, state.getHomeDir);
            end
            
            global sigmaq_cod1 sigmaq_cod2 sigmaq_ph sigmaq0_N min_nsat IAR_method flag_default_P0 flag_auto_mu mu P0 %#ok<TLEV>
            report.opt.sigmaq_cod1 = sigmaq_cod1;
            report.opt.sigmaq_cod2 = sigmaq_cod2;
            report.opt.sigmaq_ph = sigmaq_ph;
            report.opt.sigmaq0_N = sigmaq0_N;
            report.opt.min_nsat = min_nsat;
            report.opt.cs_threshold = cs_threshold;
            report.opt.IAR_method = IAR_method;
            report.opt.flag_default_P0 = flag_default_P0;
            report.opt.flag_auto_mu = flag_auto_mu;
            report.opt.mu = mu;
            report.opt.P0 = P0;
        end
        
        %----------------------------------------------------------------------------------------------
        %% FILE READING
        %----------------------------------------------------------------------------------------------
        read_files = true;
        while read_files
            read_files = false; % In case of SEID processing the read operation must be repeated twice
            
            if goGNSS.isPP(mode) % post-processing
                
                SP3 = [];
                
                % prepare the input for the load_RINEX_obs function
                filename_obs = multiple_RINEX_interface(filename_R_obs, filename_M_obs, mode);
                
                if state.isSA(mode) % absolute positioning
                    
                    %read observation RINEX file(s)
                    [pr1_R, ph1_R, pr2_R, ph2_R, dop1_R, dop2_R, snr1_R, snr2_R, ...
                        zero_time, time_GPS_diff, time_R_diff, week_R, date_R, pos_R, interval, antoff_R, antmod_R, codeC1_R, marker_R] = ...
                        load_RINEX_obs(filename_obs, state.getConstellationCollector(), processing_interval);
                    el_r = nan(size(pr1_R));
                    az_r = nan(size(pr1_R));
                    
                    time_GPS = zero_time + time_GPS_diff;
                    time_R = zero_time + time_R_diff;
                    
                    %read navigation RINEX file(s)
                    [Eph, iono, flag_return] = load_RINEX_nav(state.getEphFileName(cur_date_start, cur_date_stop), state.getConstellationCollector(), flag_SP3, iono_model, time_GPS);
                    if (flag_return)
                        return
                    end
                    
                    if (~exist('time_GPS','var') || ~any(isfinite(time_GPS)) || isempty(time_GPS))
                        fprintf('... WARNING: either there are no observations available for processing, or some epoch is not valid.\n');
                        return
                    end
                    
                    %read stations coordinates file
                    if (~exist('pos_R_crd','var') || ~any(serialize(pos_R_crd)))
                        [pos_R_crd, flag_XR, pos_M_crd, flag_XM] = load_CRD(filename_sta, marker_R, []);
                    end
                    
                    %read receiver antenna phase center offset (PCO) and variation (PCV)
                    antenna_PCV = read_antenna_PCV(filename_pco, antmod_R);
                    
                    %read satellite antenna phase center offset (NOTE: reading only L1 offset for now)
                    antmod_S = cc.getAntennaId();
                    antenna_PCV_S = read_antenna_PCV(filename_pco, antmod_S, date_R);
                    
                    % write report file    %%-> must be extented in MR case
                    if report.opt.write == 1
                        % extract quality parameters for report
                        for i = 1:size(pr1_R,3)
                            if (antenna_PCV(i).n_frequency ~= 0)
                                report.obs.pcv_yn(i)=cellstr('YES');
                            else
                                report.obs.pcv_yn(i)=cellstr('NO');
                            end
                            report.obs.antname_R(i) = cellstr(antenna_PCV(i).name);
                            report.obs.antoff_R(i,:) = antoff_R(:,:,i);
                            
                            % set ROVER initial coordinates
                            if (exist('pos_R_crd','var') && any(pos_R_crd(:)))
                                log.newLine();
                                log.addMessage('Rover apriori position set from coordinate file:');
                                log.addMessage(sprintf('     X = %12.4f m\n     Y = %12.4f m\n     Z = %12.4f m', pos_R_crd(1,1), pos_R_crd(2,1), pos_R_crd(3,1)));
                                if report.opt.write == 1
                                    report.obs.coord_R=sprintf('%-30s  %13.4f %13.4f %13.4f  approx from coordinate file', char(report.obs.filename(1)), pos_R_crd(1,1), pos_R_crd(2,1), pos_R_crd(3,1));
                                end
                                pos_R = pos_R_crd;
                            else
                                log.newLine();
                                log.addMessage('Rover apriori position set from RINEX:');
                                log.addMessage(sprintf('     X = %12.4f m\n     Y = %12.4f m\n     Z = %12.4f m', pos_R(1,1), pos_R(2,1), pos_R(3,1)));
                                if report.opt.write == 1
                                    if any(pos_R)
                                        report.obs.coord_R=sprintf('%-30s  %13.4f %13.4f %13.4f  approx from RINEX', char(report.obs.filename(1)), pos_R(1,1), pos_R(2,1), pos_R(3,1));
                                    else
                                        report.obs.coord_R=sprintf('%-30s  apriori coordinates not available          estimated from observations  ', char(report.obs.filename(1)));
                                    end
                                end
                            end
                            log.newLine();
                        end
                    end
                    
                    if (flag_SP3)
                        
                        %----------------------------------------------------------------------------------------------
                        %% LOAD SP3 DATA
                        %----------------------------------------------------------------------------------------------
                        
                        SP3 = load_SP3(state.getEphFileName(cur_date_start, cur_date_stop), state.getClkFileName(cur_date_start, cur_date_stop), time_GPS, week_R, constellations, cur_date_start, cur_date_stop);
                        
                        %store satellite antenna PCO/PCV and satellite type
                        SP3.antPCO = zeros(1,3,size(antenna_PCV_S,2));
                        SP3.satType = cell(1,size(antenna_PCV_S,2));
                        for sat = 1 : size(antenna_PCV_S,2)
                            if (antenna_PCV_S(sat).n_frequency ~= 0)
                                SP3.antPCO(:,:,sat) = antenna_PCV_S(sat).offset(:,:,1);
                                SP3.satType{1,sat} = antenna_PCV_S(sat).sat_type;
                            else
                                SP3.avail(sat) = 0;
                            end
                        end
                        
                        %compute sun and moon position
                        fprintf('Computing Sun and Moon position...');
                        [X_sun, X_moon] = sun_moon_pos(datevec(gps2utc(datenum(date_R))));
                        fprintf(' done\n');
                        
                        %store the position of Sun and Moon
                        SP3.t_sun  = time_GPS;
                        SP3.X_sun  = X_sun;
                        SP3.X_moon = X_moon;
                        
                        %----------------------------------------------------------------------------------------------
                        %% LOAD ERP DATA (EARTH ROTATION/ORIENTATION PARAMETERS)
                        %----------------------------------------------------------------------------------------------
                        
                        SP3.ERP = load_ERP(state.erp_full_name, time_GPS);
                        
                        %----------------------------------------------------------------------------------------------
                        %% LOAD DCB DATA (DIFFERENTIAL CODE BIASES)
                        %----------------------------------------------------------------------------------------------
                        
                        %NOTE: if not using SP3 ephemeris or if DCB files are not available, the
                        %      'SP3.DCB' structure will be initialized to zero/empty arrays and it will not
                        %      have any effect on the positioning
                        
                        %if (~strcmp(obs_comb, 'IONO_FREE'))
                        %try first to read already available DCB files
                        DCB = load_dcb(state.dcb_dir, week_R, time_R, codeC1_R, constellations);
                        
                        %if DCB files are not available or not sufficient, try to download them
                        if ((~any(DCB.P1C1.value(:)) || ~any(DCB.P1P2.value(:))) && constellations.GPS.enabled)
                            
                            %download
                            [file_dcb, compressed] = download_dcb([week_R(1) week_R(end)], [time_R(1) time_R(end)]);
                            
                            if (compressed)
                                return
                            end
                            
                            %try again to read DCB files
                            DCB = load_dcb(state.dcb_dir, week_R, time_R, codeC1_R, constellations);
                        end
                        
                        SP3.DCB = DCB;
                        %else
                        %SP3.DCB = [];
                        %end
                    end
                    
                    %----------------------------------------------------------------------------------------------
                    %% LOAD CRX DATA (SATELLITE PROBLEMS: MANEUVERS OR BAD OBSERVATION INTERVALS)
                    %----------------------------------------------------------------------------------------------
                    
                    %try first to read already available CRX files
                    [CRX, found]  = load_crx(state.crx_dir, week_R, time_GPS, state.getConstellationCollector());
                    %if CRX files are not available or not sufficient, try to download them
                    if (~found)
                        %download
                        file_crx = download_crx([week_R(1) week_R(end)], [time_GPS(1) time_GPS(end)]);
                        
                        %try again to read CRX files
                        [CRX, found] = load_crx(state.crx_dir, week_R, time_GPS, state.getConstellationCollector());
                    end
                    
                    %retrieve multi-constellation wavelengths
                    lambda = goGNSS.getGNSSWavelengths(Eph, SP3, nSatTot);
                    
                    %exclude for which lambda could not be computed
                    delsat = ~any(lambda,2);
                    pr1_R(delsat,:,:) = 0;
                    pr2_R(delsat,:,:) = 0;
                    ph1_R(delsat,:,:) = 0;
                    ph2_R(delsat,:,:) = 0;
                    dop1_R(delsat,:,:) = 0;
                    dop2_R(delsat,:,:) = 0;
                    snr_R(delsat,:,:) = 0; %#ok<SAGROW>
                    
                    pos_R_new    = zeros(3, length(time_GPS), size(time_R,3));
                    dtR          = zeros(length(time_GPS), 1, size(time_R,3));
                    dtRdot       = zeros(length(time_GPS), 1, size(time_R,3));
                    bad_sats_R   = zeros(nSatTot, 1, size(time_R,3));
                    status_obs_R = zeros(nSatTot, length(time_GPS), size(time_R,3));
                    bad_epochs_R = NaN(length(time_GPS), 1, size(time_R,3));
                    var_SPP_R    = NaN(length(time_GPS), 3, size(time_R,3));
                    var_dtR      = NaN(length(time_GPS), 1, size(time_R,3));
                    
                    if (~exist('pos_R_crd','var') || ~any(pos_R_crd(:)))
                        if any(pos_R)
                            flag_XR = 1;
                        else
                            flag_XR = 0;
                        end
                    end
                    
                    report.errors.few_epochs = 0;
                    report.opt.min_epoch = 0;
                    
                    if exist('min_epoch','var')
                        report.opt.min_epoch = min_epoch;
                        if size(time_R,1) < min_epoch
                            fprintf('\nERROR! The number of available epochs is lower than the minimum. The processing will not be performed.\n');
                            % write report
                            report.errors.few_epochs = 1;
                            report_generator(report);
                            return
                        end
                    end
                    
                    %if SBAS corrections are requested
                    if (flag_SBAS)
                        
                        %----------------------------------------------------------------------------------------------
                        %% LOAD SBAS DATA (EGNOS EMS FILES)
                        %----------------------------------------------------------------------------------------------
                        
                        %NOTE: if SBAS corrections are not requested by the user or not available, the
                        %      'sbas' structure will be initialized to zero/empty arrays and it will not
                        %      have any effect on the positioning
                        
                        %try first to read .ems files already available
                        [sbas] = load_ems(state.ems_dir, week_R, time_R);
                        
                        %if .ems files are not available or not sufficient, try to download them
                        if (isempty(sbas))
                            
                            %EGNOS PRNs
                            prn = [120, 124, 126];
                            
                            %download
                            for p = 1 : length(prn)
                                [file_ems] = download_ems(prn(p), [week_R(1) week_R(end)], [time_R(1) time_R(end)]);
                                if (~isempty(file_ems))
                                    break
                                end
                            end
                            
                            %try again to read .ems files
                            [sbas] = load_ems(state.ems_dir, week_R, time_R);
                        end
                        
                        %check if the survey is within the EMS grids
                        if (~isempty(sbas))
                            [ems_data_available] = check_ems_extents(time_R, pr1_R, snr1_R, nSatTot, Eph, SP3, iono, sbas, lambda, 1);
                        end
                    end
                    
                    %if SBAS corrections are not requested or not available
                    if (~flag_SBAS || isempty(sbas) || ~ems_data_available)
                        
                        %initialization
                        sbas = [];
                        
                        %if SBAS corrections are requested but not available
                        if (flag_SBAS && isempty(sbas))
                            fprintf('Switching back to standard (not SBAS-corrected) processing.\n')
                        end
                    end
                    
                    %if ocean loading correction is requested
                    if (flag_ocean)
                        ol_disp = load_BLQ(filename_blq, marker_R);
                    end
                    
                    %time adjustments (to account for sub-integer approximations in MATLAB - thanks to radiolabs.it for pointing this out!)
                    if (flag_SP3)
                        SP3.time    = SP3.time - zero_time;
                        SP3.time_hr = SP3.time_hr - zero_time;
                        SP3.t_sun   = SP3.t_sun - zero_time;
                        if ~isempty(SP3.ERP)
                            SP3.ERP.t   = SP3.ERP.t - zero_time;
                        end
                    end
                    Eph(32,:) = Eph(32,:) - zero_time;
                    Eph(33,:) = Eph(33,:) - zero_time;
                    
                    for f = 1 : size(time_R,3)
                        
                        %apply P1C1 DCBs if needed
                        if (flag_SP3 && ~isempty(SP3.DCB) && any(codeC1_R(:)))
                            avail_sat = any(lambda,2);
                            pr1_R(avail_sat,:,f) = pr1_R(avail_sat,:,f) + SP3.DCB.P1C1.value(avail_sat,ones(size(pr1_R(:,:,f),2),1))*1e-9*goGNSS.V_LIGHT.*codeC1_R(avail_sat,:,f);
                        end
                        
                        %pre-processing
                        log.addMarkedMessage(['Pre-processing rover observations (file ' filename_obs{f} ')...']);
                        log.newLine();
                        w_bar.setBarLen(length(time_GPS_diff));
                        w_bar.createNewBar('Pre-processing rover...');
                        
                        [time_R_diff(:,1,f), pr1_R(:,:,f), ph1_R(:,:,f), pr2_R(:,:,f), ph2_R(:,:,f), pos_R_new(:,:,f), dtR(:,1,f), dtRdot(:,1,f), el_r(:,:,f), az_r(:,:,f), bad_sats_R(:,1,f), bad_epochs_R(:,1,f), var_dtR(:,1,f), var_SPP_R(:,:,f), status_obs_R(:,:,f), status_cs, eclipsed, ISBs, var_ISBs] = pre_processing(time_GPS_diff, time_R_diff(:,1,f), pos_R, pr1_R(:,:,f), ph1_R(:,:,f), pr2_R(:,:,f), ph2_R(:,:,f), dop1_R(:,:,f), dop2_R(:,:,f), snr1_R(:,:,f), Eph, SP3, iono, lambda, frequencies, obs_comb, nSatTot, w_bar, flag_XR, sbas, cc, flag_full_prepro, order);
                        % the time is now changed
                        time_r(f) = GPS_Time(zero_time/86400 + GPS_Time.GPS_ZERO, time_R_diff(:,:,f)); %#ok<SAGROW>
                        date_R(:,:,f) = datevec(time_r(f).getMatlabTime);
                        
                        if report.opt.write == 1
                            report.prep.spp_threshold = SPP_threshold;
                            report.prep.flag_R = flag_XR;
                            report.prep.tot_epoch_R(f)=size(pr1_R(:,:,f),2);
                            report.prep.proc_epoch_R(f)=length(bad_epochs_R(isfinite(bad_epochs_R(:,1,f)),1,f));
                            report.prep.bad_epoch_R(f)=sum(bad_epochs_R(isfinite(bad_epochs_R(:,1,f)),1,f)==1);
                            if (~isempty(var_SPP_R(isfinite(var_SPP_R(:,1,f)),1,f)))
                                report.prep.max_varSPP_R(f)=max(var_SPP_R(isfinite(var_SPP_R(:,1,f)),1,f))^0.5;
                            else
                                report.prep.max_varSPP_R(f) = NaN;
                            end
                            report.prep.varSPP_R(f)=(sum(var_SPP_R(isfinite(var_SPP_R(:,2,f)),2,f))/sum(var_SPP_R(isfinite(var_SPP_R(:,2,f)),3,f)))^.5;
                            report.prep.tot_obs_R(f)=length(find(isfinite(status_obs_R(:,:,f))));
                            report.prep.obs_outlier_R(f)=length(find(status_obs_R(:,:,f)==-1));
                            report.prep.obs_used_R(f)=length(find(status_obs_R(:,:,f)==1));
                            report.prep.obs_undercutoff_R(f)=length(find(status_obs_R(:,:,f)==0));
                            report.prep.obs_stat_R(:,:,f)=[sum(status_obs_R(:,:,f)==0,2), sum(status_obs_R(:,:,f)==1,2), sum(status_obs_R(:,:,f)==-1,2)]; % [#under_cutoff, #used, #outlier] grouped by satellite
                            report.prep.CS_R{f}=status_cs;
                        end
                        
                        w_bar.close();
                    end
                    
                    %global residuals_fixed residuals_float outliers s02_ls %#ok<TLEV>
                    %residuals_fixed=NaN(2*length(n_freq)*nSatTot,1);
                    %residuals_float=NaN(2*length(n_freq)*nSatTot,1);
                    %outliers=zeros(2*length(n_freq)*nSatTot,1);
                    %s02_ls=NaN(length(time_GPS),1);
                    
                else %relative positioning
                    
                    %-------------------------------------------------------------------------------------------
                    %% READING MR
                    %-------------------------------------------------------------------------------------------
                    %read observation RINEX file(s)
                    [pr1_RM, ph1_RM, pr2_RM, ph2_RM, dop1_RM, dop2_RM, snr1_RM, snr2_RM, ...
                        zero_time, time_GPS_diff, time_RM_diff, week_RM, date_RM, pos_RM, interval, antoff_RM, antmod_RM, codeC1_RM, marker_RM] = ...
                        load_RINEX_obs(filename_obs, state.getConstellationCollector(), processing_interval);
                    
                    time_GPS = zero_time + time_GPS_diff;
                    time_RM = zero_time + time_RM_diff;
                    
                    [Eph, iono, flag_return] = load_RINEX_nav(state.getEphFileName(cur_date_start, cur_date_stop), state.getConstellationCollector(), flag_SP3, iono_model, time_GPS);
                    if (flag_return)
                        return
                    end
                    
                    if (~exist('time_GPS','var') || ~any(isfinite(time_GPS)) || isempty(time_GPS))
                        log.addWarning(' Either there are no observations available for processing, or some epoch is not valid.');
                        return
                    end
                    
                    pr1_R = pr1_RM(:,:,1:end-1); pr1_M = pr1_RM(:,:,end);
                    ph1_R = ph1_RM(:,:,1:end-1); ph1_M = ph1_RM(:,:,end);
                    pr2_R = pr2_RM(:,:,1:end-1); pr2_M = pr2_RM(:,:,end);
                    ph2_R = ph2_RM(:,:,1:end-1); ph2_M = ph2_RM(:,:,end);
                    dop1_R = dop1_RM(:,:,1:end-1); dop1_M = dop1_RM(:,:,end);
                    dop2_R = dop2_RM(:,:,1:end-1); dop2_M = dop2_RM(:,:,end);
                    snr1_R = snr1_RM(:,:,1:end-1); snr1_M = snr1_RM(:,:,end);
                    snr2_R = snr2_RM(:,:,1:end-1); snr2_M = snr2_RM(:,:,end);
                    time_R = time_RM(:,1,1:end-1); time_M = time_RM(:,1,end);
                    time_R_diff = time_RM_diff(:,1,1:end-1); time_M_diff = time_RM_diff(:,1,end);
                    week_R = week_RM(:,1,1:end-1); week_M = week_RM(:,1,end);
                    date_R = date_RM(:,:,1:end-1); date_M = date_RM(:,:,end);
                    pos_R = pos_RM(:,1,1:end-1); pos_M = pos_RM(:,1,end);
                    antoff_R = antoff_RM(:,1,1:end-1); antoff_M = antoff_RM(:,1,end);
                    codeC1_R = codeC1_RM(:,:,1:end-1); codeC1_M = codeC1_RM(:,:,end);
                    marker_R = marker_RM(:,1,1:end-1); marker_M = marker_RM(:,1,end);
                    el_r = nan(size(pr1_R)); el_m = nan(size(pr1_M));
                    az_r = nan(size(pr1_R)); az_m = nan(size(pr1_M));

                    % read stations coordinates file
                    if (~exist('pos_R_crd','var') || ~any(pos_R_crd(:)) || ~exist('pos_M_crd','var') || ~any(pos_M_crd(:)))
                        [pos_R_crd, flag_XR, pos_M_crd, flag_XM] = load_CRD(filename_sta, marker_R, marker_M);
                    end
                    
                    % read receiver antenna phase center offset
                    antenna_PCV = read_antenna_PCV(filename_pco, antmod_RM);
                    
                    % read satellite antenna phase center offset
                    antmod_S = cc.getAntennaId();
                    antenna_PCV_S = read_antenna_PCV(filename_pco, antmod_S, date_M);
                    
                    if report.opt.write == 1
                        % extract quality parameters for report
                        for i = 1:1 %size(pr1_R,3) that is because currenty MR is not supported
                            if (antenna_PCV(i).n_frequency ~= 0)
                                report.obs.pcv_yn(i)=cellstr('YES');
                            else
                                report.obs.pcv_yn(i)=cellstr('NO');
                            end
                            report.obs.antname_R(i) = cellstr(antenna_PCV(i).name);
                            report.obs.antoff_R(i,:) = antoff_R(:,:,i);
                        end
                        
                        report.obs.antname_M = antenna_PCV(end).name;
                        report.obs.antoff_M = antoff_M;
                        if (antenna_PCV(end).n_frequency ~= 0)
                            report.obs.pcv_yn(i+1)=cellstr('YES');
                        else
                            report.obs.pcv_yn(i+1)=cellstr('NO');
                        end
                    end
                    
                    % set ROVER initial coordinates
                    if any(pos_R_crd)
                        log.newLine();
                        log.addMessage('Rover apriori position set from coordinate file:');
                        log.addMessage(sprintf('     X = %12.4f m\n     Y = %12.4f m\n     Z = %12.4f m', pos_R_crd(1,1), pos_R_crd(2,1), pos_R_crd(3,1)));
                        if report.opt.write == 1
                            if (flag_XR ~= 2)
                                report.obs.coord_R=sprintf('%-30s  %13.4f %13.4f %13.4f  approx from coordinate file', char(report.obs.filename(1)), pos_R_crd(1,1), pos_R_crd(2,1), pos_R_crd(3,1));
                            else
                                report.obs.coord_R=sprintf('%-30s  %13.4f %13.4f %13.4f  fixed from coordinate file', char(report.obs.filename(1)), pos_R_crd(1,1), pos_R_crd(2,1), pos_R_crd(3,1));
                            end
                        end
                        pos_R = pos_R_crd;
                    else
                        log.newLine();
                        log.addMessage('Rover apriori position set from RINEX:');
                        log.addMessage(sprintf('     X = %12.4f m\n     Y = %12.4f m\n     Z = %12.4f m', pos_R(1,1), pos_R(2,1), pos_R(3,1)));
                        if report.opt.write == 1
                            if any(pos_R)
                                report.obs.coord_R=sprintf('%-30s  %13.4f %13.4f %13.4f  approx from RINEX', char(report.obs.filename(1)), pos_R(1,1), pos_R(2,1), pos_R(3,1));
                            else
                                report.obs.coord_R=sprintf('%-30s  apriori coordinates not available          estimated from observations  ', char(report.obs.filename(1)));
                            end
                        end
                    end
                    
                    % set MASTER initial coordinates
                    if (flag_ms_pos)
                        if any(pos_M_crd) % master position read from coordinate file
                            log.newLine();
                            log.addMessage('Master position fixed from coordinate file:');
                            log.addMessage(sprintf('     X = %12.4f m\n     Y = %12.4f m\n     Z = %12.4f m', pos_M_crd(1,1), pos_M_crd(2,1), pos_M_crd(3,1)));
                            if report.opt.write == 1
                                report.obs.coord_M=sprintf('%-30s  %13.4f %13.4f %13.4f  fixed from coordinate file', char(report.obs.filename(end)), pos_M_crd(1,1), pos_M_crd(2,1), pos_M_crd(3,1));
                            end
                            pos_M = pos_M_crd;
                        else % master position read from RINEX header
                            log.newLine();
                            log.addMessage('Master position fixed from RINEX:');
                            log.addMessage(sprintf('     X = %12.4f m\n     Y = %12.4f m\n     Z = %12.4f m', pos_M(1,1), pos_M(2,1), pos_M(3,1)));
                            if report.opt.write == 1
                                report.obs.coord_M=sprintf('%-30s  %13.4f %13.4f %13.4f  fixed from RINEX', char(report.obs.filename(end)), pos_M(1,1), pos_M(2,1), pos_M(3,1));
                            end
                        end
                    else % no valid pos_M_man found, so force positiong read from RINEX header
                        log.addWarning('MASTER coordinates forced fixed from RINEX');
                        log.addMessage(sprintf('                 X = %12.4f m\n                 Y = %12.4f m\n                 Z = %12.4f m', pos_M(1,1), pos_M(2,1), pos_M(3,1)));
                        if report.opt.write == 1
                            report.obs.coord_M=sprintf('%-30s  %13.4f %13.4f %13.4f  forced fixed from RINEX', char(report.obs.filename(end)), pos_M(1,1), pos_M(2,1), pos_M(3,1));
                        end
                    end
                    
                    log.newLine();
                    % apply antenna offset over the marker to master coordinates
                    pos_M = local2globalPos(antoff_M, pos_M);
                    
                    % apply antenna offset over the marker to rover apriori coordinates
                    if (any(pos_R))
                        for i = 1 : size(pr1_R,3)
                            pos_R(:,:,i) = local2globalPos(antoff_R(:,:,i), pos_R(:,:,i));
                        end
                    end
                    
                    if (flag_SP3)
                        %display message
                        SP3 = load_SP3(state.getEphFileName(cur_date_start, cur_date_stop), state.getClkFileName(cur_date_start, cur_date_stop), time_GPS, week_M, constellations, cur_date_start, cur_date_stop);
                        
                        %store satellite antenna PCO/PCV and satellite type
                        SP3.antPCO = zeros(1,3,size(antenna_PCV_S,2));
                        SP3.satType = cell(1,size(antenna_PCV_S,2));
                        for sat = 1 : size(antenna_PCV_S,2)
                            if (antenna_PCV_S(sat).n_frequency ~= 0)
                                SP3.antPCO(:,:,sat) = antenna_PCV_S(sat).offset(:,:,1);
                                SP3.satType{1,sat} = antenna_PCV_S(sat).type;
                            else
                                SP3.avail(sat) = 0;
                            end
                        end
                        
                        %compute sun and moon position
                        log.addMessage('Computing Sun and Moon position...');
                        [X_sun, X_moon] = sun_moon_pos(datevec(gps2utc(datenum(date_M))));
                        log.addStatusOk();
                        log.newLine();
                        
                        %store the position of Sun and Moon
                        SP3.t_sun  = time_GPS;
                        SP3.X_sun  = X_sun;
                        SP3.X_moon = X_moon;
                        
                        %----------------------------------------------------------------------------------------------
                        %% LOAD ERP DATA (EARTH ROTATION/ORIENTATION PARAMETERS)
                        %----------------------------------------------------------------------------------------------
                        
                        ERP = load_ERP(state.erp_full_name, time_GPS);
                        SP3.ERP = ERP;
                        
                        %----------------------------------------------------------------------------------------------
                        %% LOAD DCB DATA (DIFFERENTIAL CODE BIASES)
                        %----------------------------------------------------------------------------------------------
                        
                        %NOTE: if not using SP3 ephemeris or if DCB files are not available, the
                        %      'SP3.DCB' structure will be initialized to zero/empty arrays and it will not
                        %      have any effect on the positioning
                        
                        %if (~strcmp(obs_comb, 'IONO_FREE'))
                        %try first to read already available DCB files
                        DCB = load_dcb(state.dcb_dir, week_M, time_M, or(codeC1_R,codeC1_M(:,:,ones(1,size(codeC1_R,3)))), constellations);
                        
                        %if DCB files are not available or not sufficient, try to download them
                        if ((~any(DCB.P1C1.value(:)) || ~any(DCB.P1P2.value(:))) && constellations.GPS.enabled)
                            
                            %download
                            [file_dcb, compressed] = download_dcb([week_M(1) week_M(end)], [time_M(1) time_M(end)]);
                            
                            if (compressed)
                                return
                            end
                            
                            %try again to read DCB files
                            DCB = load_dcb(state.dcb_dir, week_M, time_M, or(codeC1_R,codeC1_M(:,:,ones(1,size(codeC1_R,3)))), constellations);
                        end
                        
                        SP3.DCB = DCB;
                        %else
                        %SP3.DCB = [];
                        %end
                    end
                    
                    %----------------------------------------------------------------------------------------------
                    %% LOAD CRX DATA (SATELLITE PROBLEMS: MANEUVERS OR BAD OBSERVATION INTERVALS)
                    %----------------------------------------------------------------------------------------------
                    
                    %try first to read already available CRX files
                    [CRX, found] = load_crx(state.crx_dir, week_M, time_GPS, state.getConstellationCollector());
                    %if CRX files are not available or not sufficient, try to download them
                    if (~found)
                        %download
                        file_crx = download_crx([week_M(1) week_M(end)], [time_GPS(1) time_GPS(end)]);
                        
                        %try again to read CRX files
                        [CRX, found] = load_crx(state.crx_dir, week_M, time_GPS, state.getConstellationCollector());
                    end
                    
                    %retrieve multi-constellation wavelengths
                    lambda = goGNSS.getGNSSWavelengths(Eph, SP3, nSatTot);
                    
                    %exclude for which lambda could not be computed
                    delsat = ~any(lambda,2);
                    pr1_R(delsat,:,:) = 0;
                    codeC1_R(delsat,:) = 0;
                    pr2_R(delsat,:,:) = 0;
                    ph1_R(delsat,:,:) = 0;
                    ph2_R(delsat,:,:) = 0;
                    dop1_R(delsat,:,:) = 0;
                    dop2_R(delsat,:,:) = 0;
                    snr_R(delsat,:,:) = 0; %#ok<SAGROW>
                    pr1_M(delsat,:) = 0;
                    codeC1_M(delsat,:) = 0;
                    pr2_M(delsat,:) = 0;
                    ph1_M(delsat,:) = 0;
                    ph2_M(delsat,:) = 0;
                    dop1_M(delsat,:) = 0;
                    dop2_M(delsat,:) = 0;
                    snr_M(delsat,:) = 0; %#ok<SAGROW>                    
                    
                    pos_R_new    = zeros(3, length(time_GPS), size(time_R,3));
                    dtR          = zeros(length(time_GPS), 1, size(time_R,3));
                    dtRdot       = zeros(length(time_GPS), 1, size(time_R,3));
                    bad_sats_R   = zeros(nSatTot, 1, size(time_R,3));
                    status_obs_R = zeros(nSatTot, length(time_GPS), size(time_R,3));
                    bad_epochs_R = NaN(length(time_GPS), 1, size(time_R,3));
                    var_SPP_R    = NaN(length(time_GPS), 3, size(time_R,3));
                    var_dtR      = NaN(length(time_GPS), 1, size(time_R,3));
                    
                    report.errors.few_epochs = 0;
                    report.opt.min_epoch = 0;
                    
                    if exist('min_epoch','var')
                        report.opt.min_epoch = min_epoch;
                        if size(time_R,1) < min_epoch
                            log.addError('The number of available epochs is lower than the minimum. The processing will not be performed.');
                            % write report
                            report.errors.few_epochs = 1;
                            report_generator(report);
                            return
                        end
                    end
                    
                    %if SBAS corrections are requested
                    if (flag_SBAS)
                        
                        %----------------------------------------------------------------------------------------------
                        %% LOAD SBAS DATA (EGNOS EMS FILES)
                        %----------------------------------------------------------------------------------------------
                        
                        %NOTE: if SBAS corrections are not requested by the user or not available, the
                        %      'sbas' structure will be initialized to zero/empty arrays and it will not
                        %      have any effect on the positioning
                        
                        %try first to read .ems files already available
                        [sbas] = load_ems(state.ems_dir, week_M, time_M);
                        
                        %if .ems files are not available or not sufficient, try to download them
                        if (isempty(sbas))
                            
                            %EGNOS PRNs
                            prn = [120, 124, 126];
                            
                            %download
                            for p = 1 : length(prn)
                                [file_ems] = download_ems(prn(p), [week_M(1) week_M(end)], [time_M(1) time_M(end)]);
                                if (~isempty(file_ems))
                                    break
                                end
                            end
                            
                            %try again to read .ems files
                            [sbas] = load_ems(state.ems_dir, week_M, time_M);
                        end
                        
                        %check if the survey is within the EMS grids
                        if (~isempty(sbas))
                            [ems_data_available] = check_ems_extents(time_M, pr1_M, snr1_M, nSatTot, Eph, SP3, iono, sbas, lambda, 1);
                        end
                    end
                    
                    %if SBAS corrections are not requested or not available
                    if (~flag_SBAS || isempty(sbas) || ~ems_data_available)
                        
                        %initialization
                        sbas = [];
                        
                        %if SBAS corrections are requested but not available
                        if (flag_SBAS && isempty(sbas))
                            log.addMessage('Switching back to standard (not SBAS-corrected) processing.')
                        end
                    end
                    
                    %if ocean loading correction is requested
                    if (flag_ocean)
                        ol_disp = load_BLQ(filename_blq, marker_RM);
                    end
                    
                    %time adjustments (to account for sub-integer approximations in MATLAB - thanks to radiolabs.it for pointing this out!)
                    if (flag_SP3)
                        SP3.time    = SP3.time - zero_time;
                        SP3.time_hr = SP3.time_hr - zero_time;
                        SP3.t_sun   = SP3.t_sun - zero_time;
                        if ~isempty(SP3.ERP)
                            SP3.ERP.t   = SP3.ERP.t - zero_time;
                        end
                    end
                    Eph(32,:) = Eph(32,:) - zero_time;
                    Eph(33,:) = Eph(33,:) - zero_time;
                    
                    log.addMessage(['Pre-processing master observations (file ' filename_obs{end} ')...']);
                    w_bar.setBarLen(length(time_GPS));
                    w_bar.createNewBar('Pre-processing master...');
                    
                    %apply P1C1 DCBs if needed
                    if (flag_SP3 && ~isempty(DCB) && any(codeC1_M(:)))
                        avail_sat = any(lambda,2);
                        pr1_M(avail_sat,:) = pr1_M(avail_sat,:) + SP3.DCB.P1C1.value(avail_sat,ones(size(pr1_M,2),1))*1e-9*goGNSS.V_LIGHT.*codeC1_M(avail_sat,:);
                    end
                    
                     if (~flag_SEID)
                        flag_XM_prep = iif(flag_XM == 0, 1, flag_XM);
                        if state.isModeBlock()
                            pp = Core_Pre_Processing(state, Eph, SP3);
                            [time_M_diff, pr1_M, ph1_M, pr2_M, ph2_M, ~, dtM, dtMdot, el_m, az_m, bad_sats_M, bad_epochs_M, var_dtM, var_SPP_M, status_obs_M, status_cs, eclipsed, ISBs, var_ISBs] = pp.execute(time_GPS_diff, time_M_diff, pos_M, pr1_M, ph1_M, pr2_M, ph2_M, dop1_M, dop2_M, snr1_M, iono, lambda, frequencies, obs_comb, nSatTot, w_bar, flag_XM_prep, sbas, flag_full_prepro, order);
                        else
                            [time_M_diff, pr1_M, ph1_M, pr2_M, ph2_M, ~, dtM, dtMdot, el_m, az_m, bad_sats_M, bad_epochs_M, var_dtM, var_SPP_M, status_obs_M, status_cs, eclipsed, ISBs, var_ISBs] = pre_processing(time_GPS_diff, time_M_diff, pos_M, pr1_M, ph1_M, pr2_M, ph2_M, dop1_M, dop2_M, snr1_M, Eph, SP3, iono, lambda, frequencies, obs_comb, nSatTot, w_bar, flag_XM_prep, sbas, cc, flag_full_prepro, order);
                        end
                    else
                        flag_XM_prep = 1;
                        [time_M_diff, pr1_M, ph1_M, pr2_M, ph2_M, ~, dtM, dtMdot, el_m, az_m, bad_sats_M, bad_epochs_M, var_dtM, var_SPP_M, status_obs_M, status_cs, eclipsed, ISBs, var_ISBs] = pre_processing(time_GPS_diff, time_M_diff, pos_M, pr1_M, ph1_M, pr2_M, ph2_M, dop1_M, dop2_M, snr1_M, Eph, SP3, iono, lambda, frequencies, 'NONE', nSatTot, w_bar, flag_XM_prep, sbas, cc, flag_full_prepro, order);
                     end
                     % the time ref is now changed
                     time_m = GPS_Time(zero_time/86400 + GPS_Time.GPS_ZERO, time_M_diff);
                     date_M = datevec(time_m.getMatlabTime);
                    
                    for f = 1 : size(time_R,3)
                        %pre-processing
                        log.addMessage(['Pre-processing rover observations (file ' filename_obs{f} ')...']);
                        w_bar.setBarLen(length(time_GPS));
                        w_bar.createNewBar('Pre-processing rover...');
                        
                        aprXR = pos_R;
                        if all(any(pos_R))
                            flag_XR = 1;
                        else
                            flag_XR = 0;
                        end
                        
                        %apply P1C1 DCBs if needed
                        if (flag_SP3 && ~isempty(SP3.DCB) && any(codeC1_R(:)))
                            avail_sat = any(lambda,2);
                            pr1_R(avail_sat,:,f) = pr1_R(avail_sat,:,f) + SP3.DCB.P1C1.value(avail_sat,ones(size(pr1_R(:,:,f),2),1))*1e-9*goGNSS.V_LIGHT.*codeC1_R(avail_sat,:,f);
                        end
                        
                        if state.isModeBlock()
                            [time_R_diff(:,1,f), pr1_R(:,:,f), ph1_R(:,:,f), pr2_R(:,:,f), ph2_R(:,:,f), pos_R_new(:,:,f), dtR(:,1,f), dtRdot(:,1,f), el_r(:,:,f), az_r(:,:,f), bad_sats_R(:,1,f), bad_epochs_R(:,1,f), var_dtR(:,1,f), var_SPP_R(:,:,f), status_obs_R(:,:,f), status_cs] = pp.execute(time_GPS_diff, time_R_diff(:,1,f), aprXR(:,:,f), pr1_R(:,:,f), ph1_R(:,:,f), pr2_R(:,:,f), ph2_R(:,:,f), dop1_R(:,:,f), dop2_R(:,:,f), snr1_R(:,:,f), iono, lambda, frequencies, obs_comb, nSatTot, w_bar, flag_XR, sbas, flag_full_prepro, order);
                        else
                            [time_R_diff(:,1,f), pr1_R(:,:,f), ph1_R(:,:,f), pr2_R(:,:,f), ph2_R(:,:,f), pos_R_new(:,:,f), dtR(:,1,f), dtRdot(:,1,f), el_r(:,:,f), az_r(:,:,f), bad_sats_R(:,1,f), bad_epochs_R(:,1,f), var_dtR(:,1,f), var_SPP_R(:,:,f), status_obs_R(:,:,f), status_cs] = pre_processing(time_GPS_diff, time_R_diff(:,1,f), aprXR(:,:,f), pr1_R(:,:,f), ph1_R(:,:,f), pr2_R(:,:,f), ph2_R(:,:,f), dop1_R(:,:,f), dop2_R(:,:,f), snr1_R(:,:,f), Eph, SP3, iono, lambda, frequencies, obs_comb, nSatTot, w_bar, flag_XR, sbas, cc, flag_full_prepro, order);
                        end
                        % the time is now changed
                        time_r(f) = GPS_Time(zero_time/86400 + GPS_Time.GPS_ZERO, time_R_diff(:,:,f)); %#ok<SAGROW>
                        date_R(:,:,f) = datevec(time_r(f).getMatlabTime);
                        
                        if report.opt.write == 1
                            report.prep.spp_threshold = SPP_threshold;
                            report.prep.flag_R = flag_XR;
                            report.prep.tot_epoch_R(f)=size(pr1_R(:,:,f),2);
                            report.prep.proc_epoch_R(f)=length(bad_epochs_R(isfinite(bad_epochs_R(:,1,f)),1,f));
                            report.prep.bad_epoch_R(f)=sum(bad_epochs_R(isfinite(bad_epochs_R(:,1,f)),1,f)==1);
                            if (~isempty(var_SPP_R(isfinite(var_SPP_R(:,1,f)),1,f)))
                                report.prep.max_varSPP_R(f)=max(var_SPP_R(isfinite(var_SPP_R(:,1,f)),1,f))^0.5;
                            else
                                report.prep.max_varSPP_R(f) = NaN;
                            end
                            report.prep.varSPP_R(f)=(sum(var_SPP_R(isfinite(var_SPP_R(:,2,f)),2,f))/sum(var_SPP_R(isfinite(var_SPP_R(:,2,f)),3,f)))^.5;
                            report.prep.tot_obs_R(f)=length(find(isfinite(status_obs_R(:,:,f))));
                            report.prep.obs_outlier_R(f)=length(find(status_obs_R(:,:,f)==-1));
                            report.prep.obs_used_R(f)=length(find(status_obs_R(:,:,f)==1));
                            report.prep.obs_undercutoff_R(f)=length(find(status_obs_R(:,:,f)==0));
                            report.prep.obs_stat_R(:,:,f)=[sum(status_obs_R(:,:,f)==0,2), sum(status_obs_R(:,:,f)==1,2), sum(status_obs_R(:,:,f)==-1,2)]; % [#under_cutoff, #used, #outlier] grouped by satellite
                            report.prep.CS_R{f}=status_cs;
                        end
                        
                        %-------------------------------------------------------------------------------------------
                        %% CYCLE SLIP / OUTLIERS DETECTION
                        %-------------------------------------------------------------------------------------------
                        
                        [ph1_R(:,:,f), ph1_M] = cycleSlipDetectSingleDiff(ph1_R(:,:,f), ph1_M);
                        if numel(frequencies) > 1
                            if ~state.isModeSEID()
                                % In seid mode ph2_M should be empty
                                [ph2_R(:,:,f), ph2_M] = cycleSlipDetectSingleDiff(ph2_R(:,:,f), ph2_M);
                            end
                        end
                        %[ph1_R(:,:,f), ph1_M] = cycle_slip_detect_single_diff(ph1_R(:,:,f), ph1_M, interval);
                        %[ph2_R(:,:,f), ph2_M] = cycle_slip_detect_single_diff(ph2_R(:,:,f), ph2_M, interval);
                        
                        w_bar.close();
                    end

                    if report.opt.write == 1
                        report.prep.tot_epoch_M=size(pr1_M,2);
                        report.prep.proc_epoch_M=length(bad_epochs_M(isfinite(bad_epochs_M)));
                        report.prep.bad_epoch_M=sum(bad_epochs_M(isfinite(bad_epochs_M))==1);
                        if (~isempty(var_SPP_M(isfinite(var_SPP_M(:,1)),1)))
                            report.prep.max_varSPP_M=max(var_SPP_M(isfinite(var_SPP_M(:,1)),1))^0.5;
                        else
                            report.prep.max_varSPP_M = NaN;
                        end
                        report.prep.varSPP_M=(sum(var_SPP_M(isfinite(var_SPP_M(:,2)),2))/sum(var_SPP_M(isfinite(var_SPP_M(:,2)),3)))^.5;
                        report.prep.tot_obs_M=length(find(isfinite(status_obs_M)));
                        report.prep.obs_outlier_M=length(find(status_obs_M==-1));
                        report.prep.obs_used_M=length(find(status_obs_M==1));
                        report.prep.obs_undercutoff_M=length(find(status_obs_M==0));
                        report.prep.obs_stat_M=[sum(status_obs_M==0,2), sum(status_obs_M==1,2), sum(status_obs_M==-1,2)];
                        report.prep.CS_M=status_cs;
                    end
                    
                    w_bar.close();
                end
                
                %         %read surveying mode
                %         if (flag_stopGOstop == 0)
                %             fid_dyn = fopen([filerootIN '_dyn_000.bin'],'r+');
                %             order = double(fread(fid_dyn,length(time_GPS),'uint8'));
                %             fclose(fid_dyn);
                %         end
                
                %TEMP
                snr_R = snr1_R;
                if (goGNSS.isDD(mode))
                    snr_M = snr1_M;
                end
                
                if (~flag_SP3)
                    %exclude satellites without ephemerides
                    delsat = setdiff(1:nSatTot,unique(Eph(30,:)));
                    %delsat = [delsat 4]; % exclude satellite 4
                    pr1_R(delsat,:,:) = 0;
                    pr2_R(delsat,:,:) = 0;
                    ph1_R(delsat,:,:) = 0;
                    ph2_R(delsat,:,:) = 0;
                    dop1_R(delsat,:,:) = 0;
                    dop2_R(delsat,:,:) = 0;
                    snr_R(delsat,:,:) = 0;
                    if (goGNSS.isDD(mode))
                        pr1_M(delsat,:,:) = 0;
                        pr2_M(delsat,:,:) = 0;
                        ph1_M(delsat,:,:) = 0;
                        ph2_M(delsat,:,:) = 0;
                        dop1_M(delsat,:,:) = 0;
                        dop2_M(delsat,:,:) = 0;
                        snr_M(delsat,:,:) = 0; %#ok<SAGROW>
                    end
                end
                
                %exclude flagged satellites (rover)
                if (exist('bad_sats_R','var'))
                    for f = 1 : size(pr1_R,3)
                        if (any(bad_sats_R(:,1,f)))
                            pos = find(bad_sats_R(:,1,f));
                            pr1_R(pos,:,f) = 0;
                            pr2_R(pos,:,f) = 0;
                            ph1_R(pos,:,f) = 0;
                            ph2_R(pos,:,f) = 0;
                            dop1_R(pos,:,f) = 0;
                            dop2_R(pos,:,f) = 0;
                            snr_R(pos,:,f) = 0;
                        end
                    end
                end
                
                %exclude flagged satellites (master)
                if (goGNSS.isDD(mode) && exist('bad_sats_M','var'))
                    if (any(bad_sats_M(:,1)))
                        pos = find(bad_sats_M(:,1));
                        pr1_M(pos,:) = 0;
                        pr2_M(pos,:) = 0;
                        ph1_M(pos,:) = 0;
                        ph2_M(pos,:) = 0;
                        dop1_M(pos,:) = 0;
                        dop2_M(pos,:) = 0;
                        snr_M(pos,:) = 0; %#ok<SAGROW>
                    end
                end
                
                %exclude flagged epochs (rover)
                if (exist('bad_epochs_R','var'))
                    for f = 1 : size(pr1_R,3)
                        if (any(bad_epochs_R(:,1,f) == 1))
                            pos = find(bad_epochs_R(:,1,f) == 1);
                            pr1_R(:,pos,f) = 0;
                            pr2_R(:,pos,f) = 0;
                            ph1_R(:,pos,f) = 0;
                            ph2_R(:,pos,f) = 0;
                            dop1_R(:,pos,f) = 0;
                            dop2_R(:,pos,f) = 0;
                            snr_R(:,pos,f) = 0;
                        end
                    end
                end
                
                %exclude flagged epochs (master)
                if (goGNSS.isDD(mode) && exist('bad_epochs_M','var'))
                    if (any(bad_epochs_M(:,1) == 1))
                        pos = find(bad_epochs_M(:,1) == 1);
                        pr1_M(:,pos) = 0;
                        pr2_M(:,pos) = 0;
                        ph1_M(:,pos) = 0;
                        ph2_M(:,pos) = 0;
                        dop1_M(:,pos) = 0;
                        dop2_M(:,pos) = 0;
                        snr_M(:,pos) = 0; %#ok<SAGROW>
                    end
                end
                
                %exclude eclipsed satellites (shadow crossing + 30 minutes; noon and midnight maneuvers)
                if (exist('eclipsed','var') && any(eclipsed(:)))
                    eclipse_map = diff(eclipsed,1,2);
                    [eclipsed_sat, eclipse_end]  = find(eclipse_map == -1);
                    %[midnight_sat, midnight_end] = find(eclipse_map == -2);
                    %[noon_sat,     noon_end]     = find(eclipse_map == -3);
                    %purge false shadow/maneuvers endings (due to missing observations)
                    for ee = length(eclipse_end) : -1 : 1
                        if (any(eclipse_map(eclipsed_sat(ee),eclipse_end(ee):end) == 1))
                            eclipsed_sat(ee) = [];
                            eclipse_end(ee)  = [];
                        end
                    end
                    extra_minutes = 30;
                    extra_epochs = extra_minutes*60/interval;
                    for e = 1 : length(eclipsed_sat)
                        idx1 = eclipse_end(e)+1;
                        idx2 = eclipse_end(e)+extra_epochs;
                        idx2 = min(idx2,size(eclipsed,2));
                        eclipsed(eclipsed_sat(e),idx1:idx2) = 1;
                    end
                    eclipsed(eclipsed>0) = 1;
                    for f = 1 : size(pr1_R,3)
                        pr1_R(:,:,f) = pr1_R(:,:,f).*~eclipsed;
                        pr2_R(:,:,f) = pr2_R(:,:,f).*~eclipsed;
                        ph1_R(:,:,f) = ph1_R(:,:,f).*~eclipsed;
                        ph2_R(:,:,f) = ph2_R(:,:,f).*~eclipsed;
                        dop1_R(:,:,f) = dop1_R(:,:,f).*~eclipsed;
                        dop2_R(:,:,f) = dop2_R(:,:,f).*~eclipsed;
                        snr_R(:,:,f) = snr_R(:,:,f).*~eclipsed;
                    end
                    if (goGNSS.isDD(mode))
                        pr1_M = pr1_M.*~eclipsed;
                        pr2_M = pr2_M.*~eclipsed;
                        ph1_M = ph1_M.*~eclipsed;
                        ph2_M = ph2_M.*~eclipsed;
                        dop1_M = dop1_M.*~eclipsed;
                        dop2_M = dop2_M.*~eclipsed;
                        snr_M = snr_M.*~eclipsed;
                    end
                end
                
                %exclude CRX-flagged satellites
                if (exist('CRX','var') && any(CRX(:)))
                    for f = 1 : size(pr1_R,3)
                        pr1_R(:,:,f) = pr1_R(:,:,f).*~CRX;
                        pr2_R(:,:,f) = pr2_R(:,:,f).*~CRX;
                        ph1_R(:,:,f) = ph1_R(:,:,f).*~CRX;
                        ph2_R(:,:,f) = ph2_R(:,:,f).*~CRX;
                        dop1_R(:,:,f) = dop1_R(:,:,f).*~CRX;
                        dop2_R(:,:,f) = dop2_R(:,:,f).*~CRX;
                        snr_R(:,:,f) = snr_R(:,:,f).*~CRX;
                    end
                    if (goGNSS.isDD(mode))
                        pr1_M = pr1_M.*~CRX;
                        pr2_M = pr2_M.*~CRX;
                        ph1_M = ph1_M.*~CRX;
                        ph2_M = ph2_M.*~CRX;
                        dop1_M = dop1_M.*~CRX;
                        dop2_M = dop2_M.*~CRX;
                        snr_M = snr_M.*~CRX;
                    end
                end
                
                %%reverse the path
                %pr1_R = pr1_R(:,end:-1:1);
                %pr1_M = pr1_M(:,end:-1:1);
                %ph1_R = ph1_R(:,end:-1:1);
                %ph1_M = ph1_M(:,end:-1:1);
                %pr2_R = pr2_R(:,end:-1:1);
                %pr2_M = pr2_M(:,end:-1:1);
                %ph2_R = ph2_R(:,end:-1:1);
                %ph2_M = ph2_M(:,end:-1:1);
                %dop1_R = dop1_R(:,end:-1:1);
                %dop1_M = dop1_M(:,end:-1:1);
                %dop2_R = dop2_R(:,end:-1:1);
                %dop2_M = dop2_M(:,end:-1:1);
                %snr_R = snr_R(:,end:-1:1);
                %snr_M = snr_M(:,end:-1:1);
                
                %time_GPS = time_GPS(end:-1:1);
                %date_R = date_R(end:-1:1,:);
                
                %if relative post-processing positioning (i.e. with master station)
                if goGNSS.isDD(mode) && goGNSS.isPP(mode)
                    %master station position management
                    %         if (flag_ms_pos) && (sum(abs(pos_M)) ~= 0)
                    %             if (size(pos_M,2) == 1)
                    %                 pos_M(1,1:length(time_GPS)) = pos_M(1);
                    %                 pos_M(2,1:length(time_GPS)) = pos_M(2);
                    %                 pos_M(3,1:length(time_GPS)) = pos_M(3);
                    %             end
                    %         else
                    pos_M(1,1:length(time_GPS)) = pos_M(1);
                    pos_M(2,1:length(time_GPS)) = pos_M(2);
                    pos_M(3,1:length(time_GPS)) = pos_M(3);
                    %             if (exist('is_batch','var') && any(pos_M_off))
                    %                 pos_M_disp = pos_M_off;
                    %             else
                    %                 pos_M_disp(1,1) = pos_M_man(1,1);
                    %                 pos_M_disp(2,1) = pos_M_man(2,1);
                    %                 pos_M_disp(3,1) = pos_M_man(3,1);
                    %             end
                    %             fprintf('Master position fixed to user-defined values:\n');
                    %             fprintf(' X=%.4f m, Y=%.4f m, Z=%.4f m\n', pos_M_disp(1,1), pos_M_disp(2,1), pos_M_disp(3,1));
                    %         end
                    
                    %if master station data are not available
                    if (~any(pr1_M(:)))
                        
                        %switch from relative to absolute positioning...
                        mode = mode - 10;
                        
                        %...and warn the user
                        if (mode_user == 1)
                            if (flag_var_dyn_model)
                                uiwait(msgbox('Warning: master data not available, forcing undifferenced mode. Variable dynamic model is not supported in undifferenced mode.','','modal'));
                            else
                                uiwait(msgbox('Warning: master data not available, forcing undifferenced mode.','','modal'));
                            end
                        else
                            if (flag_var_dyn_model)
                                log.addWarning('Master data not available, forcing undifferenced mode. Variable dynamic model is not supported in undifferenced mode.');
                            else
                                log.addWarning('Master data not available, forcing undifferenced mode.');
                            end
                        end
                    end
                end
                
                
                %global residuals_fixed residuals_float outliers s02_ls %#ok<TLEV>
                %residuals_fixed=NaN(2*length(n_freq)*nSatTot,1);
                %residuals_float=NaN(2*length(n_freq)*nSatTot,1);
                %outliers=zeros(2*length(n_freq)*nSatTot,1);
                %global min_ambfixRMS min_ambfloatRMS
                %min_ambfixRMS=NaN(length(time_GPS),1);
                %s02_ls=NaN(length(time_GPS),1);
                
            else %real-time
                
                %disable Doppler-based cycle-slips detection
                flag_doppler_cs = 0;
                
                %initialize master position variable
                if (flag_ms_pos)
                    pos_M = [];
                else
                    pos_M = pos_M_man;
                end
                
                %for the Kalman filter execution in real-time
                dop1_M = zeros(nSatTot,1);
                pr2_M  = zeros(nSatTot,1);
                pr2_R  = zeros(nSatTot,1);
                ph2_M  = zeros(nSatTot,1);
                ph2_R  = zeros(nSatTot,1);
                dop2_M = zeros(nSatTot,1);
                dop2_R = zeros(nSatTot,1);
            end
            
            %check if the dataset was surveyed with a variable dynamic model
            d = dir([state.obs_dir '\' state.obs_name{1}(1:end-4) '_dyn.bin']);
            if (goGNSS.isPP(mode) && (flag_stopGOstop || flag_var_dyn_model) && isempty(d))
                log.addWarning(' Dataset was not surveyed with a variable dynamic model:');
                log.addMessage('      Switching off variable dynamic model mode...');
                flag_var_dyn_model = 0;
            end
            
            %boolean vector for removing unused epochs in LS processing
            if (goGNSS.isPP(mode))
                unused_epochs = zeros(size(time_GPS));
            end
            
            %update variance of tropospheric delay
            global sigmaq_tropo sigmaq_tropo_gradient %#ok<TLEV>
            % these values are set in Go_State.settingsToGo
            %sigmaq_tropo = (0.001/sqrt(3600/interval))^2;
            %sigmaq_tropo_gradient = (0.0002/sqrt(3600/interval))^2;
            
            %----------------------------------------------------------------------------------------------
            %% END OF SEID (Satellite-specific Epoch-differenced Ionospheric Delay)
            %----------------------------------------------------------------------------------------------
            
            if (state.isModeSEID())
                SEID_main;
                
                if (mode == goGNSS.MODE_PP_SEID_PPP)
                    mode = goGNSS.MODE_PP_KF_CP_SA; % Switching from SEID PPP to PPP
                    state.setMode(mode); % Switching from SEID PPP to PPP
                    read_files = true; % In case of SEID processing the read operation must be repeated twice
                    flag_SEID = false;
                    clear pos_R_crd;
                end
            end
        end
        
        if (~state.isModeSEID() || mode == goGNSS.MODE_PP_SEID_PPP)
            %----------------------------------------------------------------------------------------------
            %% INITIALIZATION - SAVING SETTINGS USED FOR THE COMPUTATION
            %----------------------------------------------------------------------------------------------
            
            if state.isOutSettings()
                state.save([filerootOUT '_settings.ini']);
            end
            
            switch (mode)
                case goGNSS.MODE_PP_LS_C_SA
                    
                    %----------------------------------------------------------------------------------------------
                    %% POST-PROCESSING (ABSOLUTE POSITIONING): LEAST SQUARES ON CODE
                    %----------------------------------------------------------------------------------------------
                    
                    fid_kal = fopen([filerootOUT '_kal_000.bin'],'w+');
                    fid_sat = fopen([filerootOUT '_sat_000.bin'],'w+');
                    fid_dop = fopen([filerootOUT '_dop_000.bin'],'w+');
                    fid_conf = fopen([filerootOUT '_conf_000.bin'],'w+');
                    fid_res = fopen([filerootOUT '_res_000.bin'],'w+');
                    residuals_dummy = NaN(1,nSatTot);
                    
                    [pr1_R, ph1_R, pr2_R, ph2_R] = multi_GNSS_biases_correction(time_GPS_diff, pr1_R, ph1_R, pr2_R, ph2_R, ISBs, Eph, constellations, lambda);
                    
                    nN = nSatTot;
                    nT = 0;
                    check_on = 0;
                    check_off = 0;
                    check_pivot = 0;
                    check_cs = 0;
                    
                    plot_t = 1;
                    
                    % goGPS waiting bar
                    w_bar.setBarLen(length(time_GPS));
                    w_bar.createNewBar('Processing rover...');
                    
                    for t = 1 : length(time_GPS)
                        
                        Eph_t = rt_find_eph(Eph, time_GPS_diff(t), nSatTot);
                        
                        sbas_t = find_sbas(sbas, t);
                        
                        goGPS_LS_SA_code(time_GPS_diff(t), pr1_R(:,t), pr2_R(:,t), snr_R(:,t), Eph_t, SP3, iono, sbas_t, lambda, frequencies, state.getProcessingRate(), obs_comb, pos_R);
                        
                        if (t == 1)
                            fwrite(fid_sat, nSatTot, 'int8');
                            fwrite(fid_conf, nSatTot, 'int8');
                            fwrite(fid_res, nSatTot, 'int8');
                        end
                        
                        if ~isempty(Xhat_t_t) && ~any(isnan([Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)]))
                            Xhat_t_t_dummy = [Xhat_t_t; zeros(nN,1)];
                            Cee_dummy = [Cee zeros(o3,nN); zeros(nN,o3) zeros(nN,nN)];
                            fwrite(fid_kal, [Xhat_t_t_dummy; Cee_dummy(:)], 'double');
                            fwrite(fid_sat, [zeros(nSatTot,1); azR; zeros(nSatTot,1); elR; zeros(nSatTot,1); distR], 'double');
                            fwrite(fid_dop, [PDOP; HDOP; VDOP; 0; 0; 0], 'double');
                            fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
                            fwrite(fid_res, [residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_float(1:nSatTot); residuals_dummy'; residuals_dummy'; residuals_dummy';outliers(1:nSatTot); residuals_dummy'; residuals_dummy'; residuals_dummy'], 'double');
                            
                            if (flag_plotproc)
                                if (mode_user == 1 && t == 1), w_bar.shiftDown(); end
                                if (flag_cov == 0)
                                    if (flag_ge == 1), rtplot_googleearth (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), date_R(t,:)), end
                                    rtplot_matlab (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                                else
                                    if (flag_ge == 1), rtplot_googleearth_cov (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(t,:)), end
                                    rtplot_matlab_cov (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                                end
                                if (flag_skyplot == 1)
                                    rtplot_skyplot (plot_t, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                    rtplot_snr (snr_R(:,t), Eph_t, SP3);
                                else
                                    rttext_sat (plot_t, azR, elR, snr_R(:,t), conf_sat, pivot, Eph_t, SP3);
                                end
                                plot_t = plot_t + 1;
                                pause(0.01);
                            end
                        else
                            unused_epochs(t) = 1;
                        end
                        
                        w_bar.goTime(t);
                    end
                    w_bar.close();
                    
                    fclose(fid_kal);
                    fclose(fid_sat);
                    fclose(fid_dop);
                    fclose(fid_conf);
                    fclose(fid_res);
                    
                case goGNSS.MODE_PP_KF_C_SA
                    
                    %----------------------------------------------------------------------------------------------
                    %% POST-PROCESSING (ABSOLUTE POSITIONING): KALMAN FILTER ON CODE
                    %----------------------------------------------------------------------------------------------
                    
                    fid_kal = fopen([filerootOUT '_kal_000.bin'],'w+');
                    fid_sat = fopen([filerootOUT '_sat_000.bin'],'w+');
                    fid_dop = fopen([filerootOUT '_dop_000.bin'],'w+');
                    fid_conf = fopen([filerootOUT '_conf_000.bin'],'w+');
                    fid_res = fopen([filerootOUT '_res_000.bin'],'w+');
                    
                    nN = nSatTot;
                    nT = 0;
                    check_on = 0;
                    check_off = 0;
                    check_pivot = 0;
                    check_cs = 0;
                    
                    if ~state.isSeamlessKF()
                        kalman_initialized = 0;
                    end
                    if (~kalman_initialized)
                        while (~kalman_initialized)
                            if (isempty(time_GPS))
                                fprintf('It was not possible to initialize the Kalman filter.\n');
                                return
                            end
                            
                            Eph_t = rt_find_eph (Eph, time_GPS_diff(1), nSatTot);
                            
                            sbas_t = find_sbas(sbas, 1);
                            
                            kalman_initialized = goGPS_KF_SA_code_init(pos_R, time_GPS_diff(1), pr1_R(:,1), pr2_R(:,1), snr_R(:,1), Eph_t, SP3, iono, sbas_t, lambda, frequencies(1), state.getProcessingRate());
                            
                            if (~kalman_initialized)
                                time_GPS_diff(1) = []; time_GPS(1) = []; week_R(1) = [];
                                if not(isempty(intersect(frequencies,1)))
                                    pr1_R(:,1) = []; ph1_R(:,1) = []; dop1_R(:,1) = [];
                                end
                                if not(isempty(intersect(frequencies,2)))
                                    pr2_R(:,1) = []; ph2_R(:,1) = []; dop2_R(:,1) = [];
                                end
                                snr_R(:,1) = [];
                            end
                        end
                        
                        Xhat_t_t_dummy = [Xhat_t_t; zeros(nN,1)];
                        Cee_dummy = [Cee zeros(o3,nN); zeros(nN,o3) zeros(nN,nN)];
                        fwrite(fid_kal, [Xhat_t_t_dummy; Cee_dummy(:)], 'double');
                        fwrite(fid_sat, nSatTot, 'int8');
                        fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                        fwrite(fid_dop, [PDOP; HDOP; VDOP; KPDOP; KHDOP; KVDOP], 'double');
                        fwrite(fid_conf, [nSatTot; conf_sat; conf_cs; pivot], 'int8');
                        fwrite(fid_res, nSatTot, 'int8');
                        fwrite(fid_res, [residuals_fixed(1:nSatTot*2); residuals_fixed(nSatTot*2+1:end);residuals_float(1:nSatTot*2); residuals_float(nSatTot*2+1:end);outliers(1:nSatTot*2);outliers(nSatTot*2+1:end)], 'double');
                        
                        if (flag_plotproc)
                            if (flag_cov == 0)
                                if (flag_ge == 1), rtplot_googleearth (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), date_R(1,:)), end
                                rtplot_matlab (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), 0, 0, 0, 0, flag_ms, gs.getReferencePath());
                            else
                                if (flag_ge == 1), rtplot_googleearth_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(1,:)), end
                                rtplot_matlab_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), 0, 0, 0, 0, flag_ms, gs.getReferencePath());
                            end
                            if (flag_skyplot == 1)
                                rtplot_skyplot (1, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                rtplot_snr (snr_R(:,1), Eph_t, SP3);
                            else
                                rttext_sat (1, azR, elR, snr_R(:,1), conf_sat, pivot, Eph_t, SP3);
                            end
                        end
                        t1 = 2;
                    else
                        fwrite(fid_sat, nSatTot, 'int8');
                        fwrite(fid_conf, nSatTot, 'int8');
                        fwrite(fid_res, nSatTot, 'int8');
                        t1 = 1;
                    end
                    
                    w_bar.setBarLen(length(time_GPS));
                    w_bar.createNewBar('Processing...');
                    
                    for t = t1 : length(time_GPS)
                        residuals_fixed=NaN(4*nSatTot,1);
                        residuals_float=NaN(4*nSatTot,1);
                        outliers=zeros(4*nSatTot,1);
                        
                        Eph_t = rt_find_eph (Eph, time_GPS_diff(t), nSatTot);
                        
                        goGPS_KF_SA_code_loop(time_GPS_diff(t), pr1_R(:,t), pr2_R(:,t), snr_R(:,t), Eph_t, SP3, iono, sbas_t, lambda, frequencies(1), state.getProcessingRate());
                        
                        Xhat_t_t_dummy = [Xhat_t_t; zeros(nN,1)];
                        Cee_dummy = [Cee zeros(o3,nN); zeros(nN,o3) zeros(nN,nN)];
                        fwrite(fid_kal, [Xhat_t_t_dummy; Cee_dummy(:)], 'double');
                        fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                        fwrite(fid_dop, [PDOP; HDOP; VDOP; KPDOP; KHDOP; KVDOP], 'double');
                        fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
                        fwrite(fid_res, [residuals_fixed(1:nSatTot*2); residuals_fixed(nSatTot*2+1:end);residuals_float(1:nSatTot*2); residuals_float(nSatTot*2+1:end);outliers(1:nSatTot*2);outliers(nSatTot*2+1:end)], 'double');
                        
                        if (flag_plotproc)
                            if (mode_user == 1 && t == 2), w_bar.shiftDown(); end
                            if (flag_cov == 0)
                                if (flag_ge == 1), rtplot_googleearth (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), date_R(t,:)), end
                                rtplot_matlab (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                            else
                                if (flag_ge == 1), rtplot_googleearth_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(t,:)), end
                                rtplot_matlab_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                            end
                            if (flag_skyplot == 1)
                                rtplot_skyplot (t, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                rtplot_snr (snr_R(:,t), Eph_t, SP3);
                            else
                                rttext_sat (t, azR, elR, snr_R(:,t), conf_sat, pivot, Eph_t, SP3);
                            end
                            pause(0.01);
                        end
                        
                        w_bar.goTime(t);
                    end
                    
                    w_bar.close();
                    
                    fclose(fid_kal);
                    fclose(fid_sat);
                    fclose(fid_dop);
                    fclose(fid_conf);
                    fclose(fid_res);
                    
                case goGNSS.MODE_PP_LS_CP_SA
                    
                    %----------------------------------------------------------------------------------------------
                    %% POST-PROCESSING (ABSOLUTE POSITIONING): LEAST SQUARES ON CODE AND PHASE   (DISABLED IN GUI)
                    %----------------------------------------------------------------------------------------------
                    
                    fid_kal = fopen([filerootOUT '_kal_000.bin'],'w+');
                    fid_sat = fopen([filerootOUT '_sat_000.bin'],'w+');
                    fid_dop = fopen([filerootOUT '_dop_000.bin'],'w+');
                    fid_conf = fopen([filerootOUT '_conf_000.bin'],'w+');
                    fid_res = fopen([filerootOUT '_res_000.bin'],'w+');
                    
                    nN = nSatTot;
                    nT = 0;
                    check_on = 0;
                    check_off = 0;
                    check_pivot = 0;
                    check_cs = 0;
                    
                    plot_t = 1;
                    
                    
                    % goGPS waiting bar
                    w_bar.setBarLen(length(time_GPS));
                    w_bar.createNewBar('Processing...');
                    
                    for t = 1 : length(time_GPS)
                        
                        Eph_t = rt_find_eph (Eph, time_GPS_diff(t), nSatTot);
                        
                        sbas_t = find_sbas(sbas, t);
                        
                        goGPS_LS_SA_code_phase(time_GPS_diff(t), pr1_R(:,t), pr2_R(:,t), ph1_R(:,t), ph2_R(:,t), snr_R(:,t), Eph_t, SP3, iono, sbas_t, lambda, frequencies(1), state.getProcessingRate());
                        
                        if (t == 1)
                            fwrite(fid_sat, nSatTot, 'int8');
                            fwrite(fid_conf, nSatTot, 'int8');
                            fwrite(fid_res, nSatTot, 'int8');
                        end
                        
                        if ~isempty(Xhat_t_t) && ~any(isnan([Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)]))
                            fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
                            fwrite(fid_sat, [zeros(nSatTot,1); azR; zeros(nSatTot,1); elR; zeros(nSatTot,1); distR], 'double');
                            fwrite(fid_dop, [PDOP; HDOP; VDOP; 0; 0; 0], 'double');
                            fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
                            residuals_dummy = NaN(1,nSatTot);
                            fwrite(fid_res, [residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'], 'double');
                            
                            if (flag_plotproc)
                                if (mode_user == 1 && t == 1), w_bar.shiftDown(); end
                                if (flag_cov == 0)
                                    if (flag_ge == 1), rtplot_googleearth (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), date_R(t,:)), end
                                    rtplot_matlab (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                                else
                                    if (flag_ge == 1), rtplot_googleearth_cov (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(t,:)), end
                                    rtplot_matlab_cov (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                                end
                                if (flag_skyplot == 1)
                                    rtplot_skyplot (plot_t, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                    rtplot_snr (snr_R(:,t), Eph_t, SP3);
                                else
                                    rttext_sat (plot_t, azR, elR, snr_R(:,t), conf_sat, pivot, Eph_t, SP3);
                                end
                                plot_t = plot_t + 1;
                                pause(0.01);
                            end
                        else
                            unused_epochs(t) = 1;
                        end
                        
                        w_bar.goTime(t);
                    end
                    
                    w_bar.close();
                    
                    fclose(fid_kal);
                    fclose(fid_sat);
                    fclose(fid_dop);
                    fclose(fid_conf);
                    fclose(fid_res);
                    
                case goGNSS.MODE_PP_LS_CP_VEL
                    
                    %----------------------------------------------------------------------------------------------
                    %% VARIOMETRIC APPROACH FOR VELOCITY ESTIMATION STAND-ALONE ENGINE
                    %----------------------------------------------------------------------------------------------
                    
                    fid_kal = fopen([filerootOUT '_kal_000.txt'],'w');
                    fid_sat = fopen([filerootOUT '_sat_000.bin'],'w+');
                    fid_dop = fopen([filerootOUT '_dop_000.bin'],'w+');
                    fid_conf = fopen([filerootOUT '_conf_000.bin'],'w+');
                    fid_res = fopen([filerootOUT '_res_000.bin'],'w+');
                    
                    nN = nSatTot;
                    nT = 0;
                    check_on = 0;
                    check_off = 0;
                    check_pivot = 0;
                    check_cs = 0;
                    
                    plot_t = 1;
                    time_step = state.getVariometricTimeStep();   % time step to perform the difference between phase observations
                    fprintf('TimeStep used is %d epochs\n', time_step);
                    % External loop to show bar update every 15 epochs
                    stepUpdate = 15;
                    w_bar.setBarLen((length(time_GPS)-(time_step))/stepUpdate);
                    w_bar.createNewBar('Variometric approach running...');
                    ind=0;
                    for tExt = 1:stepUpdate:(length(time_GPS)-(time_step))
                        for t = tExt:min(tExt+stepUpdate-1,length(time_GPS)-(time_step))
                            Eph_t = rt_find_eph (Eph, time_GPS_diff(t), nSatTot);
                            Eph_t1 = rt_find_eph (Eph, time_GPS_diff(t+time_step), nSatTot);
                            
                            sbas_t = find_sbas(sbas, t);
                            sbas_t1 = find_sbas(sbas, t+time_step);
                            
                            goGPS_LS_SA_variometric(time_GPS_diff(t), time_GPS_diff(t+time_step), pr1_R(:,t), pr1_R(:,t+time_step), pr2_R(:,t), pr2_R(:,t+time_step), ph1_R(:,t), ph1_R(:,t+time_step), ph2_R(:,t), ph2_R(:,t+time_step), snr_R(:,t), snr_R(:,t+time_step), Eph_t, Eph_t1, [], [], iono, sbas_t, sbas_t1, lambda, frequencies(1), state.getProcessingRate());
                            Xhat_t_t(1:6)=-Xhat_t_t(1:6)./(interval.*time_step);
                            if ~isempty(Xhat_t_t) && ~any(isnan([Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)]))
                                Xhat_t_t_dummy = [Xhat_t_t]; %#ok<NBRAK>
                                ind=ind+1;
                                vel_pos(ind,:) = Xhat_t_t; %#ok<SAGROW>
                                Cee_dummy = Cee;
                                %         fprintf(fid_kal,'%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n', Xhat_t_t );
                                if (t == 1)
                                    fwrite(fid_sat, nSatTot, 'int8');
                                    fwrite(fid_conf, nSatTot, 'int8');
                                    fwrite(fid_res, nSatTot, 'int8');
                                end
                                fwrite(fid_sat, [zeros(nSatTot,1); azR; zeros(nSatTot,1); elR; zeros(nSatTot,1); distR], 'double');
                                fwrite(fid_dop, [PDOP; HDOP; VDOP; 0; 0; 0], 'double');
                                fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
                                residuals_dummy = NaN(1,nSatTot);
                                fwrite(fid_res, [residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'], 'double');
                                
                                %if (flag_plotproc)
                                %    if (flag_cov == 0)
                                %        if (flag_ge == 1), rtplot_googleearth (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), date_R(t,:)), end;
                                %        rtplot_matlab (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                                %    else
                                %        if (flag_ge == 1), rtplot_googleearth_cov (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(t,:)), end;
                                %        rtplot_matlab_cov (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                                %    end
                                %    if (flag_skyplot == 1)
                                %        rtplot_skyplot (plot_t, azR, elR, conf_sat, pivot);
                                %        rtplot_snr (snr_R(:,t));
                                %    else
                                %        rttext_sat (plot_t, azR, elR, snr_R(:,t), conf_sat, pivot);
                                %    end
                                %    plot_t = plot_t + 1;
                                %    %pause(0.01);
                                %end
                            else
                                unused_epochs(t) = 1;
                            end
                        end
                        w_bar.goTime(tExt/stepUpdate);
                    end
                    goDX=cumsum(vel_pos(:,1)).*(interval.*time_step);
                    goDY=cumsum(vel_pos(:,3)).*(interval.*time_step);
                    goDZ=cumsum(vel_pos(:,5)).*(interval.*time_step);
                    %jumps=0.*goDX;
                    jumps=(find(vel_pos(:,2)==-9999));
                    jumps=[0; jumps; length(goDX)-1]; %#ok<AGROW>
                    mX = zeros(length(jumps)-1,1);
                    mY = zeros(length(jumps)-1,1);
                    mZ = zeros(length(jumps)-1,1);
                    lastX=0;
                    lastY=0;
                    lastZ=0;
                    for epo=1:length(jumps)-1
                        if (jumps(epo+1)-jumps(epo))~=1
                            mX(epo) = mean(vel_pos(jumps(epo)+1:jumps(epo+1)-1,7));
                            mY(epo) = mean(vel_pos(jumps(epo)+1:jumps(epo+1)-1,8));
                            mZ(epo) = mean(vel_pos(jumps(epo)+1:jumps(epo+1)-1,9));
                            lastX = goDX(jumps(epo)+1:jumps(epo+1)-1) - mean(goDX(jumps(epo)+1:jumps(epo+1)-1)) + mX(epo);
                            lastY = goDY(jumps(epo)+1:jumps(epo+1)-1) - mean(goDY(jumps(epo)+1:jumps(epo+1)-1)) + mY(epo);
                            lastZ = goDZ(jumps(epo)+1:jumps(epo+1)-1) - mean(goDZ(jumps(epo)+1:jumps(epo+1)-1)) + mZ(epo);
                            goDX(jumps(epo)+1:jumps(epo+1)-1)= lastX;
                            goDY(jumps(epo)+1:jumps(epo+1)-1)= lastY;
                            goDZ(jumps(epo)+1:jumps(epo+1)-1)= lastZ;
                        elseif (jumps(epo) ~= 0)
                            mX(epo) = mean(vel_pos(jumps(epo)+1:jumps(epo+1)-1,7));
                            mY(epo) = mean(vel_pos(jumps(epo)+1:jumps(epo+1)-1,8));
                            mZ(epo) = mean(vel_pos(jumps(epo)+1:jumps(epo+1)-1,9));
                            goDX(jumps(epo)) = lastX(end);
                            goDY(jumps(epo)) = lastY(end);
                            goDZ(jumps(epo)) = lastZ(end);
                        end
                    end
                    flag = find(diff(goDX) > 100);
                    goDX(flag)=goDX(flag-1);
                    goDY(flag)=goDY(flag-1);
                    goDZ(flag)=goDZ(flag-1);
                    %goDX=goDX-goDX(1)+(vel_pos(1,7));
                    %goDY=goDY-goDY(1)+(vel_pos(1,8));
                    %goDZ=goDZ-goDZ(1)+(vel_pos(1,9));
                    [phiX, lamX, hX] = cart2geod(goDX,goDY,goDZ);
                    vel_pos(:,7)=goDX;
                    vel_pos(:,8)=goDY;
                    vel_pos(:,9)=goDZ;
                    vel_pos(:,10)=phiX.*180/pi;
                    vel_pos(:,11)=lamX.*180/pi;
                    vel_pos(:,12)=hX;
                    
                    for epo=1:length(goDX)
                        
                        xENU(epo,:) = global2localPos(vel_pos(epo,7:9)', vel_pos(1,7:9)'); %#ok<SAGROW>
                        vENU(epo,:) = global2localVel(vel_pos(epo,1:2:5)', [phiX(epo), lamX(epo)]'.*180/pi); %#ok<SAGROW>
                        velpos(epo,:)=[vel_pos(epo,:), vENU(epo,:), time_GPS_diff(epo)]; %#ok<SAGROW>
                        
                        fprintf(fid_kal,'%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n',velpos(epo,:));
                    end
                    
                    % figure
                    % plot(xENU(:,1))
                    % hold on
                    % plot(xENU(:,2),'r')
                    % plot(xENU(:,3),'g')
                    % title('Displacement (blue=E; red=N; green=U)')
                    
                    if (mode_user == 1)
                        figure
                        plot(vENU(:,1))
                        hold on
                        plot(vENU(:,2),'r')
                        plot(vENU(:,3),'g')
                        title('Velocity (blue=E; red=N; green=U)')
                    end
                    
                    w_bar.close();
                    
                    fclose(fid_kal);
                    fclose(fid_sat);
                    fclose(fid_dop);
                    fclose(fid_conf);
                    fclose(fid_res);
                    
                    time_GPS_diff = time_GPS_diff(1:end-time_step);
                    time_GPS = time_GPS(1:end-time_step);
                    week_R = week_R(1:end-time_step);
                    
                case goGNSS.MODE_PP_KF_CP_SA
                    
                    %----------------------------------------------------------------------------------------------
                    %% POST-PROCESSING (ABSOLUTE POSITIONING): KALMAN FILTER ON CODE AND PHASE (PPP)
                    %----------------------------------------------------------------------------------------------
                    
                    cs_threshold = 1e30;
                    
                    fid_kal = fopen([filerootOUT '_kal_000.bin'],'w+');
                    fid_sat = fopen([filerootOUT '_sat_000.bin'],'w+');
                    fid_dop = fopen([filerootOUT '_dop_000.bin'],'w+');
                    fid_conf = fopen([filerootOUT '_conf_000.bin'],'w+');
                    fid_res = fopen([filerootOUT '_res_000.bin'],'w+');
                    fid_trp = fopen([filerootOUT '_trp_000.bin'],'w+');
                    
                    if ~state.isSeamlessKF()
                        kalman_initialized = 0;
                    end
                    if (~kalman_initialized)
                        
                        [pr1_R, ph1_R, pr2_R, ph2_R] = multi_GNSS_biases_correction(time_GPS_diff, pr1_R, ph1_R, pr2_R, ph2_R, ISBs, Eph, constellations, lambda);
                        ISBs_init = ISBs;
                        
                        while (~kalman_initialized)
                            if (isempty(time_GPS))
                                fprintf('It was not possible to initialize the Kalman filter.\n');
                                return
                            end
                            
                            Eph_t = rt_find_eph (Eph, time_GPS_diff(1), nSatTot);
                            
                            sbas_t = find_sbas(sbas, 1);
                            
                            kalman_initialized = goGPS_KF_SA_code_phase_init(pos_R, time_GPS_diff(1), pr1_R(:,1), ph1_R(:,1), dop1_R(:,1), pr2_R(:,1), ph2_R(:,1), dop2_R(:,1), snr_R(:,1), Eph_t, SP3, iono, sbas_t, lambda, frequencies, obs_comb, state.getProcessingRate(), flag_XR, flag_tropo, flag_tropo_gradient);
                            
                            if (~kalman_initialized)
                                time_GPS_diff(1) = []; time_GPS(1) = []; week_R(1) = [];
                                pr1_R(:,1) = []; ph1_R(:,1) = []; dop1_R(:,1) = [];
                                pr2_R(:,1) = []; ph2_R(:,1) = []; dop2_R(:,1) = [];
                                snr_R(:,1) = [];
                            end
                        end
                        
                        fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
                        fwrite(fid_sat, nSatTot, 'int8');
                        fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                        fwrite(fid_dop, [PDOP; HDOP; VDOP; KPDOP; KHDOP; KVDOP], 'double');
                        fwrite(fid_conf, [nSatTot; conf_sat; conf_cs; pivot], 'int8');
                        fwrite(fid_res, nSatTot, 'int8');
                        fwrite(fid_res, [residuals_fixed(1:nSatTot*2); residuals_fixed(nSatTot*2+1:end);residuals_float(1:nSatTot*2); residuals_float(nSatTot*2+1:end);outliers(1:nSatTot*2);outliers(nSatTot*2+1:end)], 'double');
                        fwrite(fid_trp, nSatTot, 'int8');
                        fwrite(fid_trp, [apriori_ZHD; STDs], 'double');
                        
                        if (flag_plotproc)
                            if (flag_cov == 0)
                                if (flag_ge == 1), rtplot_googleearth (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), date_R(1,:)), end
                                rtplot_matlab (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), 0, 0, 0, 0, flag_ms, gs.getReferencePath(), flag_amb);
                            else
                                if (flag_ge == 1), rtplot_googleearth_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(1,:)), end
                                rtplot_matlab_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), 0, 0, 0, 0, flag_ms, gs.getReferencePath(), flag_amb);
                            end
                            if (flag_amb == 1)
                                rtplot_amb (1, window, Xhat_t_t(o3+1:o3+nSatTot), sqrt(diag(Cee(o3+1:o3+nSatTot,o3+1:o3+nSatTot))), conf_cs)
                            else
                                if (flag_skyplot == 1)
                                    rtplot_skyplot (1, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                    rtplot_snr (snr_R(:,1), Eph_t, SP3);
                                else
                                    rttext_sat (1, azR, elR, snr_R(:,1), conf_sat, pivot, Eph_t, SP3);
                                end
                            end
                        end
                        t1 = 2;
                    else
                        fwrite(fid_sat, nSatTot, 'int8');
                        fwrite(fid_conf, nSatTot, 'int8');
                        fwrite(fid_res, nSatTot, 'int8');
                        fwrite(fid_trp, nSatTot, 'int8');
                        t1 = 1;
                        [pr1_R, ph1_R, pr2_R, ph2_R] = multi_GNSS_biases_correction(time_GPS_diff, pr1_R, ph1_R, pr2_R, ph2_R, ISBs_init, Eph, constellations, lambda);
                    end
                    
                    % goGPS waiting bar
                    w_bar.setBarLen(length(time_GPS));
                    w_bar.createNewBar('Processing...');
                    
                    for t = t1 : length(time_GPS)
                        residuals_fixed=NaN(4*nSatTot,1);
                        residuals_float=NaN(4*nSatTot,1);
                        outliers=zeros(4*nSatTot,1);
                        
                        Eph_t = rt_find_eph (Eph, time_GPS_diff(t), nSatTot);
                        
                        sbas_t = find_sbas(sbas, t);
                        
                        [check_on, check_off, check_pivot, check_cs] = goGPS_KF_SA_code_phase_loop(time_GPS_diff(t), pr1_R(:,t), ph1_R(:,t), dop1_R(:,t), pr2_R(:,t), ph2_R(:,t), dop2_R(:,t), snr_R(:,t), Eph_t, SP3, iono, sbas_t, lambda, frequencies, obs_comb, state.getProcessingRate(), flag_tropo, flag_tropo_gradient, antenna_PCV, antenna_PCV_S);
                        
                        fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
                        fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                        fwrite(fid_dop, [PDOP; HDOP; VDOP; KPDOP; KHDOP; KVDOP], 'double');
                        fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
                        fwrite(fid_res, [residuals_fixed(1:nSatTot*2); residuals_fixed(nSatTot*2+1:end);residuals_float(1:nSatTot*2); residuals_float(nSatTot*2+1:end);outliers(1:nSatTot*2);outliers(nSatTot*2+1:end)], 'double');
                        fwrite(fid_trp, [apriori_ZHD; STDs], 'double');
                        
                        if (flag_plotproc)
                            if (mode_user == 1 && t == 2), w_bar.shiftDown(); end
                            if (flag_cov == 0)
                                if (flag_ge == 1), rtplot_googleearth (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), date_R(t,:)), end
                                rtplot_matlab (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath(), flag_amb);
                            else
                                if (flag_ge == 1), rtplot_googleearth_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(t,:)), end
                                rtplot_matlab_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath(), flag_amb);
                            end
                            if (flag_amb == 1)
                                rtplot_amb (t, window, Xhat_t_t(o3+1:o3+nSatTot), sqrt(diag(Cee(o3+1:o3+nSatTot,o3+1:o3+nSatTot))), conf_cs);
                                pause(0.1);
                            else
                                if (flag_skyplot == 1)
                                    rtplot_skyplot (t, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                    rtplot_snr (snr_R(:,t), Eph_t, SP3);
                                else
                                    rttext_sat (t, azR, elR, snr_R(:,t), conf_sat, pivot, Eph_t, SP3);
                                end
                                pause(0.01);
                            end
                        end
                        
                        w_bar.goTime(t);
                    end
                    
                    w_bar.close();
                    
                    fclose(fid_kal);
                    fclose(fid_sat);
                    fclose(fid_dop);
                    fclose(fid_conf);
                    fclose(fid_res);
                    fclose(fid_trp);
                    
                case goGNSS.MODE_PP_LS_C_DD
                    
                    %----------------------------------------------------------------------------------------------
                    %% POST-PROCESSING (RELATIVE POSITIONING): LEAST SQUARES ON CODE DOUBLE DIFFERENCES
                    %----------------------------------------------------------------------------------------------
                    
                    fid_kal = fopen([filerootOUT '_kal_000.bin'],'w+');
                    fid_sat = fopen([filerootOUT '_sat_000.bin'],'w+');
                    fid_dop = fopen([filerootOUT '_dop_000.bin'],'w+');
                    fid_conf = fopen([filerootOUT '_conf_000.bin'],'w+');
                    fid_res = fopen([filerootOUT '_res_000.bin'],'w+');
                    
                    nN = nSatTot;
                    nT = 0;
                    check_on = 0;
                    check_off = 0;
                    check_pivot = 0;
                    check_cs = 0;
                    
                    plot_t = 1;
                    
                    % goGPS waiting bar
                    w_bar.setBarLen(length(time_GPS));
                    w_bar.createNewBar('Processing...');
                    
                    for t = 1 : length(time_GPS)
                        
                        Eph_t = rt_find_eph (Eph, time_GPS_diff(t), nSatTot);
                        
                        goGPS_LS_DD_code(time_GPS_diff(t), pos_R, pos_M(:,t), pr1_R(:,t), pr1_M(:,t), pr2_R(:,t), pr2_M(:,t), snr_R(:,t), snr_M(:,t), Eph_t, SP3, iono, lambda, frequencies(1), state.getProcessingRate());
                        
                        if (t == 1)
                            fwrite(fid_sat, nSatTot, 'int8');
                            fwrite(fid_conf, nSatTot, 'int8');
                            fwrite(fid_res, nSatTot, 'int8');
                        end
                        
                        if ~isempty(Xhat_t_t) && ~any(isnan([Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)]))
                            Xhat_t_t_dummy = [Xhat_t_t; zeros(nN,1)];
                            Cee_dummy = [Cee zeros(o3,nN); zeros(nN,o3) zeros(nN,nN)];
                            fwrite(fid_kal, [Xhat_t_t_dummy; Cee_dummy(:)], 'double');
                            fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                            fwrite(fid_dop, [PDOP; HDOP; VDOP; 0; 0; 0], 'double');
                            fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
                            residuals_dummy = NaN(1,nSatTot);
                            fwrite(fid_res, [residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'], 'double');
                            
                            if (flag_plotproc)
                                if (mode_user == 1 && t == 1), w_bar.shiftDown(); end
                                if (flag_cov == 0)
                                    if (flag_ge == 1), rtplot_googleearth (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), date_R(t,:)), end
                                    rtplot_matlab (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                                else
                                    if (flag_ge == 1), rtplot_googleearth_cov (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(t,:)), end
                                    rtplot_matlab_cov (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                                end
                                if (flag_skyplot == 1)
                                    rtplot_skyplot (plot_t, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                    rtplot_snr (snr_R(:,t), Eph_t, SP3);
                                else
                                    rttext_sat (plot_t, azR, elR, snr_R(:,t), conf_sat, pivot, Eph_t, SP3);
                                end
                                plot_t = plot_t + 1;
                                pause(0.01);
                            end
                            if (state.isModeBlock())
                                if (~isempty(A_epo))
                                    n_obs_epoch(t) = length(y0_epo);
                                    y0_all(epoch_track+1:epoch_track+n_obs_epoch(t)) = y0_epo;
                                    b_all( epoch_track+1:epoch_track+n_obs_epoch(t)) =  b_epo;
                                    A_all( epoch_track+1:epoch_track+n_obs_epoch(t),1:3) = A_epo(:,1:3);  %#ok<SPRIX>
                                    Q_all( epoch_track+1:epoch_track+n_obs_epoch(t),epoch_track+1:epoch_track+n_obs_epoch(t)) = Q_epo;  %#ok<SPRIX>
                                    satpr_track(:,t) = 0;
                                    satpr_track(conf_sat==-1 | conf_sat==+1,t) = 1;
                                    sat_track(epoch_track+1:epoch_track+n_obs_epoch(t),3) = -1;
                                    sat_track(epoch_track+1:epoch_track+n_obs_epoch(t),2) = find(satpr_track(:,t) == 1);
                                    sat_track(epoch_track+1:epoch_track+n_obs_epoch(t),1) = t;
                                    pivot_track(t) = pivot;
                                    epoch_index(t) = epoch_track + n_obs_epoch(t);
                                    epoch_track = epoch_index(t);
                                end
                            end
                        else
                            unused_epochs(t) = 1;
                        end
                        
                        w_bar.goTime(t);
                    end
                    
                    w_bar.close();
                    
                    fclose(fid_kal);
                    fclose(fid_sat);
                    fclose(fid_dop);
                    fclose(fid_conf);
                    fclose(fid_res);
                    
                case goGNSS.MODE_PP_KF_C_DD
                    
                    %----------------------------------------------------------------------------------------------
                    %% POST-PROCESSING (RELATIVE POSITIONING): KALMAN FILTER ON CODE DOUBLE DIFFERENCES
                    %----------------------------------------------------------------------------------------------
                    
                    fid_kal = fopen([filerootOUT '_kal_000.bin'],'w+');
                    fid_sat = fopen([filerootOUT '_sat_000.bin'],'w+');
                    fid_dop = fopen([filerootOUT '_dop_000.bin'],'w+');
                    fid_conf = fopen([filerootOUT '_conf_000.bin'],'w+');
                    fid_res = fopen([filerootOUT '_res_000.bin'],'w+');
                    
                    nN = nSatTot;
                    nT = 0;
                    check_on = 0;
                    check_off = 0;
                    check_pivot = 0;
                    check_cs = 0;
                    
                    if ~state.isSeamlessKF()
                        kalman_initialized = 0;
                    end
                    if (~kalman_initialized)
                        while (~kalman_initialized)
                            if (isempty(time_GPS))
                                fprintf('It was not possible to initialize the Kalman filter.\n');
                                return
                            end
                            
                            Eph_t = rt_find_eph (Eph, time_GPS_diff(1), nSatTot);
                            
                            kalman_initialized = goGPS_KF_DD_code_init(pos_R, pos_M(:,1), time_GPS_diff(1), pr1_R(:,1), pr1_M(:,1), pr2_R(:,1), pr2_M(:,1), snr_R(:,1), snr_M(:,1), Eph_t, SP3, iono, lambda, frequencies(1), state.getProcessingRate());
                            
                            if (~kalman_initialized)
                                pos_M(:,1) = []; time_GPS_diff(1) = []; time_GPS(1) = []; week_R(1) = [];
                                pr1_R(:,1) = []; pr1_M(:,1) = []; ph1_R(:,1) = []; ph1_M(:,1) = []; dop1_R(:,1) = []; dop1_M(:,1) = [];
                                pr2_R(:,1) = []; pr2_M(:,1) = []; ph2_R(:,1) = []; ph2_M(:,1) = []; dop2_R(:,1) = []; dop2_M(:,1) = [];
                                snr_R(:,1) = []; snr_M(:,1) = []; dtMdot(1) = [];
                            end
                        end
                        
                        Xhat_t_t_dummy = [Xhat_t_t; zeros(nN,1)];
                        Cee_dummy = [Cee zeros(o3,nN); zeros(nN,o3) zeros(nN,nN)];
                        fwrite(fid_kal, [Xhat_t_t_dummy; Cee_dummy(:)], 'double');
                        fwrite(fid_sat, nSatTot, 'int8');
                        fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                        fwrite(fid_dop, [PDOP; HDOP; VDOP; KPDOP; KHDOP; KVDOP], 'double');
                        fwrite(fid_conf, [nSatTot; conf_sat; conf_cs; pivot], 'int8');
                        fwrite(fid_res, nSatTot, 'int8');
                        fwrite(fid_res, [residuals_fixed(1:nSatTot); residuals_fixed(nSatTot+1:end); residuals_float(1:nSatTot); residuals_float(nSatTot+1:end); outliers(1:nSatTot); outliers(nSatTot+1:end)], 'double');
                        
                        if (flag_plotproc)
                            if (flag_cov == 0)
                                if (flag_ge == 1), rtplot_googleearth (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), date_R(1,:)), end
                                rtplot_matlab (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), 0, 0, 0, 0, flag_ms, gs.getReferencePath());
                            else
                                if (flag_ge == 1), rtplot_googleearth_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(1,:)), end
                                rtplot_matlab_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), 0, 0, 0, 0, flag_ms, gs.getReferencePath());
                            end
                            if (flag_skyplot == 1)
                                rtplot_skyplot (1, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                rtplot_snr (snr_R(:,1), Eph_t, SP3);
                            else
                                rttext_sat (1, azR, elR, snr_R(:,1), conf_sat, pivot, Eph_t, SP3);
                            end
                        end
                        t1 = 2;
                    else
                        fwrite(fid_sat, nSatTot, 'int8');
                        fwrite(fid_conf, nSatTot, 'int8');
                        fwrite(fid_res, nSatTot, 'int8');
                        t1 = 1;
                    end
                    
                    % goGPS waiting bar
                    w_bar.setBarLen(length(time_GPS));
                    w_bar.createNewBar('Processing...');
                    
                    for t = t1 : length(time_GPS)
                        
                        Eph_t = rt_find_eph (Eph, time_GPS_diff(t), nSatTot);
                        
                        [check_on, check_off, check_pivot, check_cs] = goGPS_KF_DD_code_loop(pos_M(:,t), time_GPS_diff(t), pr1_R(:,t), pr1_M(:,t), pr2_R(:,t), pr2_M(:,t), snr_R(:,t), snr_M(:,t), Eph_t, SP3, iono, lambda, frequencies(1), state.getProcessingRate());
                        
                        Xhat_t_t_dummy = [Xhat_t_t; zeros(nN,1)];
                        Cee_dummy = [Cee zeros(o3,nN); zeros(nN,o3) zeros(nN,nN)];
                        fwrite(fid_kal, [Xhat_t_t_dummy; Cee_dummy(:)], 'double');
                        fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                        fwrite(fid_dop, [PDOP; HDOP; VDOP; KPDOP; KHDOP; KVDOP], 'double');
                        fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
                        residuals_dummy = NaN(1,nSatTot);
                        fwrite(fid_res, [residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'], 'double');
                        
                        if (flag_plotproc)
                            if (mode_user == 1 && t == 2), w_bar.shiftDown(); end
                            if (flag_cov == 0)
                                if (flag_ge == 1), rtplot_googleearth (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), date_R(t,:)), end
                                rtplot_matlab (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                            else
                                if (flag_ge == 1), rtplot_googleearth_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(t,:)), end
                                rtplot_matlab_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                            end
                            if (flag_skyplot == 1)
                                rtplot_skyplot (t, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                rtplot_snr (snr_R(:,t), Eph_t, SP3);
                            else
                                rttext_sat (t, azR, elR, snr_R(:,t), conf_sat, pivot, Eph_t, SP3);
                            end
                            pause(0.01);
                        end
                        
                        w_bar.goTime(t);
                    end
                    
                    w_bar.close();
                    
                    fclose(fid_kal);
                    fclose(fid_sat);
                    fclose(fid_dop);
                    fclose(fid_conf);
                    fclose(fid_res);
                    
                case goGNSS.MODE_PP_LS_CP_DD_L
                    
                    %----------------------------------------------------------------------------------------------
                    %% POST-PROCESSING (RELATIVE POSITIONING): LEAST SQUARES ON CODE AND PHASE DOUBLE DIFFERENCES
                    %----------------------------------------------------------------------------------------------
                    
                    fid_kal = fopen([filerootOUT '_kal_000.bin'],'w+');
                    fid_sat = fopen([filerootOUT '_sat_000.bin'],'w+');
                    fid_dop = fopen([filerootOUT '_dop_000.bin'],'w+');
                    fid_conf = fopen([filerootOUT '_conf_000.bin'],'w+');
                    fid_res = fopen([filerootOUT '_res_000.bin'],'w+');
                    
                    nN = nSatTot;
                    nT = 0;
                    check_on = 0;
                    check_off = 0;
                    check_pivot = 0;
                    check_cs = 0;
                    
                    plot_t = 1;
                    
                    % goGPS waiting bar
                    w_bar.setBarLen(length(time_GPS_diff));
                    w_bar.createNewBar('Processing...');
                    
                    for t = 1 : length(time_GPS_diff)
                        
                        Eph_t = rt_find_eph (Eph, time_GPS_diff(t), nSatTot);
                        
                        goGPS_LS_DD_code_phase(time_GPS_diff(t), pos_M(:,t), pr1_R(:,t), pr1_M(:,t), pr2_R(:,t), pr2_M(:,t), ph1_R(:,t), ph1_M(:,t), ph2_R(:,t), ph2_M(:,t), snr_R(:,t), snr_M(:,t), Eph_t, SP3, iono, lambda, frequencies(1), state.getProcessingRate(), flag_IAR);
                        
                        if ~isempty(Xhat_t_t) && ~any(isnan([Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)]))
                            fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
                            if (t == 1)
                                fwrite(fid_sat, nSatTot, 'int8');
                                fwrite(fid_conf, nSatTot, 'int8');
                                fwrite(fid_res, nSatTot, 'int8');
                            end
                            fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                            fwrite(fid_dop, [PDOP; HDOP; VDOP; 0; 0; 0], 'double');
                            fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
                            residuals_dummy = NaN(1,nSatTot);
                            fwrite(fid_res, [residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'], 'double');
                            
                            if (flag_plotproc)
                                if (mode_user == 1 && t == 1), w_bar.shiftDown(); end
                                if (flag_cov == 0)
                                    if (flag_ge == 1), rtplot_googleearth (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), date_R(t,:)), end
                                    rtplot_matlab (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                                else
                                    if (flag_ge == 1), rtplot_googleearth_cov (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(t,:)), end
                                    rtplot_matlab_cov (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                                end
                                if (flag_skyplot == 1)
                                    rtplot_skyplot (plot_t, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                    rtplot_snr (snr_R(:,t), Eph_t, SP3);
                                else
                                    rttext_sat (plot_t, azR, elR, snr_R(:,t), conf_sat, pivot, Eph_t, SP3);
                                end
                                plot_t = plot_t + 1;
                                pause(0.01);
                            end
                        else
                            unused_epochs(t) = 1;
                        end
                        
                        w_bar.goTime(t);
                    end
                    
                    w_bar.close();
                    
                    fclose(fid_kal);
                    fclose(fid_sat);
                    fclose(fid_dop);
                    fclose(fid_conf);
                    fclose(fid_res);
                    
                case goGNSS.MODE_PP_LS_CP_DD_MR
                    
                    %------------------------------------------------------------------------------------------------------------
                    %% POST-PROCESSING (RELATIVE POSITIONING): LEAST SQUARES ON CODE AND PHASE DOUBLE DIFFERENCES (MULTI-RECEIVER)
                    %------------------------------------------------------------------------------------------------------------
                    
                    fid_kal = fopen([filerootOUT '_kal_000.bin'],'w+');
                    fid_sat = fopen([filerootOUT '_sat_000.bin'],'w+');
                    fid_dop = fopen([filerootOUT '_dop_000.bin'],'w+');
                    fid_conf = fopen([filerootOUT '_conf_000.bin'],'w+');
                    fid_res = fopen([filerootOUT '_res_000.bin'],'w+');
                    
                    nN = nSatTot;
                    nT = 0;
                    check_on = 0;
                    check_off = 0;
                    check_pivot = 0;
                    check_cs = 0;
                    
                    plot_t = 1;
                    
                    % goGPS waiting bar
                    w_bar.setBarLen(length(time_GPS_diff));
                    w_bar.createNewBar('Processing...');
                    
                    for t = 1 : length(time_GPS_diff)
                        
                        Eph_t = rt_find_eph (Eph, time_GPS_diff(t), nSatTot);
                        
                        goGPS_LS_DD_code_phase_MR(time_GPS_diff(t), multi_antenna_rf, pos_M(:,t), squeeze(pr1_R(:,t,:)), pr1_M(:,t), squeeze(pr2_R(:,t,:)), pr2_M(:,t), squeeze(ph1_R(:,t,:)), ph1_M(:,t), squeeze(ph2_R(:,t,:)), ph2_M(:,t), squeeze(snr_R(:,t,:)), snr_M(:,t), Eph_t, SP3, iono, lambda, frequencies(1), state.getProcessingRate(), state.getProcessingRate(), flag_IAR);
                        
                        if ~isempty(Xhat_t_t) && ~any(isnan([Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)]))
                            fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
                            if (t == 1)
                                fwrite(fid_sat, nSatTot, 'int8');
                                fwrite(fid_conf, nSatTot, 'int8');
                                fwrite(fid_res, nSatTot, 'int8');
                            end
                            fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                            fwrite(fid_dop, [PDOP; HDOP; VDOP; 0; 0; 0], 'double');
                            fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
                            residuals_dummy = NaN(1,nSatTot);
                            fwrite(fid_res, [residuals_dummy; residuals_dummy; residuals_dummy; residuals_dummy; residuals_dummy; residuals_dummy], 'double');
                            
                            if (flag_plotproc)
                                if (mode_user == 1 && t == 1), w_bar.shiftDown(); end
                                if (flag_cov == 0)
                                    if (flag_ge == 1), rtplot_googleearth (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), date_R(t,:)), end
                                    rtplot_matlab (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                                else
                                    if (flag_ge == 1), rtplot_googleearth_cov (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(t,:)), end
                                    rtplot_matlab_cov (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                                end
                                if (flag_skyplot == 1)
                                    rtplot_skyplot (plot_t, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                    rtplot_snr (snr_R(:,t), Eph_t, SP3);
                                else
                                    rttext_sat (plot_t, azR, elR, snr_R(:,t), conf_sat, pivot, Eph_t, SP3);
                                end
                                plot_t = plot_t + 1;
                                pause(0.01);
                            end
                        else
                            unused_epochs(t) = 1;
                        end
                        
                        w_bar.goTime(t);
                    end
                    w_bar.close();
                    
                    fclose(fid_kal);
                    fclose(fid_sat);
                    fclose(fid_dop);
                    fclose(fid_conf);
                    fclose(fid_res);
                    
                case goGNSS.MODE_PP_LS_C_SA_MR
                    
                    %----------------------------------------------------------------------------------------------
                    %% POST-PROCESSING (ABSOLUTE POSITIONING): LEAST SQUARES ON CODE (MULTI-RECEIVER, AVERAGE)
                    %----------------------------------------------------------------------------------------------
                    
                    fid_kal = fopen([filerootOUT '_kal_000.bin'],'w+');
                    fid_sat = fopen([filerootOUT '_sat_000.bin'],'w+');
                    fid_dop = fopen([filerootOUT '_dop_000.bin'],'w+');
                    fid_conf = fopen([filerootOUT '_conf_000.bin'],'w+');
                    fid_res = fopen([filerootOUT '_res_000.bin'],'w+');
                    
                    nN = nSatTot;
                    nT = 0;
                    check_on = 0;
                    check_off = 0;
                    check_pivot = 0;
                    check_cs = 0;
                    
                    plot_t = 1;
                    
                    % goGPS waiting bar
                    w_bar.setBarLen(length(time_GPS_diff));
                    w_bar.createNewBar('Processing...');
                    
                    for t = 1 : length(time_GPS_diff)
                        
                        Eph_t = rt_find_eph (Eph, time_GPS_diff(t), nSatTot);
                        
                        sbas_t = find_sbas(sbas, t);
                        
                        goGPS_LS_SA_code_MR(time_GPS_diff(t), squeeze(pr1_R(:,t,:)), squeeze(pr2_R(:,t,:)), squeeze(snr_R(:,t,:)), Eph_t, SP3, iono, sbas_t, lambda, frequencies(1), state.getProcessingRate());
                        
                        if (t == 1)
                            fwrite(fid_sat, nSatTot, 'int8');
                            fwrite(fid_conf, nSatTot, 'int8');
                            fwrite(fid_res, nSatTot, 'int8');
                        end
                        
                        if ~isempty(Xhat_t_t) && ~any(isnan([Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)]))
                            Xhat_t_t_dummy = [Xhat_t_t; zeros(nN,1)];
                            Cee_dummy = [Cee zeros(o3,nN); zeros(nN,o3) zeros(nN,nN)];
                            fwrite(fid_kal, [Xhat_t_t_dummy; Cee_dummy(:)], 'double');
                            fwrite(fid_sat, [zeros(nSatTot,1); azR(:,1); zeros(nSatTot,1); elR(:,1); zeros(nSatTot,1); distR(:,1)], 'double');
                            fwrite(fid_dop, [PDOP; HDOP; VDOP; 0; 0; 0], 'double');
                            fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
                            residuals_dummy = NaN(1,nSatTot);
                            fwrite(fid_res, [residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'; residuals_dummy'], 'double');
                            
                            if (flag_plotproc)
                                if (mode_user == 1 && t == 1), w_bar.shiftDown(); end
                                if (flag_cov == 0)
                                    if (flag_ge == 1), rtplot_googleearth (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), date_R(t,:)), end
                                    rtplot_matlab (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                                else
                                    if (flag_ge == 1), rtplot_googleearth_cov (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(t,:)), end
                                    rtplot_matlab_cov (plot_t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], zeros(3,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath());
                                end
                                if (flag_skyplot == 1)
                                    rtplot_skyplot (plot_t, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                    rtplot_snr (snr_R(:,t), Eph_t, SP3);
                                else
                                    rttext_sat (plot_t, azR, elR, snr_R(:,t), conf_sat, pivot, Eph_t, SP3);
                                end
                                plot_t = plot_t + 1;
                                pause(0.01);
                            end
                        else
                            unused_epochs(t) = 1;
                        end
                        
                        w_bar.goTime(t);
                    end
                    w_bar.close();
                    
                    fclose(fid_kal);
                    fclose(fid_sat);
                    fclose(fid_dop);
                    fclose(fid_conf);
                    fclose(fid_res);
                    
                case goGNSS.MODE_PP_BLK_CP_DD_STATIC
                    %----------------------------------------------------------------------------------------------
                    %% GO BLOCK) POST-PROCESSING (RELATIVE POSITIONING): BLOCK ON CODE AND PHASE DOUBLE DIFFERENCES
                    %----------------------------------------------------------------------------------------------
                    
                    % set a priori coordinates
                    if (exist('pos_R_crd','var') && any(pos_R_crd))
                        pos_R = pos_R_crd;
                    else
                        pos_R = median (pos_R_new(:,:,1),2);
                    end
                    
                    if state.isSeamlessHR()
                        go_block = Core_Block_Old (numel(time_GPS), sum(serialize(pr1_R(:,:,1) ~= 0)), sum(serialize(ph1_R(:,:,1) ~= 0)));
                        go_block.prepare (time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV);
                        go_block.solve(state.getSolutionRate());
                    else
                        
                        %-----------------------------------------------------------------------------------------------------------
                        %% GO BLOCK HR) POST-PROCESSING (RELATIVE POSITIONING): BLOCK ON CODE AND PHASE DOUBLE DIFFERENCES HIGH RATE
                        %-----------------------------------------------------------------------------------------------------------
                        go_block = Core_Block_Old.goMultiHighRate(time_GPS_diff, pos_R, pos_M, pr1_R, pr1_M, pr2_R, pr2_M, ph1_R, ph1_M, ph2_R, ph2_M, snr_R, snr_M,  Eph, SP3, iono, lambda, antenna_PCV, state.getSolutionRate());
                        unused_epochs = [];
                    end
                    
                    if state.isOutBlockObj()
                        block_struct = go_block.struct();
                        save([filerootOUT '_block_obj.mat'], 'block_struct');
                    end
                case goGNSS.MODE_PP_KF_CP_DD
                    if (mode_vinc == 0)
                        %--------------------------------------------------------------------------------------------------------------------
                        %% POST-PROCESSING (RELATIVE POSITIONING): KALMAN FILTER ON CODE AND PHASE DOUBLE DIFFERENCES WITHOUT LINE CONSTRAINT
                        %--------------------------------------------------------------------------------------------------------------------
                        
                        if (flag_var_dyn_model == 0)
                            
                            fid_kal = fopen([filerootOUT '_kal_000.bin'],'w+');
                            fid_sat = fopen([filerootOUT '_sat_000.bin'],'w+');
                            fid_dop = fopen([filerootOUT '_dop_000.bin'],'w+');
                            fid_conf = fopen([filerootOUT '_conf_000.bin'],'w+');
                            fid_res = fopen([filerootOUT '_res_000.bin'],'w+');
                            
                            % apply SBAS corrections on code
                            if ~isempty(sbas)
                                sbas_prc=sbas.prc';
                                index = find(pr1_R~=0);
                                pr1_R(index) = pr1_R(index) + sbas_prc(index);
                                
                                index = find(pr1_M~=0);
                                pr1_M(index) = pr1_M(index) + sbas_prc(index);
                                clear sbas_prc
                            end
                            
                            if ~state.isSeamlessKF()
                                kalman_initialized = 0;
                            end
                            if (~kalman_initialized)
                                if state.getForwardBackwardKF() < 0
                                    id_init = numel(time_GPS_diff);
                                else
                                    id_init = 1;
                                end
                                
                                while (~kalman_initialized)
                                    if (numel(time_GPS_diff) < id_init)
                                        fprintf('It was not possible to initialize the Kalman filter.\n');
                                        return
                                    end
                                    
                                    Eph_t = rt_find_eph (Eph, time_GPS_diff(id_init), nSatTot);
                                    
                                    sbas_t = find_sbas(sbas, 1);
                                    
                                    kalman_initialized = goGPS_KF_DD_code_phase_init(pos_R, pos_M(:,id_init), time_GPS_diff(id_init), pr1_R(:,id_init), pr1_M(:,id_init), ph1_R(:,id_init), ph1_M(:,id_init), dop1_R(:,id_init), dop1_M(:,id_init), pr2_R(:,id_init), pr2_M(:,id_init), ph2_R(:,id_init), ph2_M(:,id_init), dop2_R(:,id_init), dop2_M(:,id_init), snr_R(:,id_init), snr_M(:,id_init), Eph_t, SP3, iono, lambda, frequencies, state.getProcessingRate(), dtMdot(id_init), flag_IAR, flag_XR, flag_tropo, sbas_t);
                                    
                                    if (~kalman_initialized)
                                        pos_M(:,id_init) = []; time_GPS_diff(id_init) = []; time_GPS(id_init) = []; week_R(id_init) = [];
                                        pr1_R(:,id_init) = []; pr1_M(:,id_init) = []; ph1_R(:,id_init) = []; ph1_M(:,id_init) = []; dop1_R(:,id_init) = []; dop1_M(:,id_init) = [];
                                        pr2_R(:,id_init) = []; pr2_M(:,id_init) = []; ph2_R(:,id_init) = []; ph2_M(:,id_init) = []; dop2_R(:,id_init) = []; dop2_M(:,id_init) = [];
                                        snr_R(:,id_init) = []; snr_M(:,id_init) = []; dtMdot(id_init) = [];
                                    end
                                end
                                
                                fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
                                fwrite(fid_sat, nSatTot, 'int8');
                                fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                                fwrite(fid_dop, [PDOP; HDOP; VDOP; KPDOP; KHDOP; KVDOP], 'double');
                                fwrite(fid_conf, [nSatTot; conf_sat; conf_cs; pivot], 'int8');
                                fwrite(fid_res, nSatTot, 'int8');
                                fwrite(fid_res, [residuals_fixed(1:nSatTot*2); residuals_fixed(nSatTot*2+1:end);residuals_float(1:nSatTot*2); residuals_float(nSatTot*2+1:end);outliers(1:nSatTot*2);outliers(nSatTot*2+1:end)], 'double');
                                
                                if (flag_plotproc)
                                    if (flag_cov == 0)
                                        if (flag_ge == 1), rtplot_googleearth (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), date_R(1,:)), end
                                        rtplot_matlab (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), 0, 0, 0, 0, flag_ms, gs.getReferencePath(), flag_amb);
                                    else
                                        if (flag_ge == 1), rtplot_googleearth_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(1,:)), end
                                        rtplot_matlab_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), 0, 0, 0, 0, flag_ms, gs.getReferencePath(), flag_amb);
                                    end
                                    if (flag_amb == 1)
                                        rtplot_amb (1, window, Xhat_t_t(o3+1:o3+nSatTot), sqrt(diag(Cee(o3+1:o3+nSatTot,o3+1:o3+nSatTot))), conf_cs)
                                    else
                                        if (flag_skyplot == 1)
                                            rtplot_skyplot (1, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                            rtplot_snr (snr_R(:,1), Eph_t, SP3);
                                        else
                                            rttext_sat (1, azR, elR, snr_R(:,1), conf_sat, pivot, Eph_t, SP3);
                                        end
                                    end
                                end
                                t1 = 2;
                            else
                                fwrite(fid_sat, nSatTot, 'int8');
                                fwrite(fid_conf, nSatTot, 'int8');
                                fwrite(fid_res, nSatTot, 'int8');
                                t1 = 1;
                            end
                            
                            % goGPS waiting bar
                            w_bar.setBarLen((1+state.isForwardBackwardKF()) * length(time_GPS_diff));
                            if state.getForwardBackwardKF() < 0
                                w_bar.createNewBar('Processing backward...');
                            else
                                w_bar.createNewBar('Processing forward...');
                            end
                            
                            % forward - backward filter
                            switch state.getForwardBackwardKF()
                                case -1, time_steps = [((length(time_GPS_diff) - t1 + 1) : -1 : 1) (1 : length(time_GPS_diff))];
                                case  1, time_steps = [(t1 : length(time_GPS_diff)) ((length(time_GPS_diff)) : -1 : 1)];
                                case  0, time_steps = t1 : length(time_GPS_diff);
                            end
                            t0 = 0;
                            is_forward = state.getForwardBackwardKF() > -1;
                            for t = time_steps
                                if t == t0
                                    if state.getForwardBackwardKF() > 0
                                        w_bar.createNewBar('Processing backward...');
                                    else
                                        w_bar.createNewBar('Processing forward...');
                                    end
                                    is_forward = ~is_forward;
                                end
                                t0 = t;
                                residuals_fixed=NaN(4*nSatTot,1);
                                residuals_float=NaN(4*nSatTot,1);
                                outliers=zeros(4*nSatTot,1);
                                
                                Eph_t = rt_find_eph (Eph, time_GPS_diff(t), nSatTot);
                                
                                sbas_t = find_sbas(sbas, t);
                                
                                [check_on, check_off, check_pivot, check_cs] = goGPS_KF_DD_code_phase_loop(pos_M(:,t), time_GPS_diff(t), pr1_R(:,t), pr1_M(:,t), ph1_R(:,t), ph1_M(:,t), dop1_R(:,t), dop1_M(:,t), pr2_R(:,t), pr2_M(:,t), ph2_R(:,t), ph2_M(:,t), dop2_R(:,t), dop2_M(:,t), snr_R(:,t), snr_M(:,t), Eph_t, SP3, iono, lambda, frequencies, obs_comb, state.getProcessingRate(), dtMdot(t), flag_IAR, flag_tropo, antenna_PCV, sbas_t);
                                
                                fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
                                fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                                fwrite(fid_dop, [PDOP; HDOP; VDOP; KPDOP; KHDOP; KVDOP], 'double');
                                fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
                                fwrite(fid_res, [residuals_fixed(1:nSatTot*2); residuals_fixed(nSatTot*2+1:end);residuals_float(1:nSatTot*2); residuals_float(nSatTot*2+1:end);outliers(1:nSatTot*2);outliers(nSatTot*2+1:end)], 'double');
                                
                                if (flag_plotproc)
                                    if (mode_user == 1 && t == 2), w_bar.shiftDown(); end
                                    if (flag_cov == 0)
                                        if (flag_ge == 1), rtplot_googleearth (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), date_R(t,:)), end
                                        rtplot_matlab (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath(), flag_amb);
                                    else
                                        if (flag_ge == 1), rtplot_googleearth_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(t,:)), end
                                        rtplot_matlab_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath(), flag_amb);
                                    end
                                    if (flag_amb == 1)
                                        rtplot_amb (t, window, Xhat_t_t(o3+1:o3+nSatTot), sqrt(diag(Cee(o3+1:o3+nSatTot,o3+1:o3+nSatTot))), conf_cs);
                                        pause(0.1);
                                    else
                                        if (flag_skyplot == 1)
                                            rtplot_skyplot (t, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                            rtplot_snr (snr_R(:,t), Eph_t, SP3);
                                        else
                                            rttext_sat (t, azR, elR, snr_R(:,t), conf_sat, pivot, Eph_t, SP3);
                                        end
                                        pause(0.01);
                                    end
                                end
                                switch state.getForwardBackwardKF()
                                    case -1, w_bar.goTime(~is_forward * (length(time_GPS_diff) - t + 1) + is_forward * (length(time_GPS_diff) + t));
                                    case  1, w_bar.goTime(~is_forward * (length(time_GPS_diff) * 2 - t + 1) + is_forward * t);
                                    case  0, w_bar.goTime(t);
                                end
                            end
                            w_bar.close();
                            
                            fclose(fid_kal);
                            fclose(fid_sat);
                            fclose(fid_dop);
                            fclose(fid_conf);
                            fclose(fid_res);
                        else
                            
                            fid_dyn = fopen([state.obs_dir '\' state.obs_name{1}(1:end-4) '_dyn.bin'],'r+');
                            fid_kal = fopen([filerootOUT '_kal_000.bin'],'w+');
                            fid_sat = fopen([filerootOUT '_sat_000.bin'],'w+');
                            fid_dop = fopen([filerootOUT '_dop_000.bin'],'w+');
                            fid_conf = fopen([filerootOUT '_conf_000.bin'],'w+');
                            fid_res = fopen([filerootOUT '_res_000.bin'],'w+');
                            
                            if ~state.isSeamlessKF()
                                kalman_initialized = 0;
                            end
                            if (~kalman_initialized)
                                while (~kalman_initialized)
                                    if (isempty(time_GPS_diff))
                                        fprintf('It was not possible to initialize the Kalman filter.\n');
                                        return
                                    end
                                    
                                    Eph_t = rt_find_eph (Eph, time_GPS_diff(1), nSatTot);
                                    
                                    flag_dyn = 1;
                                    order = fread(fid_dyn,1,'uint8');
                                    
                                    kalman_initialized = goGPS_KF_DD_code_phase_init_model(pos_R, pos_M(:,1), time_GPS_diff(1), pr1_R(:,1), pr1_M(:,1), ph1_R(:,1), ph1_M(:,1), dop1_R(:,1), dop1_M(:,1), pr2_R(:,1), pr2_M(:,1), ph2_R(:,1), ph2_M(:,1), dop2_R(:,1), dop2_M(:,1), snr_R(:,1), snr_M(:,1), Eph_t, SP3, iono, lambda, order, frequencies, state.getProcessingRate(), dtMdot(1), flag_IAR);
                                    
                                    if (~kalman_initialized)
                                        pos_M(:,1) = []; time_GPS_diff(1) = []; time_GPS(1) = []; week_R(1) = [];
                                        pr1_R(:,1) = []; pr1_M(:,1) = []; ph1_R(:,1) = []; ph1_M(:,1) = []; dop1_R(:,1) = []; dop1_M(:,1) = [];
                                        pr2_R(:,1) = []; pr2_M(:,1) = []; ph2_R(:,1) = []; ph2_M(:,1) = []; dop2_R(:,1) = []; dop2_M(:,1) = [];
                                        snr_R(:,1) = []; snr_M(:,1) = []; dtMdot(1) = [];
                                    end
                                end
                                
                                if (flag_stopGOstop == 1)
                                    index = 1;
                                    X_init = Xhat_t_t([1 o1+1 o2+1]);
                                    X_ENU = global2localPos(Xhat_t_t([1 o1+1 o2+1]), X_init);
                                    E0(index,1) = X_ENU(1,1); %#ok<SAGROW>
                                    N0(index,1) = X_ENU(2,1); %#ok<SAGROW>
                                    Cee_ENU = global2localCov(Cee([1 o1+1 o2+1],[1 o1+1 o2+1],:), X_init);
                                    sigmaq_E0(index,1) = Cee_ENU(1,1); %#ok<SAGROW>
                                    sigmaq_N0(index,1) = Cee_ENU(2,2); %#ok<SAGROW>
                                    sigma_EN0(index,1) = Cee_ENU(1,2); %#ok<SAGROW>
                                    mDIR = 0; qDIR = 0; angleDIR = 0;
                                    sigma_angleDIR = 0;
                                    P1 = [E0(index), N0(index)]; P2 = P1;
                                    P1_ENU = [P1(1); P1(2); X_ENU(3,1)];
                                    P2_ENU = [P2(1); P2(2); X_ENU(3,1)];
                                    P1_GLB = local2globalPos(P1_ENU, X_init);
                                    P2_GLB = local2globalPos(P2_ENU, X_init);
                                    [P1_UTM_E, P1_UTM_N] = cart2plan(P1_GLB(1), P1_GLB(2), P1_GLB(3));
                                    [P2_UTM_E, P2_UTM_N] = cart2plan(P2_GLB(1), P2_GLB(2), P2_GLB(3));
                                end
                                
                                fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
                                fwrite(fid_sat, nSatTot, 'int8');
                                fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                                fwrite(fid_dop, [PDOP; HDOP; VDOP; KPDOP; KHDOP; KVDOP], 'double');
                                fwrite(fid_conf, [nSatTot; conf_sat; conf_cs; pivot], 'int8');
                                fwrite(fid_res, nSatTot, 'int8');
                                fwrite(fid_res, [residuals_fixed(1:nSatTot*2); residuals_fixed(nSatTot*2+1:end);residuals_float(1:nSatTot*2); residuals_float(nSatTot*2+1:end);outliers(1:nSatTot*2);outliers(nSatTot*2+1:end)], 'double');
                                
                                if (flag_plotproc)
                                    if (flag_cov == 0)
                                        if (flag_ge == 1), rtplot_googleearth (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), date_R(1,:)), end
                                        if (flag_stopGOstop == 1)
                                            rtplot_matlab_stopGOstop (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), [P1_UTM_E, P1_UTM_N], [P2_UTM_E, P2_UTM_N], flag_ms, gs.getReferencePath(), flag_dyn, flag_amb);
                                        else
                                            rtplot_matlab (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), 0, 0, 0, 0, flag_ms, gs.getReferencePath(), flag_amb);
                                        end
                                    else
                                        if (flag_ge == 1), rtplot_googleearth_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(1,:)), end
                                        if (flag_stopGOstop == 1)
                                            rtplot_matlab_cov_stopGOstop (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), [P1_UTM_E, P1_UTM_N], [P2_UTM_E, P2_UTM_N], flag_ms, gs.getReferencePath(), flag_dyn, flag_amb);
                                        else
                                            rtplot_matlab_cov (1, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,1), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), 0, 0, 0, 0, flag_ms, gs.getReferencePath(), flag_amb);
                                        end
                                    end
                                    if (flag_amb == 1)
                                        rtplot_amb (1, window, Xhat_t_t(o3+1:o3+nSatTot), sqrt(diag(Cee(o3+1:o3+nSatTot,o3+1:o3+nSatTot))), conf_cs)
                                    else
                                        if (flag_skyplot == 1)
                                            rtplot_skyplot (1, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                            rtplot_snr (snr_R(:,1), Eph_t, SP3);
                                        else
                                            rttext_sat (1, azR, elR, snr_R(:,1), conf_sat, pivot, Eph_t, SP3);
                                        end
                                    end
                                end
                                t1 = 2;
                            else
                                fwrite(fid_sat, nSatTot, 'int8');
                                fwrite(fid_conf, nSatTot, 'int8');
                                fwrite(fid_res, nSatTot, 'int8');
                                t1 = 1;
                            end
                            
                            % goGPS waiting bar
                            w_bar.setBarLen(length(time_GPS_diff));
                            w_bar.createNewBar('Processing...');
                            
                            for t = t1 : length(time_GPS_diff)
                                residuals_fixed=NaN(4*nSatTot,1);
                                residuals_float=NaN(4*nSatTot,1);
                                outliers=zeros(4*nSatTot,1);
                                
                                Eph_t = rt_find_eph (Eph, time_GPS_diff(t), nSatTot);
                                
                                order0 = order;
                                order = fread(fid_dyn,1,'uint8');
                                
                                [check_on, check_off, check_pivot, check_cs] = goGPS_KF_DD_code_phase_loop_model(pos_M(:,t), time_GPS_diff(t), pr1_R(:,t), pr1_M(:,t), ph1_R(:,t), ph1_M(:,t), dop1_R(:,t), dop1_M(:,t), pr2_R(:,t), pr2_M(:,t), ph2_R(:,t), ph2_M(:,t), dop2_R(:,t), dop2_M(:,t), snr_R(:,t), snr_M(:,t), Eph_t, SP3, iono, lambda, order, frequencies, state.getProcessingRate(), dtMdot(t), flag_IAR, antenna_PCV);
                                
                                if (flag_stopGOstop == 1)
                                    if (order > order0)
                                        flag_dyn = 2;
                                    end
                                    if (order < order0)
                                        flag_dyn = 3;
                                        index = index+1;
                                    end
                                    if (order == 2)
                                        index = index+1;
                                    end
                                    X_ENU = global2localPos(Xhat_t_t([1 o1+1 o2+1]), X_init);
                                    E0(index,1) = X_ENU(1,1);
                                    N0(index,1) = X_ENU(2,1);
                                    Cee_ENU = global2localCov(Cee([1 o1+1 o2+1],[1 o1+1 o2+1],:), X_init);
                                    sigmaq_E0(index,1) = Cee_ENU(1,1);
                                    sigmaq_N0(index,1) = Cee_ENU(2,2);
                                    sigma_EN0(index,1) = Cee_ENU(1,2);
                                    if (index == 1)
                                        mDIR = 0; qDIR = 0; angleDIR = 0;
                                        P1 = [E0(index), N0(index)]; P2 = P1;
                                        P1_ENU = [P1(1); P1(2); X_ENU(3,1)];
                                        P2_ENU = [P2(1); P2(2); X_ENU(3,1)];
                                    else
                                        [mDIR, qDIR, sigmaq_mDIR, sigmaq_qDIR] = LSinterp(E0, N0, sigmaq_E0, sigmaq_N0, sigma_EN0);
                                        m1 = -(sigmaq_N0(1)   - mDIR*sigma_EN0(1))   / (mDIR*sigmaq_E0(1)   - sigma_EN0(1));
                                        m2 = -(sigmaq_N0(end) - mDIR*sigma_EN0(end)) / (mDIR*sigmaq_E0(end) - sigma_EN0(end));
                                        X1 = (m1*E0(1)   + qDIR - N0(1))   / (m1-mDIR);
                                        X2 = (m2*E0(end) + qDIR - N0(end)) / (m2-mDIR);
                                        P1 = [X1, mDIR*X1+qDIR ];   % projecting according to error covariance
                                        P2 = [X2, mDIR*X2+qDIR ];
                                        P1_ENU = [P1(1); P1(2); X_ENU(3,1)];
                                        P2_ENU = [P2(1); P2(2); X_ENU(3,1)];
                                        angleDIR = atan2(P2(1)-P1(1),P2(2)-P1(2)) * 180/pi;
                                        sigma_angleDIR = 1/(1+mDIR^2) * sqrt(sigmaq_mDIR) * 180/pi;
                                        % sigma_angleDIR = atan(sqrt(sigmaq_mDIR));
                                    end
                                    P1_GLB = local2globalPos(P1_ENU, X_init);
                                    P2_GLB = local2globalPos(P2_ENU, X_init);
                                    [P1_UTM_E, P1_UTM_N] = cart2plan(P1_GLB(1), P1_GLB(2), P1_GLB(3));
                                    [P2_UTM_E, P2_UTM_N] = cart2plan(P2_GLB(1), P2_GLB(2), P2_GLB(3));
                                    pause(0.05)
                                end
                                
                                fwrite(fid_kal, [Xhat_t_t; Cee(:)], 'double');
                                fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                                fwrite(fid_dop, [PDOP; HDOP; VDOP; KPDOP; KHDOP; KVDOP], 'double');
                                fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
                                fwrite(fid_res, [residuals_fixed(1:nSatTot*2); residuals_fixed(nSatTot*2+1:end);residuals_float(1:nSatTot*2); residuals_float(nSatTot*2+1:end);outliers(1:nSatTot*2);outliers(nSatTot*2+1:end)], 'double');
                                
                                if (flag_plotproc)
                                    if (mode_user == 1 && t == 2), w_bar.shiftDown(); end
                                    if (flag_cov == 0)
                                        if (flag_ge == 1), rtplot_googleearth (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), date_R(t,:)), end
                                        if (flag_stopGOstop == 1)
                                            rtplot_matlab_stopGOstop (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), [P1_UTM_E, P1_UTM_N], [P2_UTM_E, P2_UTM_N], flag_ms, gs.getReferencePath(), flag_dyn, flag_amb);
                                        else
                                            rtplot_matlab (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath(), flag_amb);
                                        end
                                    else
                                        if (flag_ge == 1), rtplot_googleearth_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), date_R(t,:)), end
                                        if (flag_stopGOstop == 1)
                                            rtplot_matlab_cov_stopGOstop (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), [P1_UTM_E, P1_UTM_N], [P2_UTM_E, P2_UTM_N], flag_ms, gs.getReferencePath(), flag_dyn, flag_amb);
                                        else
                                            rtplot_matlab_cov (t, [Xhat_t_t(1); Xhat_t_t(o1+1); Xhat_t_t(o2+1)], pos_M(:,t), Cee([1 o1+1 o2+1],[1 o1+1 o2+1]), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath(), flag_amb);
                                        end
                                    end
                                    if (flag_amb == 1)
                                        rtplot_amb (t, window, Xhat_t_t(o3+1:o3+nSatTot), sqrt(diag(Cee(o3+1:o3+nSatTot,o3+1:o3+nSatTot))), conf_cs);
                                        pause(0.1);
                                    else
                                        if (flag_skyplot == 1)
                                            rtplot_skyplot (t, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                            rtplot_snr (snr_R(:,t), Eph_t, SP3);
                                        else
                                            rttext_sat (t, azR, elR, snr_R(:,t), conf_sat, pivot, Eph_t, SP3);
                                        end
                                        pause(0.01);
                                    end
                                end
                                w_bar.goTime(t);
                            end
                            
                            w_bar.close();
                            
                            if (flag_stopGOstop == 1)
                                %azimuth computation
                                if (angleDIR < 0)
                                    angleDIR = angleDIR + 360;
                                end
                                
                                %conversion to sexagesimal degrees
                                angleDIR_deg = floor(angleDIR);
                                min_dec = (angleDIR-angleDIR_deg)*60;
                                angleDIR_min = floor(min_dec);
                                angleDIR_sec = (min_dec - angleDIR_min)*60;
                                %---------------------------------
                                sigma_angleDIR_deg = floor(sigma_angleDIR);
                                min_dec = (sigma_angleDIR-sigma_angleDIR_deg)*60;
                                sigma_angleDIR_min = floor(min_dec);
                                sigma_angleDIR_sec = (min_dec - sigma_angleDIR_min)*60;
                                
                                fprintf('\n')
                                fprintf('Estimated azimuth = %d deg %d min %6.3f sec\n', angleDIR_deg, angleDIR_min, angleDIR_sec);
                                fprintf('Standard deviation  = %d deg %d min %6.3f sec\n', sigma_angleDIR_deg, sigma_angleDIR_min, sigma_angleDIR_sec);
                                fprintf('\n')
                            end
                            
                            fclose(fid_dyn);
                            fclose(fid_kal);
                            fclose(fid_sat);
                            fclose(fid_dop);
                            fclose(fid_conf);
                            fclose(fid_res);
                        end
                    elseif (mode_vinc == 1)
                        %--------------------------------------------------------------------------------------------------------------------
                        %% POST-PROCESSING (RELATIVE POSITIONING): KALMAN FILTER ON CODE AND PHASE DOUBLE DIFFERENCES WITH LINE CONSTRAINT
                        %--------------------------------------------------------------------------------------------------------------------
                        
                        fid_kal = fopen([filerootOUT '_kal_000.bin'],'w+');
                        fid_sat = fopen([filerootOUT '_sat_000.bin'],'w+');
                        fid_dop = fopen([filerootOUT '_dop_000.bin'],'w+');
                        fid_conf = fopen([filerootOUT '_conf_000.bin'],'w+');
                        fid_res = fopen([filerootOUT '_res_000.bin'],'w+');
                        
                        %repeat more than once the reference loop
                        %(this constrained mode works only for circuits)
                        ref = gs.getReferencePath();
                        ref_loop = [ref.path; ref.path];
                        
                        if ~state.isSeamlessKF()
                            kalman_initialized = 0;
                        end
                        if (~kalman_initialized)
                            while (~kalman_initialized)
                                if (isempty(time_GPS_diff))
                                    fprintf('It was not possible to initialize the Kalman filter.\n');
                                    return
                                end
                                
                                Eph_t = rt_find_eph (Eph, time_GPS_diff(1), nSatTot);
                                
                                kalman_initialized = goGPS_KF_DD_code_phase_init_vinc(pos_R, pos_M(:,1), time_GPS_diff(1), pr1_R(:,1), pr1_M(:,1), ph1_R(:,1), ph1_M(:,1), dop1_R(:,1), dop1_M(:,1), pr2_R(:,1), pr2_M(:,1), ph2_R(:,1), ph2_M(:,1), dop2_R(:,1), dop2_M(:,1), snr_R(:,1), snr_M(:,1), Eph_t, SP3, iono, lambda, frequencies, state.getProcessingRate(), ref_loop, dtMdot(1));
                                
                                if (~kalman_initialized)
                                    pos_M(:,1) = []; time_GPS_diff(1) = []; time_GPS(1) = []; week_R(1) = [];
                                    pr1_R(:,1) = []; pr1_M(:,1) = []; ph1_R(:,1) = []; ph1_M(:,1) = []; dop1_R(:,1) = []; dop1_M(:,1) = [];
                                    pr2_R(:,1) = []; pr2_M(:,1) = []; ph2_R(:,1) = []; ph2_M(:,1) = []; dop2_R(:,1) = []; dop2_M(:,1) = [];
                                    snr_R(:,1) = []; snr_M(:,1) = []; dtMdot(1) = [];
                                end
                            end
                            
                            fwrite(fid_kal, [Xhat_t_t; Yhat_t_t; Cee(:)], 'double');
                            fwrite(fid_sat, nSatTot, 'int8');
                            fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                            fwrite(fid_dop, [PDOP; HDOP; VDOP; 0; 0; 0], 'double');
                            fwrite(fid_conf, [nSatTot; conf_sat; conf_cs; pivot], 'int8');
                            fwrite(fid_res, nSatTot, 'int8');
                            fwrite(fid_res, [residuals_fixed(1:nSatTot*2); residuals_fixed(nSatTot*2+1:end);residuals_float(1:nSatTot*2); residuals_float(nSatTot*2+1:end);outliers(1:nSatTot*2);outliers(nSatTot*2+1:end)], 'double');
                            
                            if (flag_plotproc)
                                if (flag_ge == 1), rtplot_googleearth (1, [Yhat_t_t(1); Yhat_t_t(2); Yhat_t_t(3)], pos_M(:,1), date_R(1,:)), end
                                rtplot_matlab (1, [Yhat_t_t(1); Yhat_t_t(2); Yhat_t_t(3)], pos_M(:,1), 0, 0, 0, 0, flag_ms, gs.getReferencePath(), flag_amb);
                                if (flag_amb == 1)
                                    rtplot_amb (1, window, Xhat_t_t(o1+1:o1+nSatTot), sqrt(diag(Cee(o1+1:o1+nSatTot,o1+1:o1+nSatTot))), conf_cs);
                                else
                                    if (flag_skyplot == 1)
                                        rtplot_skyplot (1, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                        rtplot_snr (snr_R(:,1), Eph_t, SP3);
                                    else
                                        rttext_sat (1, azR, elR, snr_R(:,1), conf_sat, pivot, Eph_t, SP3);
                                    end
                                end
                            end
                            t1 = 2;
                        else
                            fwrite(fid_sat, nSatTot, 'int8');
                            fwrite(fid_conf, nSatTot, 'int8');
                            fwrite(fid_res, nSatTot, 'int8');
                            t1 = 1;
                        end
                        
                        % goGPS waiting bar
                        w_bar.setBarLen(length(time_GPS_diff));
                        w_bar.createNewBar('Processing...');
                        
                        for t = t1 : length(time_GPS_diff)
                            residuals_fixed=NaN(4*nSatTot,1);
                            residuals_float=NaN(4*nSatTot,1);
                            outliers=zeros(4*nSatTot,1);
                            
                            Eph_t = rt_find_eph (Eph, time_GPS_diff(t), nSatTot);
                            
                            [check_on, check_off, check_pivot, check_cs] = goGPS_KF_DD_code_phase_loop_vinc(pos_M(:,t), time_GPS_diff(t), pr1_R(:,t), pr1_M(:,t), ph1_R(:,t), ph1_M(:,t), dop1_R(:,t), dop1_M(:,t), pr2_R(:,t), pr2_M(:,t), ph2_R(:,t), ph2_M(:,t), dop2_R(:,t), dop2_M(:,t), snr_R(:,t), snr_M(:,t), Eph_t, SP3, iono, lambda, frequencies, state.getProcessingRate(), ref_loop, dtMdot(t));
                            
                            fwrite(fid_kal, [Xhat_t_t; Yhat_t_t; Cee(:)], 'double');
                            fwrite(fid_sat, [azM; azR; elM; elR; distM; distR], 'double');
                            fwrite(fid_dop, [PDOP; HDOP; VDOP; 0; 0; 0], 'double');
                            fwrite(fid_conf, [conf_sat; conf_cs; pivot], 'int8');
                            fwrite(fid_res, [residuals_fixed(1:nSatTot*2); residuals_fixed(nSatTot*2+1:end);residuals_float(1:nSatTot*2); residuals_float(nSatTot*2+1:end);outliers(1:nSatTot*2);outliers(nSatTot*2+1:end)], 'double');
                            
                            if (flag_plotproc)
                                if (mode_user == 1 && t == 2), w_bar.shiftDown(); end
                                if (flag_ge == 1), rtplot_googleearth (t, [Yhat_t_t(1); Yhat_t_t(2); Yhat_t_t(3)], pos_M(:,t), date_R(t,:)), end
                                rtplot_matlab (t, [Yhat_t_t(1); Yhat_t_t(2); Yhat_t_t(3)], pos_M(:,t), check_on, check_off, check_pivot, check_cs, flag_ms, gs.getReferencePath(), flag_amb);
                                if (flag_amb == 1)
                                    rtplot_amb (t, window, Xhat_t_t(o1+1:o1+nSatTot), sqrt(diag(Cee(o1+1:o1+nSatTot,o1+1:o1+nSatTot))), conf_cs);
                                    pause(0.1);
                                else
                                    if (flag_skyplot == 1)
                                        rtplot_skyplot (t, azR, elR, conf_sat, pivot, Eph_t, SP3);
                                        rtplot_snr (snr_R(:,t), Eph_t, SP3);
                                    else
                                        rttext_sat (t, azR, elR, snr_R(:,t), conf_sat, pivot, Eph_t, SP3);
                                    end
                                    pause(0.01);
                                end
                            end
                            
                            w_bar.goTime(t);
                        end
                        
                        w_bar.close();
                        
                        fclose(fid_kal);
                        fclose(fid_sat);
                        fclose(fid_dop);
                        fclose(fid_conf);
                        fclose(fid_res);
                    end
                case goGNSS.MODE_RT_R_MON
                    %----------------------------------------------------------------------------------------------
                    %% REAL-TIME: ROVER MONITORING
                    %----------------------------------------------------------------------------------------------
                    
                    goGPS_rover_monitor(filerootOUT, protocol_idx, flag_var_dyn_model, flag_stopGOstop, state.getCaptureRate(), constellations);
                    
                case goGNSS.MODE_RT_M_MON
                    %----------------------------------------------------------------------------------------------
                    %% REAL-TIME: MASTER MONITORING
                    %----------------------------------------------------------------------------------------------
                    
                    goGPS_master_monitor(filerootOUT, flag_NTRIP);
                    
                case goGNSS.MODE_RT_RM_MON
                    %----------------------------------------------------------------------------------------------
                    %% REAL-TIME: ROVER AND MASTER MONITORING
                    %----------------------------------------------------------------------------------------------
                    
                    goGPS_realtime_monitor(filerootOUT, protocol_idx, flag_NTRIP, flag_ms_pos, flag_var_dyn_model, flag_stopGOstop, pos_M, constellations);
                    
                case goGNSS.MODE_RT_NAV
                    %----------------------------------------------------------------------------------------------
                    %% REAL-TIME: KALMAN FILTER ON PHASE AND CODE DOUBLE DIFFERENCES WITH/WITHOUT A CONSTRAINT
                    %----------------------------------------------------------------------------------------------
                    
                    goGPS_realtime(filerootOUT, protocol_idx, mode_vinc, flag_ms, flag_ge, flag_cov, flag_NTRIP, flag_ms_pos, flag_skyplot, flag_plotproc, flag_var_dyn_model, flag_stopGOstop, gs.getReferencePath(), pos_M, dop1_M, pr2_M, pr2_R, ph2_M, ph2_R, dop2_M, dop2_R, constellations);
            end
            
            if (state.isModeLS()) %remove unused epochs from time_GPS_diff (for LS modes)
                if (any(unused_epochs))
                    time_GPS_diff(logical(unused_epochs)) = [];
                    time_GPS(logical(unused_epochs)) = [];
                    week_R(logical(unused_epochs)) = [];
                end
            end
        end
        
        if (~state.isModeSEID())
            %----------------------------------------------------------------------------------------------
            %% OUTPUT DATA FILE READING
            %----------------------------------------------------------------------------------------------
            
            %if any positioning was done (either post-processing or real-time)
            if goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)
                
                %time adjustments (to account for sub-integer approximations in MATLAB - thanks to radiolabs.it for pointing this out!)
                time_GPS = time_GPS_diff + zero_time;
                time_R   = time_R_diff   + zero_time;
                if goGNSS.isDD(mode)
                    time_M   = time_M_diff + zero_time;
                end
                % Eph(32,:) = Eph(32,:) + zero_time; % not needed
                % Eph(33,:) = Eph(33,:) + zero_time; % not needed
                if (flag_SP3)
                    SP3.time    = SP3.time + zero_time;
                    SP3.time_hr = SP3.time_hr + zero_time;
                    SP3.t_sun   = SP3.t_sun + zero_time;
                    if ~isempty(SP3.ERP)
                        SP3.ERP.t   = SP3.ERP.t - zero_time;
                    end
                end
                
                %---------------------------------
                
                %reading of the files with Kalman filter results
                [Xhat_t_t_OUT, Yhat_t_t_OUT, Cee_OUT, azM, azR, elM, elR, distM, distR, ...
                    conf_sat_OUT, conf_cs, pivot_OUT, PDOP, HDOP, VDOP, KPDOP, KHDOP, KVDOP, ...
                    RES_CODE1_FIXED, RES_CODE2_FIXED, RES_PHASE1_FIXED, RES_PHASE2_FIXED,...
                    RES_CODE1_FLOAT, RES_CODE2_FLOAT, RES_PHASE1_FLOAT, RES_PHASE2_FLOAT,...
                    outliers_CODE1, outliers_CODE2, outliers_PHASE1, outliers_PHASE2, apriori_ZHD_OUT, STDs] = load_goGPSoutput(filerootOUT, mode, mode_vinc);
                
                if (state.isModeBlock())
                    [pos_KAL, Xhat_t_t_OUT, conf_sat_OUT, Cee_OUT, pivot_OUT, nsat, fixed_amb] = go_block.getLegacyOutput();
                end
                
                %variable saving for final graphical representations
                if (~state.isModeBlock())
                    nSol = size(Xhat_t_t_OUT,2);
                    pos_KAL = zeros(3,nSol);
                    pos_REF = zeros(3,nSol);
                    stat_SC = zeros(2,nSol);
                    if (any(fixed_solution))
                        fixed_amb = fixed_solution;
                    else
                        fixed_amb = zeros(1,nSol);
                        succ_rate = zeros(1,nSol);
                    end
                    estim_amb = zeros(nSatTot,nSol);
                    sigma_amb = zeros(nSatTot,nSol);
                    estim_tropo = zeros(nSatTot,nSol);
                    estim_gradN = zeros(nSatTot,nSol);
                    estim_gradE = zeros(nSatTot,nSol);
                    for i = 1 : nSol
                        if (mode == goGNSS.MODE_PP_KF_CP_DD && mode_vinc == 1)
                            pos_KAL(:,i) = [Yhat_t_t_OUT(1,i); Yhat_t_t_OUT(2,i); Yhat_t_t_OUT(3,i)];
                            estim_amb(:,i) = Xhat_t_t_OUT(o1+1:o1+nSatTot,i);
                            sigma_amb(:,i) = sqrt(diag(Cee_OUT(o1+1:o1+nSatTot,o1+1:o1+nSatTot,i)));
                        else
                            pos_KAL(:,i) = [Xhat_t_t_OUT(1,i); Xhat_t_t_OUT(o1+1,i); Xhat_t_t_OUT(o2+1,i)];
                            estim_amb(:,i) = Xhat_t_t_OUT(o3+1:o3+nSatTot,i);
                            sigma_amb(:,i) = sqrt(diag(Cee_OUT(o3+1:o3+nSatTot,o3+1:o3+nSatTot,i)));
                        end
                        
                        if (exist('antoff_R','var'))
                            pos_KAL(:,i) = local2globalPos(-antoff_R(:,1,1), pos_KAL(:,i));
                        end
                        %if tropospheric delay was estimated in PPP
                        if (state.isSA(mode) && flag_tropo)
                            if (~flag_tropo_gradient)
                                estim_tropo = Xhat_t_t_OUT(end-nC-2,:);
                            else
                                estim_tropo = Xhat_t_t_OUT(end-nC-2,:);
                                estim_gradN = Xhat_t_t_OUT(end-nC-1,:);
                                estim_gradE = Xhat_t_t_OUT(end-nC  ,:);
                            end
                        else
                            estim_tropo = zeros(nSol, 1);
                            estim_gradN = zeros(nSol, 1);
                            estim_gradE = zeros(nSol, 1);
                        end
                    end
                else
                    nSol = 1;
                    estim_tropo = zeros(nSol, 1);
                    succ_rate = zeros(1,nSol);
                    
                    if (exist('antoff_R','var'))
                        pos_KAL(:,1) = local2globalPos(-antoff_R(:,1,1), Xhat_t_t_OUT(:,1));
                    end
                end
                
                
                switch state.getForwardBackwardKF()
                    case -1
                        for i = 1 : nSol/2
                            %if relative positioning (i.e. with master station)
                            if goGNSS.isDD(mode) && goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)
                                pos_REF(:,nSol/2-i+1) = [pos_M(1,i); pos_M(2,i); pos_M(3,i)];
                                pos_REF(:,nSol/2+i) = [pos_M(1,i); pos_M(2,i); pos_M(3,i)];
                            else
                                pos_REF(:,i) = pos_KAL(:,1);
                                pos_REF(:,nSol-i+1) = pos_KAL(:,1);
                            end
                        end
                    case 1
                        for i = 1 : nSol/2
                            %if relative positioning (i.e. with master station)
                            if goGNSS.isDD(mode) && goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)
                                pos_REF(:,i) = [pos_M(1,i); pos_M(2,i); pos_M(3,i)];
                                pos_REF(:,nSol-i+1) = [pos_M(1,i); pos_M(2,i); pos_M(3,i)];
                            else
                                pos_REF(:,i) = pos_KAL(:,1);
                                pos_REF(:,nSol-i+1) = pos_KAL(:,1);
                            end
                        end
                    case 0
                        for i = 1 : nSol
                            %if relative positioning (i.e. with master station)
                            if goGNSS.isDD(mode) && goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)
                                pos_REF(:,i) = [pos_M(1,i); pos_M(2,i); pos_M(3,i)];
                            else
                                pos_REF(:,i) = pos_KAL(:,1);
                            end
                        end
                end
                
                %----------------------------------------------------------------------------------------------
                %% OUTPUT FILE SAVING (TEXT FILE)
                %----------------------------------------------------------------------------------------------
                
                nsat = sum(abs(conf_sat_OUT),1);
                
                %if any positioning was done (either post-processing or real-time)
                if goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)
                    %display information
                    fprintf('Writing output file...\n');
                    %cartesian coordinates (X,Y,Z)
                    X_KAL = pos_KAL(1,:)';
                    Y_KAL = pos_KAL(2,:)';
                    Z_KAL = pos_KAL(3,:)';
                    
                    if (any(pos_KAL))
                        %coordinate transformation (geodetic)
                        [phi_KAL, lam_KAL, h_KAL] = cart2geod(X_KAL, Y_KAL, Z_KAL);
                        phi_KAL = phi_KAL * 180/pi;
                        lam_KAL = lam_KAL * 180/pi;
                        
                        %coordinate transformation (UTM)
                        [EAST_UTM, NORTH_UTM, h_UTM, utm_zone] = cart2plan(X_KAL, Y_KAL, Z_KAL);
                    else
                        phi_KAL = nan(size(pos_KAL,2));
                        lam_KAL = nan(size(pos_KAL,2));
                        h_KAL = nan(size(pos_KAL,2));
                        EAST_UTM = nan(size(pos_KAL,2));
                        NORTH_UTM = nan(size(pos_KAL,2));
                        h_UTM = nan(size(pos_KAL,2));
                        utm_zone = nan(size(pos_KAL,2));
                    end
                    
                    %if relative positioning (i.e. with master station)
                    if goGNSS.isDD(mode) || (mode == goGNSS.MODE_RT_NAV)
                        X_ENU = global2localPos(pos_KAL, pos_REF);
                        EAST_KAL  = X_ENU(1,:)';
                        NORTH_KAL = X_ENU(2,:)';
                        UP_KAL    = X_ENU(3,:)';
                        
                        EAST  = EAST_KAL;
                        NORTH = NORTH_KAL;
                    else
                        EAST_KAL  = zeros(size(EAST_UTM));
                        NORTH_KAL = zeros(size(NORTH_UTM));
                        UP_KAL    = zeros(size(h_UTM));
                        
                        EAST  = EAST_UTM;
                        NORTH = NORTH_UTM;
                    end
                    
                    %if no Kalman filter is used or if the positioning is constrained
                    if (mode_vinc == 1) || (mode == goGNSS.MODE_PP_LS_C_SA) || (mode == goGNSS.MODE_PP_LS_CP_SA) || (mode == goGNSS.MODE_PP_LS_C_DD) || (mode == goGNSS.MODE_PP_LS_CP_DD_L)
                        %initialization to -9999 (no data available)
                        KHDOP(1:nSol) = -9999;
                    end
                    N = [];
                    h_ortho(1:nSol) = -9999;
                    
                     %******TIME FORMATTING
                     %Convert GPS Time & GPS Week to Time of Week(tow)
                    [tow] = weektime2tow(week_R(:,1,1), time_GPS);
                    
                    %******DATE FORMATTING
                    %Convert GPS Time to UTC Time & Day of Year(DOY)
                    [date_R, doy] = gps2date(week_R(:,1,1), tow);%UTC TIME IN [YEAR MM DD HR MN SS] 

                    %CONVERT GREGORIAN /UTC DATE TO MODIFIED JULIAN DATE(MJD) FOR LATER USE
                    MJD = date2mjd(date_R);

                    for i = 1 : nSol
                        if (geoid.ncols ~= 0)
                            %geoid undulation interpolation
                            N = grid_bilin_interp(lam_KAL(i), phi_KAL(i), geoid.grid, geoid.ncols, geoid.nrows, geoid.cellsize, geoid.Xll, geoid.Yll, -9999);
                            %orthometric height
                            h_ortho(i) = h_KAL(i) - N;
                        end
                    end
                    
                    
                    %---------------------------------------------------------------------+
                    %% ****COMPUTE ZENITH TROPOSPHERIC DELAYs FROM TROPO MODELS FOR OUTPUT
                    %---------------------------------------------------------------------+
                   
                    %Call the "tropo_error_correction.m" function
                    [~,~,~,ZTD,ZHD,ZWD] = tropo_error_correction(date_R,[phi_KAL,lam_KAL, h_KAL],[] ) ;
                    
                    
                    %IF TROPO DELAY IS ESTIMATED USING PPP 
                    if state.isTropoOn()
                       
                       %CHECK FOR NAN
                       na = isnan(estim_tropo(:));

                       estim_tropo(na)=0;%REPLACE NAN ENTRIES WITH ZEROS(0)
                    
                       %TRANSPOSE ROW MATRIX(1xn) TO COLUMN MATRIX(nx1)
                       if size(estim_tropo) ~= size(phi_KAL)  
                          estim_tropo = estim_tropo';
                       end  
                        
                       ZTDppp = estim_tropo(:);
                       
                       if any([~isempty(estim_tropo),all(estim_tropo(:) ~=0),any(find(estim_tropo(:) ~=0)),any(estim_tropo(:) ~=0)])
                       
                          ZWDppp = ZTDppp - ZHD; 
                       
                       end
                       
                    end
                        
                    %IF USER WANTS TO EXTRACT VMF GRIDDED ZENITH DELAYS IN ADDITION
                    if getappdata(0,'getvmfZTDs_state')== 1
        
                       %*********GET STRUCT VMF GRID FILES
                       VMFgrids = getappdata(0,'VMFgrids');

                       [ZHDvmf,ZWDvmf,ZTDvmf] = SearchReadVMFgrids(date_R,[phi_KAL,lam_KAL, h_KAL],VMFgrids);
      
                       if any([isempty(ZHDvmf),isempty(ZWDvmf),isempty(ZTDvmf)])
            
                          ZHDvmf = zeros(size(date_R,1),1);
           
                          [ZWDvmf,ZTDvmf] = deal(ZHDvmf); 

                       end   
                       
                       
                    end  

                    %-----------------------------------------------------------------------------------------------
                    %% **************************PWV RETRIEVAL
                    %-----------------------------------------------------------------------------------------------
                    %CHECK IF USER WANTS TO COMPUTE/RETRIEVE PRECIPITABLE WATER VAPOUR(PWV)
                    if getappdata(0,'pwv_state') == 1 %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
          
                       if state.isTropoOn() 
                           
                          if any([~isempty(estim_tropo),all(estim_tropo(:) ~=0),any(find(estim_tropo(:) ~=0)),any(estim_tropo(:) ~=0)])
                           
                             %Call "the compute_PWV.m" Sub-routine
                             [PWV,Tm,PI] = compute_PWV(date_R,[phi_KAL,lam_KAL, h_KAL],ZTDppp,ZHD);
                       
                          else  
                              %Call "the compute_PWV.m" Sub-routine
                              [PWV,Tm,PI] = compute_PWV(date_R,[phi_KAL,lam_KAL, h_KAL],ZTD,ZHD);
       
                          end 
                          
                       else
                           %Call "the compute_PWV.m" Sub-routine
                           [PWV,Tm,PI] = compute_PWV(date_R,[phi_KAL,lam_KAL, h_KAL],ZTD,ZHD);
       
                       end
                       
                    end %//if getappdata(0,'pwv_state') == 1
                     
                    %           -----------------------------------------------
                    %*********GET MODELLING OPTIONS FROM goGPS GUI FOR OUTPUT
                    %           -----------------------------------------------
                    %(1.0)GET OPTION FOR TROPO MODELLING
                    option_Cmodels = getappdata(0,'option_Cmodels');%option button for combine tropo models
                    option_Smodels = getappdata(0,'option_Smodels');%option button for separate tropo models
     
                    %(2.0)*****************GET TROPO MODELs
                    %(2.1)COMBINE TROPO MODEL
                    TROPO_Cmodel = getappdata(0,'TROPO_Cmodel');
     
                    %(2.2)SEPARATE TROPO MODELs
                    dryModel = getappdata(0,'dryModel');%get selected Hydrostatic model(string)from pop_ITropo_hydrostatic
                    wetModel = getappdata(0,'wetModel');%get selected Wet model(string)from pop_ITropo_wet
   
                    DryModel = getappdata(0,'DryModel');%Hydrostatic model for goGPS output file
                    WetModel = getappdata(0,'WetModel');%Wet model for goGPS output file
     
                    %(2.2.1)GET OPTIONs FOR USE OF MODELs / OTHER MAPPING FUNCTION
                    option_model_MF     = getappdata(0,'option_model_MF');%option button for models MF
                    option_different_MF = getappdata(0,'option_different_MF');%option button for different MF models

                    %(2.2.2)GET SELECTED MAPPING FUNCTION MODELs
                    MFh_Model = getappdata(0,'MFh_Model');%Hydrostatic mapping function
                    MFw_Model = getappdata(0,'MFw_Model');%Wet mapping function
     
                    MF_Model  = getappdata(0,'MF_Model');%Combine models mapping function
     
                    %(3.0)GET SOURCE METEOROLOGICAL PARAMETERS
                    metPARA_typeC   = getappdata(0,'metPARA_typeC');%MET STRING FOR COMBINE TROPO MODELS
                    metPARA_typeD   = getappdata(0,'metPARA_typeD');%MET STRING FOR DRY TROPO MODELS
                    metPARA_typeW   = getappdata(0,'metPARA_typeW');%MET STRING FOR WET TROPO MODELS
                    metPARA_typeW_D = getappdata(0,'metPARA_typeW_D');%MET STRING FOR WET TROPO MODELS
     
                   %GET SELECTED WEIGHTED MEAN TEMPERATURE(Tm) MODEL TYPE
                   TmModel_type = getappdata(0,'TmModel_type');
     
                   %======================WRITE SELECTED OPTIONS ========================
                   %                      -----------------------
                   if option_Cmodels == 1 & option_Smodels == 0 %if user selects option for combine tropo model
        
                      tropomodel = sprintf('  Tropospheric Delay Correction model     : %s \n',TROPO_Cmodel); 
                      MF_model   = sprintf('  Mapping Function model                  : %s \n',MF_Model);
                      METsource  = sprintf('  Source of Meteorological Parameters     : %s \n',metPARA_typeC);
          
                      %*****PWV
                      if getappdata(0,'pwv_state') == 1 %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                         TmModel = sprintf('  Weighted Mean Temperature(Tm) Model     : %s \n\n',TmModel_type);
                      else 
                          TmModel = sprintf('   \n');
                      end 
          
                   elseif option_Smodels == 1 & option_Cmodels == 0 %if user selects option for separate tropo model
            
                          DRYmodel = sprintf('  Zenith Hydrostatic Delay(ZHD)Model      : %s \n',DryModel);
                          WETmodel = sprintf('  Zenith Wet Delay(ZWD) Model             : %s \n',WetModel);
             
                          if option_model_MF == 1 & option_different_MF == 0 %if user selects option for use of tropo models  MF 
               
                             hMF_model = sprintf('  Hydrostatic Mapping Function model      : %s \n',MFh_Model); 
                             wMF_model = sprintf('  Wet Mapping Function model              : %s \n',MFw_Model); 
                  
                          elseif option_different_MF == 1 & option_model_MF == 0 %if user selects option for use of different MF model
                  
                                 hMF_model = sprintf('  Hydrostatic Mapping Function model      : %s \n',MFh_Model); 
                                 wMF_model = sprintf('  Wet Mapping Function model              : %s \n',MFw_Model); 
                  
                          end   
          
                          if any([any([strcmpi(dryModel,wetModel),strfind(dryModel,wetModel)]),metPARA_typeW_D])
                   
                             METsourceWD=sprintf('  Source of Meteorological Parameters     : %s \n',metPARA_typeD); 
           
                          else   
                              METsourceD=sprintf('  Source of Meteorological Parameters[ZHD]: %s \n',metPARA_typeD);
                              METsourceW=sprintf('  Source of Meteorological Parameters[ZWD]: %s \n',metPARA_typeW);
                          end  
           
                          %*****PWV
                          if getappdata(0,'pwv_state') == 1 %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                             TmModel = sprintf('  Weighted Mean Temperature(Tm) Model     : %s \n\n',TmModel_type);
                          else   
                              TmModel = sprintf('   \n');
                          end   
           
                   end    
     
                    %-------------------------------------------------------------------------------------------------------------------+
                    %-------------------------------------------------------------------------------------------------------------------+

                    %*******FILE SAVING AS TEXT
                    if (strcmp(fsep_char,'default'))
                        
                        %WHEN VMF ZENITH DELAYS, PWV & PPP DELAYS NOT ICLUDED IN THE OUTPUT
                       if all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
              
                          head_str='    Date        GPS time        GPS week       Doy        Latitude       Longitude     h (ellips.)     ECEF X           ECEF Y          ECEF Z          UTM North        UTM East       h (orthom.)     UTM zone        Num. Sat.     HDOP        KHDOP     Local North      Local East          Local H           Ambiguity fix       Success rate       ZHD            ZWD            ZTD \n';                                                                                                                                                 
                          head_str_unit='[yyyy/mm/dd]  [hr:mn:ss.ss]                                [deg]            [deg]         [m]            [m]              [m]             [m]                [m]              [m]            [m]                                                                                                                                                          [m]            [m]            [m]    \n'; 
                          head_str_='--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n';
                          row_str = '%04d/%02d/%02d    %02d:%02d:%06.2f %10d %13.2f %16.8f %16.8f %10.4f %16.4f %16.4f %16.4f %16.4f %16.4f %12.4f %11s%16d  %11.3f %11.3f    %12.4f%17.4f%17.4f %16d %21.4f %14.4f %14.4f %14.4f\n';  

                          %WHEN ONLY VMF ZENITH DELAYS ARE INCLUDED
                       elseif all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==0])

                              head_str='    Date        GPS time        GPS week       Doy        Latitude       Longitude     h (ellips.)     ECEF X           ECEF Y          ECEF Z          UTM North        UTM East       h (orthom.)     UTM zone        Num. Sat.     HDOP        KHDOP     Local North      Local East          Local H           Ambiguity fix       Success rate       ZHD            ZWD            ZTD            ZHDvmf            ZWDvmf            ZTDvmf \n';                                                                                                                                                 
                              head_str_unit='[yyyy/mm/dd]  [hr:mn:ss.ss]                                [deg]            [deg]         [m]            [m]              [m]             [m]                [m]              [m]            [m]                                                                                                                                                          [m]            [m]            [m]             [m]               [m]               [m]     \n'; 
                              head_str_='---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n';
                              row_str = '%04d/%02d/%02d    %02d:%02d:%06.2f %10d %13.2f %16.8f %16.8f %10.4f %16.4f %16.4f %16.4f %16.4f %16.4f %12.4f %11s%16d  %11.3f %11.3f    %12.4f%17.4f%17.4f %16d %21.4f %14.4f %14.4f %14.4f  %14.4f    %14.4f    %14.4f\n';  

                             %WHEN ONLY PWVs ARE INCLUDED
                       elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0])
                
                              head_str='    Date        GPS time        GPS week       Doy        Latitude       Longitude     h (ellips.)     ECEF X           ECEF Y          ECEF Z          UTM North        UTM East       h (orthom.)     UTM zone        Num. Sat.     HDOP        KHDOP     Local North      Local East          Local H           Ambiguity fix       Success rate       ZHD            ZWD            ZTD            Tm            PI            PWV \n';                                                                                                                                                 
                              head_str_unit='[yyyy/mm/dd]  [hr:mn:ss.ss]                                [deg]            [deg]         [m]            [m]              [m]             [m]                [m]              [m]            [m]                                                                                                                                                          [m]            [m]            [m]            [K]                         [mm]    \n'; 
                              head_str_='---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n';
                              row_str = '%04d/%02d/%02d    %02d:%02d:%06.2f %10d %13.2f %16.8f %16.8f %10.4f %16.4f %16.4f %16.4f %16.4f %16.4f %12.4f %11s%16d  %11.3f %11.3f    %12.4f%17.4f%17.4f %16d %21.4f %14.4f %14.4f %14.4f  %14.3f%14.5f%14.4f\n';  

                          %WHEN PWV & VMF ZENITH DELAYS ARE INCLUDED
                       elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==0]) %WHEN BOTH PWV & VMF ZTDs ARE INCLUDED TO THE OUTPUT
                              
                              head_str='    Date        GPS time        GPS week       Doy        Latitude       Longitude     h (ellips.)     ECEF X           ECEF Y          ECEF Z          UTM North        UTM East       h (orthom.)     UTM zone        Num. Sat.     HDOP        KHDOP     Local North      Local East          Local H           Ambiguity fix       Success rate       ZHD            ZWD            ZTD            ZHDvmf            ZWDvmf            ZTDvmf            Tm            PI           PWV  \n';                                                                                                                                                 
                              head_str_unit='[yyyy/mm/dd]  [hr:mn:ss.ss]                                [deg]            [deg]         [m]            [m]              [m]             [m]                [m]              [m]            [m]                                                                                                                                                          [m]            [m]            [m]             [m]               [m]               [m]             [K]                         [mm]         \n'; 
                              head_str_='------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ \n';
                              row_str = '%04d/%02d/%02d    %02d:%02d:%06.2f %10d %13.2f %16.8f %16.8f %10.4f %16.4f %16.4f %16.4f %16.4f %16.4f %12.4f %11s%16d  %11.3f %11.3f    %12.4f%17.4f%17.4f %16d %21.4f %14.4f %14.4f %14.4f  %14.4f    %14.4f    %14.4f   %14.4f%14.5f %14.5f\n';  

                         %WHEN VMF ZENITH DELAYS & PWV ARE NOT ICLUDED BUT PPP DELAYS ARE INCLUDED IN THE OUTPUT
                       elseif all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==1 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                    
                              head_str='    Date        GPS time        GPS week       Doy        Latitude       Longitude     h (ellips.)     ECEF X           ECEF Y          ECEF Z          UTM North        UTM East       h (orthom.)     UTM zone        Num. Sat.     HDOP        KHDOP     Local North      Local East          Local H           Ambiguity fix       Success rate      ZHD            ZWD            ZTD           ZWDppp         ZTDppp \n';                                                                                                                                                 
                              head_str_unit='[yyyy/mm/dd]  [hr:mn:ss.ss]                                [deg]            [deg]         [m]            [m]              [m]             [m]                [m]              [m]            [m]                                                                                                                                                         [m]            [m]            [m]             [m]            [m]    \n'; 
                              head_str_='--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n';
                              row_str = '%04d/%02d/%02d    %02d:%02d:%06.2f %10d %13.2f %16.8f %16.8f %10.4f %16.4f %16.4f %16.4f %16.4f %16.4f %12.4f %11s%16d  %11.3f %11.3f %12.4f %17.4f %17.4f %16d %21.4f %14.4f %14.4f %14.4f %14.4f %14.4f\n';  

                           %WHEN PWV IS NOT ICLUDED BUT PPP & VMF ZENITH DELAYS &  ARE INCLUDED IN THE OUTPUT
                       elseif  all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==1])
                           
                               head_str='    Date        GPS time        GPS week       Doy        Latitude       Longitude     h (ellips.)     ECEF X           ECEF Y          ECEF Z          UTM North        UTM East       h (orthom.)     UTM zone        Num. Sat.      HDOP       KHDOP     Local North      Local East         Local H          Ambiguity fix       Success rate       ZHD            ZWD            ZTD            ZWDppp         ZTDppp          ZHDvmf            ZWDvmf            ZTDvmf \n';                                                                                                                                                 
                               head_str_unit='[yyyy/mm/dd]  [hr:mn:ss.ss]                                [deg]            [deg]         [m]            [m]              [m]             [m]              [m]              [m]            [m]                                                                                                                                                          [m]            [m]            [m]              [m]            [m]            [m]               [m]               [m]    \n'; 
                               head_str_='------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ \n';
                               row_str = '%04d/%02d/%02d    %02d:%02d:%06.2f %10d %13.2f %16.8f %16.8f %10.4f %16.4f %16.4f %16.4f %16.4f %16.4f %12.4f %11s %16d %11.3f %11.3f %12.4f %17.4f %17.4f %16d %21.4f %13.4f %14.4f %14.4f%15.4f %14.4f %15.4f %17.4f %17.4f\n';
          
                           %WHEN PWV & PPP ZENITH DELAYS ARE INCLUDED BUT VMF ZTDs ARE NOT ICLUDED  IN THE OUTPUT      
                       elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==1])
                      
                              head_str='    Date        GPS time        GPS week       Doy        Latitude       Longitude     h (ellips.)     ECEF X           ECEF Y          ECEF Z          UTM North        UTM East       h (orthom.)     UTM zone        Num. Sat.      HDOP       KHDOP     Local North      Local East         Local H          Ambiguity fix       Success rate       ZHD            ZWD            ZTD           ZWDppp         ZTDppp        Tm              PI            PWV \n';                                                                                                                                                 
                              head_str_unit='[yyyy/mm/dd]  [hr:mn:ss.ss]                                [deg]            [deg]         [m]            [m]              [m]             [m]              [m]              [m]            [m]                                                                                                                                                          [m]            [m]            [m]             [m]            [m]         [K]                           [mm]    \n'; 
                              head_str_='------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ \n';
                              row_str = '%04d/%02d/%02d    %02d:%02d:%06.2f %10d %13.2f %16.8f %16.8f %10.4f %16.4f %16.4f %16.4f %16.4f %16.4f %12.4f %11s %16d %11.3f %11.3f %12.4f%17.4f %17.4f %16d %21.4f %14.4f %14.4f %14.4f%14.4f %14.4f %14.3f %14.5f %14.4f\n';
        
                       else %WHEN ALL(PWV,VMF & PPP ZENITH DELAYS) ARE INCLUDED IN THE OUTPUT
                           
                           head_str='    Date        GPS time        GPS week       Doy        Latitude       Longitude     h (ellips.)     ECEF X           ECEF Y          ECEF Z          UTM North        UTM East       h (orthom.)     UTM zone        Num. Sat.      HDOP       KHDOP     Local North      Local East         Local H          Ambiguity fix       Success rate       ZHD             ZWD            ZTD            ZWDppp         ZTDppp         ZHDvmf        ZWDvmf         ZTDvmf         Tm              PI            PWV \n';                                                                                                                                                 
                           head_str_unit='[yyyy/mm/dd]  [hr:mn:ss.ss]                                [deg]            [deg]         [m]            [m]              [m]             [m]              [m]              [m]            [m]                                                                                                                                                          [m]             [m]            [m]              [m]            [m]            [m]           [m]            [m]          [K]                           [mm]    \n'; 
                           head_str_='------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n';
                           row_str = '%04d/%02d/%02d    %02d:%02d:%06.2f %10d %13.2f %16.8f %16.8f %10.4f %16.4f %16.4f %16.4f %16.4f %16.4f %12.4f %11s%16d  %11.3f %11.3f %12.4f %17.4f%17.4f  %16d %21.4f%14.4f  %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.3f %14.5f %14.4f\n';
    
                       end
         
                    else 
                         %*****PWV
                        if getappdata(0,'pwv_state') == 0 %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                           head_str = strcat('Date',fsep_char,'GPS time',fsep_char,'GPS week',fsep_char,'GPS tow',fsep_char,'Latitude',fsep_char,'Longitude',fsep_char,'h (ellips.)',fsep_char,'ECEF X',fsep_char,'ECEF Y',fsep_char,'ECEF Z',fsep_char,'UTM North',fsep_char,'UTM East',fsep_char,'h (orthom.)',fsep_char,'UTM zone',fsep_char,'Num. Sat.',fsep_char,'HDOP',fsep_char,'KHDOP',fsep_char,'Local North',fsep_char,'Local East',fsep_char,'Local H',fsep_char,'Ambiguity fix',fsep_char,'Success rate',fsep_char,'ZHD',fsep_char,'ZWD',fsep_char,'ZTD\n');
                           row_str = strcat('%02d/%02d/%02d',fsep_char,'%02d:%02d:%f',fsep_char,'%d',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%s',fsep_char,'%d',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%d',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f\n'); 
          
                        else %IF PWV IS RETRIEVED 
                            head_str = strcat('Date',fsep_char,'GPS time',fsep_char,'GPS week',fsep_char,'GPS tow',fsep_char,'Latitude',fsep_char,'Longitude',fsep_char,'h (ellips.)',fsep_char,'ECEF X',fsep_char,'ECEF Y',fsep_char,'ECEF Z',fsep_char,'UTM North',fsep_char,'UTM East',fsep_char,'h (orthom.)',fsep_char,'UTM zone',fsep_char,'Num. Sat.',fsep_char,'HDOP',fsep_char,'KHDOP',fsep_char,'Local North',fsep_char,'Local East',fsep_char,'Local H',fsep_char,'Ambiguity fix',fsep_char,'Success rate',fsep_char,'ZHD',fsep_char,'ZWD',fsep_char,'ZTD',fsep_char,'Tm',fsep_char,'PI',fsep_char,'PWV\n');
                            row_str = strcat('%02d/%02d/%02d',fsep_char,'%02d:%02d:%f',fsep_char,'%d',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%s',fsep_char,'%d',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%d',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f\n'); 
                        end 
          
                    end  %//if (strcmp(fsep_char,'default'))
    
                 
                    if state.isOutPosition()
                        
                        if state.isForwardBackwardKF()

                            fid_out_f = fopen([filerootOUT '_position_forward.txt'], 'wt');
                            fprintf(fid_out_f,'  %s \n' ,          'goGPS v0.5.2 beta 1');
                            fprintf(fid_out_f,' << %s >>\n' ,         'TROPOSPHERIC DELAY MODELLING OPTIONS');
                            fprintf(fid_out_f,'    ====================================\n');

                            fid_out_b = fopen([filerootOUT '_position_backward.txt'], 'wt');
                            fprintf(fid_out_b,'  %s \n' ,          'goGPS v0.5.2 beta 1');
                            fprintf(fid_out_b,' << %s >>\n' ,         'TROPOSPHERIC DELAY MODELLING OPTIONS');
                            fprintf(fid_out_b,'    ====================================\n');
                            
                            if option_Cmodels == 1 & option_Smodels == 0 %if user selects option for combine tropo model

                               fprintf(fid_out_f,tropomodel);%Tropospheric Delay Correction model
                               fprintf(fid_out_f,MF_model); %Mapping Function(MF) model
                               fprintf(fid_out_f,METsource); %source Meteorological parameters
                               fprintf(fid_out_f,TmModel); %Mean Temperature(Tm) model

                               fprintf(fid_out_b,tropomodel);%Tropospheric Delay Correction model
                               fprintf(fid_out_b,MF_model); %Mapping Function(MF) model
                               fprintf(fid_out_b,METsource); %source Meteorological parameters
                               fprintf(fid_out_b,TmModel); %Mean Temperature(Tm) model
                               
                            elseif option_Smodels == 1 & option_Cmodels == 0 %if user selects option for separate tropo model

                                  fprintf(fid_out_f,DRYmodel);%Hydrostatic Tropospheric Delay Correction model
                                  fprintf(fid_out_f,WETmodel);%Wet Tropospheric Delay Correction model
       
                                  fprintf(fid_out_b,DRYmodel);%Hydrostatic Tropospheric Delay Correction model
                                  fprintf(fid_out_b,WETmodel);%Wet Tropospheric Delay Correction model
       
                                  %GET MAPPING FUNCTION(MF) OPTION
                                  if option_different_MF==1 & option_model_MF==0 %if user selects option for use of different MF model
                                     
                                      fprintf(fid_out_f,hMF_model);%Hydrostatic MF model
                                      fprintf(fid_out_f,wMF_model);%Wet MF model
                                      
                                      fprintf(fid_out_b,hMF_model);%Hydrostatic MF model
                                      fprintf(fid_out_b,wMF_model);%Wet MF model
          
                                  elseif option_model_MF == 1 & option_different_MF == 0 %if user selects option for use of tropo models  MF
                                         fprintf(fid_out_f,hMF_model);%Hydrostatic MF model
                                         fprintf(fid_out_f,wMF_model);%Wet MF model
                                         
                                         fprintf(fid_out_b,hMF_model);%Hydrostatic MF model
                                         fprintf(fid_out_b,wMF_model);%Wet MF model
                                  end      
       
                                  if any([any([strcmpi(dryModel,wetModel),strfind(dryModel,wetModel)]),metPARA_typeW_D])
                                     fprintf(fid_out_f,METsourceWD);%source Meteorological parameters 
                                     fprintf(fid_out_b,METsourceWD);%source Meteorological parameters 
                                  else   
                                      fprintf(fid_out_f,METsourceD);%source Meteorological parameters 
                                      fprintf(fid_out_f,METsourceW);%source Meteorological parameters 
                                      
                                      fprintf(fid_out_b,METsourceD);%source Meteorological parameters 
                                      fprintf(fid_out_b,METsourceW);%source Meteorological parameters 
                                  end  
             
                                  fprintf(fid_out_f,TmModel); %Mean Temperature(Tm) model
                                  fprintf(fid_out_b,TmModel); %Mean Temperature(Tm) model
             
                            end  %//if option_Cmodels == 1 & option_Smodels == 0 

                             fprintf(fid_out_f, head_str);
                             fprintf(fid_out_f,head_str_unit);
                             fprintf(fid_out_f,head_str_);
                             
                             fprintf(fid_out_b, head_str);
                             fprintf(fid_out_b,head_str_unit);
                             fprintf(fid_out_b,head_str_);
                            
                             
                             %****************FILE WRITING 
                             
                            for i = 1 : nSol / 2
                                %file writing
                                if (pivot_OUT(i) ~= 0)
                                    if state.getForwardBackwardKF > 0
                                        
                                       %WHEN VMF ZENITH DELAYS, PWV & PPP DELAYS NOT ICLUDED IN THE OUTPUT
                                       if all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                                          
                                          fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i));

                                           %WHEN ONLY VMF ZENITH DELAYS ARE INCLUDED
                                       elseif all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==0])
                                          
                                              fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i),ZHDvmf(i), ZWDvmf(i), ZTDvmf(i));
                 
                                         %WHEN ONLY PWVs ARE INCLUDED
                                       elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0])
                                              
                                              fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i), Tm(i), PI(i), PWV(i));
                      
                                              %WHEN BOTH PWV & VMF ZTDs ARE INCLUDED TO THE OUTPUT
                                       elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==0]) 
                                             
                                              fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i),ZHDvmf(i), ZWDvmf(i), ZTDvmf(i), Tm(i), PI(i), PWV(i));
                 
                                           %WHEN VMF ZENITH DELAYS & PWV ARE NOT ICLUDED BUT PPP DELAYS ARE INCLUDED IN THE OUTPUT
                                       elseif all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==1 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                    
                                              fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i), ZWDppp(i), ZTDppp(i));
                                           
                                        %WHEN PWV IS NOT ICLUDED BUT PPP & VMF ZENITH DELAYS &  ARE INCLUDED IN THE OUTPUT
                                       elseif  all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==1])
                              
                                               fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i), ZWDppp(i), ZTDppp(i),ZHDvmf(i), ZWDvmf(i), ZTDvmf(i));
                                           
                                        %WHEN PWV & PPP ZENITH DELAYS ARE INCLUDED BUT VMF ZTDs ARE NOT ICLUDED  IN THE OUTPUT      
                                       elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==1])
                         
                                              fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i), ZWDppp(i), ZTDppp(i), Tm(i), PI(i), PWV(i));
                                           
                                       else %WHEN ALL(PWV & PPP ZENITH DELAYS,PWV) ARE INCLUDED IN THE OUTPUT
                                           
                                           fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i), ZWDppp(i), ZTDppp(i),ZHDvmf(i), ZWDvmf(i), ZTDvmf(i), Tm(i), PI(i), PWV(i));
                                           
                                       end %//if all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)

                                        ib = nSol + 1 - i;
                                    else
                                        ii = nSol/2 + i;
                                        ib = nSol/2 + 1 - i;
                                        
                                         %WHEN VMF ZENITH DELAYS, PWV & PPP DELAYS NOT ICLUDED IN THE OUTPUT
                                         if all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                                          
                                            fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ii), lam_KAL(ii), h_KAL(ii), X_KAL(ii), Y_KAL(ii), Z_KAL(ii), NORTH_UTM(ii), EAST_UTM(ii), h_ortho(ii), utm_zone(ii,:), nsat(ii), HDOP(ii), KHDOP(ii), NORTH_KAL(ii), EAST_KAL(ii), UP_KAL(ii), fixed_amb(ii), succ_rate(ii), ZHD(ii), ZWD(ii), ZTD(ii));

                                              %WHEN ONLY VMF ZENITH DELAYS ARE INCLUDED
                                         elseif all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==0])
                                          
                                                fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ii), lam_KAL(ii), h_KAL(ii), X_KAL(ii), Y_KAL(ii), Z_KAL(ii), NORTH_UTM(ii), EAST_UTM(ii), h_ortho(ii), utm_zone(ii,:), nsat(ii), HDOP(ii), KHDOP(ii), NORTH_KAL(ii), EAST_KAL(ii), UP_KAL(ii), fixed_amb(ii), succ_rate(ii), ZHD(ii), ZWD(ii), ZTD(ii),ZHDvmf(ii), ZWDvmf(ii), ZTDvmf(ii));
                 
                                              %WHEN ONLY PWVs ARE INCLUDED
                                         elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0])
                                              
                                                fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ii), lam_KAL(ii), h_KAL(ii), X_KAL(ii), Y_KAL(ii), Z_KAL(ii), NORTH_UTM(ii), EAST_UTM(ii), h_ortho(ii), utm_zone(ii,:), nsat(ii), HDOP(ii), KHDOP(ii), NORTH_KAL(ii), EAST_KAL(ii), UP_KAL(ii), fixed_amb(ii), succ_rate(ii), ZHD(ii), ZWD(ii), ZTD(ii), Tm(ii), PI(ii), PWV(ii));
                      
                                              %WHEN BOTH PWV & VMF ZTDs ARE INCLUDED TO THE OUTPUT
                                         elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==0]) 
                                             
                                                fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ii), lam_KAL(ii), h_KAL(ii), X_KAL(ii), Y_KAL(ii), Z_KAL(ii), NORTH_UTM(ii), EAST_UTM(ii), h_ortho(ii), utm_zone(ii,:), nsat(ii), HDOP(ii), KHDOP(ii), NORTH_KAL(ii), EAST_KAL(ii), UP_KAL(ii), fixed_amb(ii), succ_rate(ii), ZHD(ii), ZWD(ii), ZTD(ii),ZHDvmf(ii), ZWDvmf(ii), ZTDvmf(ii), Tm(ii), PI(ii), PWV(ii));
                 
                                              %WHEN VMF ZENITH DELAYS & PWV ARE NOT ICLUDED BUT PPP DELAYS ARE INCLUDED IN THE OUTPUT
                                         elseif all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==1 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                    
                                                fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ii), lam_KAL(ii), h_KAL(ii), X_KAL(ii), Y_KAL(ii), Z_KAL(ii), NORTH_UTM(ii), EAST_UTM(ii), h_ortho(ii), utm_zone(ii,:), nsat(ii), HDOP(ii), KHDOP(ii), NORTH_KAL(ii), EAST_KAL(ii), UP_KAL(ii), fixed_amb(ii), succ_rate(ii), ZHD(ii), ZWD(ii), ZTD(ii), ZWDppp(ii), ZTDppp(ii));
                                           
                                             %WHEN PWV IS NOT ICLUDED BUT PPP & VMF ZENITH DELAYS &  ARE INCLUDED IN THE OUTPUT
                                         elseif  all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==1])
                              
                                                 fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ii), lam_KAL(ii), h_KAL(ii), X_KAL(ii), Y_KAL(ii), Z_KAL(ii), NORTH_UTM(ii), EAST_UTM(ii), h_ortho(ii), utm_zone(ii,:), nsat(ii), HDOP(ii), KHDOP(ii), NORTH_KAL(ii), EAST_KAL(ii), UP_KAL(ii), fixed_amb(ii), succ_rate(ii), ZHD(ii), ZWD(ii), ZTD(ii), ZWDppp(ii), ZTDppp(ii),ZHDvmf(ii), ZWDvmf(ii), ZTDvmf(ii));
                                           
                                               %WHEN PWV & PPP ZENITH DELAYS ARE INCLUDED BUT VMF ZTDs ARE NOT ICLUDED  IN THE OUTPUT      
                                         elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==1])
                         
                                                fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ii), lam_KAL(ii), h_KAL(ii), X_KAL(ii), Y_KAL(ii), Z_KAL(ii), NORTH_UTM(ii), EAST_UTM(ii), h_ortho(ii), utm_zone(ii,:), nsat(ii), HDOP(ii), KHDOP(ii), NORTH_KAL(ii), EAST_KAL(ii), UP_KAL(ii), fixed_amb(ii), succ_rate(ii), ZHD(ii), ZWD(ii), ZTD(ii), ZWDppp(ii), ZTDppp(ii), Tm(ii), PI(ii), PWV(ii));
                                           
                                         else %WHEN ALL(PWV & PPP ZENITH DELAYS,PWV) ARE INCLUDED IN THE OUTPUT
                                           
                                              fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ii), lam_KAL(ii), h_KAL(ii), X_KAL(ii), Y_KAL(ii), Z_KAL(ii), NORTH_UTM(ii), EAST_UTM(ii), h_ortho(ii), utm_zone(ii,:), nsat(ii), HDOP(ii), KHDOP(ii), NORTH_KAL(ii), EAST_KAL(ii), UP_KAL(ii), fixed_amb(ii), succ_rate(ii), ZHD(ii), ZWD(ii), ZTD(ii), ZWDppp(ii), ZTDppp(ii),ZHDvmf(ii), ZWDvmf(ii), ZTDvmf(ii), Tm(ii), PI(ii), PWV(ii));
                                           
                                         end %//if all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)


                                    end %//if state.getForwardBackwardKF > 0
                                    
                                    
                                    %WHEN VMF ZENITH DELAYS, PWV & PPP DELAYS NOT ICLUDED IN THE OUTPUT
                                    if all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                                          
                                       fprintf(fid_out_b, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ib), lam_KAL(ib), h_KAL(ib), X_KAL(ib), Y_KAL(ib), Z_KAL(ib), NORTH_UTM(ib), EAST_UTM(ib), h_ortho(ib), utm_zone(ib,:), nsat(ib), HDOP(ib), KHDOP(ib), NORTH_KAL(ib), EAST_KAL(ib), UP_KAL(ib), fixed_amb(ib), succ_rate(ib), ZHD(ib), ZWD(ib), ZTD(ib));

                                         %WHEN ONLY VMF ZENITH DELAYS ARE INCLUDED
                                    elseif all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==0])
                                          
                                           fprintf(fid_out_b, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ib), lam_KAL(ib), h_KAL(ib), X_KAL(ib), Y_KAL(ib), Z_KAL(ib), NORTH_UTM(ib), EAST_UTM(ib), h_ortho(ib), utm_zone(ib,:), nsat(ib), HDOP(ib), KHDOP(ib), NORTH_KAL(ib), EAST_KAL(ib), UP_KAL(ib), fixed_amb(ib), succ_rate(ib), ZHD(ib), ZWD(ib), ZTD(ib),ZHDvmf(ib), ZWDvmf(ib), ZTDvmf(ib));
                 
                                          %WHEN ONLY PWVs ARE INCLUDED
                                    elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0])
                                              
                                           fprintf(fid_out_b, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ib), lam_KAL(ib), h_KAL(ib), X_KAL(ib), Y_KAL(ib), Z_KAL(ib), NORTH_UTM(ib), EAST_UTM(ib), h_ortho(ib), utm_zone(ib,:), nsat(ib), HDOP(ib), KHDOP(ib), NORTH_KAL(ib), EAST_KAL(ib), UP_KAL(ib), fixed_amb(ib), succ_rate(ib), ZHD(ib), ZWD(ib), ZTD(ib), Tm(ib), PI(ib), PWV(ib));
                      
                                              %WHEN BOTH PWV & VMF ZTDs ARE INCLUDED TO THE OUTPUT
                                         elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==0]) 
                                             
                                                fprintf(fid_out_b, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ib), lam_KAL(ib), h_KAL(ib), X_KAL(ib), Y_KAL(ib), Z_KAL(ib), NORTH_UTM(ib), EAST_UTM(ib), h_ortho(ib), utm_zone(ib,:), nsat(ib), HDOP(ib), KHDOP(ib), NORTH_KAL(ib), EAST_KAL(ib), UP_KAL(ib), fixed_amb(ib), succ_rate(ib), ZHD(ib), ZWD(ib), ZTD(ib),ZHDvmf(ib), ZWDvmf(ib), ZTDvmf(ib), Tm(ib), PI(ib), PWV(ib));
                 
                                          %WHEN VMF ZENITH DELAYS & PWV ARE NOT ICLUDED BUT PPP DELAYS ARE INCLUDED IN THE OUTPUT
                                    elseif all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==1 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                    
                                           fprintf(fid_out_b, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ib), lam_KAL(ib), h_KAL(ib), X_KAL(ib), Y_KAL(ib), Z_KAL(ib), NORTH_UTM(ib), EAST_UTM(ib), h_ortho(ib), utm_zone(ib,:), nsat(ib), HDOP(ib), KHDOP(ib), NORTH_KAL(ib), EAST_KAL(ib), UP_KAL(ib), fixed_amb(ib), succ_rate(ib), ZHD(ib), ZWD(ib), ZTD(ib), ZWDppp(ib), ZTDppp(ib));
                                           
                                          %WHEN PWV IS NOT ICLUDED BUT PPP & VMF ZENITH DELAYS &  ARE INCLUDED IN THE OUTPUT
                                    elseif  all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==1])
                              
                                            fprintf(fid_out_b, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ib), lam_KAL(ib), h_KAL(ib), X_KAL(ib), Y_KAL(ib), Z_KAL(ib), NORTH_UTM(ib), EAST_UTM(ib), h_ortho(ib), utm_zone(ib,:), nsat(ib), HDOP(ib), KHDOP(ib), NORTH_KAL(ib), EAST_KAL(ib), UP_KAL(ib), fixed_amb(ib), succ_rate(ib), ZHD(ib), ZWD(ib), ZTD(ib), ZWDppp(ib), ZTDppp(ib),ZHDvmf(ib), ZWDvmf(ib), ZTDvmf(ib));
                                           
                                          %WHEN PWV & PPP ZENITH DELAYS ARE INCLUDED BUT VMF ZTDs ARE NOT ICLUDED  IN THE OUTPUT      
                                    elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==1])
                         
                                           fprintf(fid_out_b, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ib), lam_KAL(ib), h_KAL(ib), X_KAL(ib), Y_KAL(ib), Z_KAL(ib), NORTH_UTM(ib), EAST_UTM(ib), h_ortho(ib), utm_zone(ib,:), nsat(ib), HDOP(ib), KHDOP(ib), NORTH_KAL(ib), EAST_KAL(ib), UP_KAL(ib), fixed_amb(ib), succ_rate(ib), ZHD(ib), ZWD(ib), ZTD(ib), ZWDppp(ib), ZTDppp(ib), Tm(ib), PI(ib), PWV(ib));
                                           
                                    else %WHEN ALL(PWV & PPP ZENITH DELAYS,PWV) ARE INCLUDED IN THE OUTPUT
                                           
                                         fprintf(fid_out_b, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(ib), lam_KAL(ib), h_KAL(ib), X_KAL(ib), Y_KAL(ib), Z_KAL(ib), NORTH_UTM(ib), EAST_UTM(ib), h_ortho(ib), utm_zone(ib,:), nsat(ib), HDOP(ib), KHDOP(ib), NORTH_KAL(ib), EAST_KAL(ib), UP_KAL(ib), fixed_amb(ib), succ_rate(ib), ZHD(ib), ZWD(ib), ZTD(ib), ZWDppp(ib), ZTDppp(ib),ZHDvmf(ib), ZWDvmf(ib), ZTDvmf(ib), Tm(ii), PI(ii), PWV(ii));
                                           
                                    end  %//if all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)

                                else 
                                    
                                    %WHEN VMF ZENITH DELAYS, PWV & PPP DELAYS NOT ICLUDED IN THE OUTPUT
                                    if any([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==1 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                                          
                                       fprintf(fid_out_f, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i));
                                       fprintf(fid_out_f, '\n');
                                       fprintf(fid_out_b, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i));
                                       fprintf(fid_out_b, '\n'); 
                                    end  

                                end %//if (pivot_OUT(i) ~= 0)
                                
                            end %//for i = 1 : nSol / 2
                            
                            fclose(fid_out_f);
                            fclose(fid_out_b);
                            
                        else
                            fid_out = fopen([filerootOUT '_position.txt'], 'wt');
                            fprintf(fid_out,'  %s \n' ,          'goGPS v0.5.2 beta 1');
                            fprintf(fid_out,' << %s >>\n' ,         'TROPOSPHERIC DELAY MODELLING OPTIONS');
                            fprintf(fid_out,'    ====================================\n');

      
                            if option_Cmodels == 1 & option_Smodels == 0 %if user selects option for combine tropo model

                               fprintf(fid_out,tropomodel);%Tropospheric Delay Correction model
                               fprintf(fid_out,MF_model); %Mapping Function(MF) model
                               fprintf(fid_out,METsource); %source Meteorological parameters
                               fprintf(fid_out,TmModel); %Mean Temperature(Tm) model

                            elseif option_Smodels == 1 & option_Cmodels == 0 %if user selects option for separate tropo model

                                   fprintf(fid_out,DRYmodel);%Hydrostatic Tropospheric Delay Correction model
                                   fprintf(fid_out,WETmodel);%Wet Tropospheric Delay Correction model
       
                                   %GET MAPPING FUNCTION(MF) OPTION
                                   if option_different_MF==1 & option_model_MF==0 %if user selects option for use of different MF model
                                      fprintf(fid_out,hMF_model);%Hydrostatic MF model
                                      fprintf(fid_out,wMF_model);%Wet MF model
          
                                   elseif option_model_MF == 1 & option_different_MF == 0 %if user selects option for use of tropo models  MF
                                          fprintf(fid_out,hMF_model);%Hydrostatic MF model
                                          fprintf(fid_out,wMF_model);%Wet MF model
                                   end      
       
                                  if any([any([strcmpi(dryModel,wetModel),strfind(dryModel,wetModel)]),metPARA_typeW_D])
                                     fprintf(fid_out,METsourceWD);%source Meteorological parameters 
                                  else   
                                      fprintf(fid_out,METsourceD);%source Meteorological parameters 
                                      fprintf(fid_out,METsourceW);%source Meteorological parameters 
                                  end  
             
                                  fprintf(fid_out,TmModel); %Mean Temperature(Tm) model
             
                            end %//if option_Cmodels == 1 & option_Smodels == 0 

                            fprintf(fid_out, head_str);
                            fprintf(fid_out,head_str_unit);
                            fprintf(fid_out,head_str_);

                            for i = 1 : nSol
                                %file writing
                                if (pivot_OUT(i) ~= 0 || state.isModeBlock())
                                    
                                    if state.isModeBlock()
                                        
                                       %WHEN VMF ZENITH DELAYS, PWV & PPP DELAYS NOT ICLUDED IN THE OUTPUT
                                       if all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                                          
                                          fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat, 0, 0, NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i));

                                           %WHEN ONLY VMF ZENITH DELAYS ARE INCLUDED
                                       elseif all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==0])
                                          
                                              fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat, 0, 0, NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i),ZHDvmf(i), ZWDvmf(i), ZTDvmf(i));
                 
                                         %WHEN ONLY PWVs ARE INCLUDED
                                       elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0])
                                              
                                              fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat, 0, 0, NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i), Tm(i), PI(i), PWV(i));
                      
                                              %WHEN BOTH PWV & VMF ZTDs ARE INCLUDED TO THE OUTPUT
                                       elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==0]) 
                                             
                                              fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat, 0, 0, NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i),ZHDvmf(i), ZWDvmf(i), ZTDvmf(i), Tm(i), PI(i), PWV(i));
                 
                                           %WHEN VMF ZENITH DELAYS & PWV ARE NOT ICLUDED BUT PPP DELAYS ARE INCLUDED IN THE OUTPUT
                                       elseif all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==1 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                    
                                              fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat, 0, 0, NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i), ZWDppp(i), ZTDppp(i));
                                           
                                        %WHEN PWV IS NOT ICLUDED BUT PPP & VMF ZENITH DELAYS &  ARE INCLUDED IN THE OUTPUT
                                       elseif  all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==1])
                              
                                               fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat, 0, 0, NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i), ZWDppp(i), ZTDppp(i),ZHDvmf(i), ZWDvmf(i), ZTDvmf(i));
                                           
                                        %WHEN PWV & PPP ZENITH DELAYS ARE INCLUDED BUT VMF ZTDs ARE NOT ICLUDED  IN THE OUTPUT      
                                       elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==1])
                         
                                              fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat, 0, 0, NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i), ZWDppp(i), ZTDppp(i), Tm(i), PI(i), PWV(i));
                                           
                                       else  %WHEN ALL(PWV & PPP ZENITH DELAYS,PWV) ARE INCLUDED IN THE OUTPUT
                                           
                                           fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat, 0, 0, NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i), ZWDppp(i), ZTDppp(i),ZHDvmf(i), ZWDvmf(i), ZTDvmf(i), Tm(i), PI(i), PWV(i));
                                           
                                       end  %//if all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
 
                                    else
                                        
                                        %WHEN VMF ZENITH DELAYS, PWV & PPP DELAYS NOT ICLUDED IN THE OUTPUT
                                        if all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                                          
                                           fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i));

                                           %WHEN ONLY VMF ZENITH DELAYS ARE INCLUDED
                                        elseif all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==0])
                                          
                                               fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i),ZHDvmf(i), ZWDvmf(i), ZTDvmf(i));
                 
                                         %WHEN ONLY PWVs ARE INCLUDED
                                        elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0])
                                              
                                               fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i), Tm(i), PI(i), PWV(i));
                      
                                              %WHEN BOTH PWV & VMF ZTDs ARE INCLUDED TO THE OUTPUT
                                        elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==0]) 
                                             
                                               fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i),ZHDvmf(i), ZWDvmf(i), ZTDvmf(i), Tm(i), PI(i), PWV(i));
                 
                                           %WHEN VMF ZENITH DELAYS & PWV ARE NOT ICLUDED BUT PPP DELAYS ARE INCLUDED IN THE OUTPUT
                                        elseif all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==1 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                    
                                               fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i), ZWDppp(i), ZTDppp(i));
                                           
                                        %WHEN PWV IS NOT ICLUDED BUT PPP & VMF ZENITH DELAYS &  ARE INCLUDED IN THE OUTPUT
                                        elseif  all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==1])
                              
                                                fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i), ZWDppp(i), ZTDppp(i),ZHDvmf(i), ZWDvmf(i), ZTDvmf(i));
                                           
                                        %WHEN PWV & PPP ZENITH DELAYS ARE INCLUDED BUT VMF ZTDs ARE NOT ICLUDED  IN THE OUTPUT      
                                        elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==1])
                         
                                               fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i), ZWDppp(i), ZTDppp(i), Tm(i), PI(i), PWV(i));
                                           
                                        else %WHEN ALL(PWV & PPP ZENITH DELAYS,PWV) ARE INCLUDED IN THE OUTPUT
                                           
                                             fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), doy(i), phi_KAL(i), lam_KAL(i), h_KAL(i), X_KAL(i), Y_KAL(i), Z_KAL(i), NORTH_UTM(i), EAST_UTM(i), h_ortho(i), utm_zone(i,:), nsat(i), HDOP(i), KHDOP(i), NORTH_KAL(i), EAST_KAL(i), UP_KAL(i), fixed_amb(i), succ_rate(i), ZHD(i), ZWD(i), ZTD(i), ZWDppp(i), ZTDppp(i),ZHDvmf(i), ZWDvmf(i), ZTDvmf(i), Tm(i), PI(i), PWV(i));
                                           
                                        end  %//if all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)

                                    end  %//if state.isModeBlock()
                                    
                                else
                                     if any([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==1 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                                    
                                        fprintf(fid_out, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), tow(i));
                                        fprintf(fid_out, '\n');
                                     end
                                end %//if (pivot_OUT(i) ~= 0 || state.isModeBlock())
                                
                            end %//for i = 1 : nSol
                            
                            fclose(fid_out);
                            
                        end %//if state.isForwardBackwardKF()
                        
                        
                       %*********WRITE OUTPUT TO CSV
                       %WHEN VMF ZENITH DELAYS, PWV & PPP DELAYS ARE NOT ICLUDED IN THE OUTPUT
                       if all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
          
                          %CREATE COLUMN HEADERS FOR ALL EXTRACTED DATA
                          header1='MJD,GPS week,Doy,Latitude,Longitude,h (ellips.),ECEF X,ECEF Y,ECEF Z,UTM North,UTM East,h (orthom.),ZHD,ZWD,ZTD '; 
                          header2=',,,[deg],[deg],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m]';
         
                          %DATA
                          CSVout= [MJD week_R doy [phi_KAL,lam_KAL, h_KAL],[X_KAL, Y_KAL, Z_KAL], round([NORTH_UTM, EAST_UTM, h_ortho'],4),round(ZHD,4), round(ZWD,4), round(ZTD,4)];     
              
                           %WHEN ONLY VMF ZENITH DELAYS ARE ICLUDED IN THE OUTPUT
                       elseif all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==0])
                 
                              %CREATE COLUMN HEADERS FOR ALL EXTRACTED DATA
                              header1='MJD,GPS week,Doy,Latitude,Longitude,h (ellips.),ECEF X,ECEF Y,ECEF Z,UTM North,UTM East,h (orthom.),ZHD,ZWD,ZTD,ZHDvmf,ZWDvmf,ZTDvmf '; 
                              header2=',,,[deg],[deg],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m]';
         
                              %DATA
                              CSVout= [MJD week_R doy [phi_KAL,lam_KAL, h_KAL],[X_KAL, Y_KAL, Z_KAL], round([NORTH_UTM, EAST_UTM, h_ortho'],4),round(ZHD,4), round(ZWD,4), round(ZTD,4),round(ZHDvmf,4), round(ZWDvmf,4), round(ZTDvmf,4)];     
             
                              %WHEN ONLY PWVs ARE ICLUDED IN THE OUTPUT  
                       elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==0])
              
                              %CREATE COLUMN HEADERS FOR ALL EXTRACTED DATA
                              header1='MJD,GPS week,Doy,Latitude,Longitude,h (ellips.),ECEF X,ECEF Y,ECEF Z,UTM North,UTM East,h (orthom.),ZHD,ZWD,ZTD,Tm,PI,PWV '; 
                              header2=',,,[deg],[deg],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[K],,[mm]';
         
                              %DATA
                              CSVout= [MJD week_R doy [phi_KAL,lam_KAL, h_KAL],[X_KAL, Y_KAL, Z_KAL], round([NORTH_UTM, EAST_UTM, h_ortho'],4),round(ZHD,4), round(ZWD,4), round(ZTD,4),round(Tm,3),round(PI,5),round(PWV,4)];     
                  
                            %WHEN BOTH VMF ZENITH DELAYS & PWV ARE ICLUDED IN THE OUTPUT
                       elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==0]) %WHEN BOTH PWV & VMF ZTDs ARE INCLUDED TO THE OUTPUT
               
                              %CREATE COLUMN HEADERS FOR ALL EXTRACTED DATA
                              header1='MJD,GPS week,Doy,Latitude,Longitude,h (ellips.),ECEF X,ECEF Y,ECEF Z,UTM North,UTM East,h (orthom.),ZHD,ZWD,ZTD,ZHDvmf,ZWDvmf,ZTDvmf,Tm,PI,PWV '; 
                              header2=',,,[deg],[deg],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[K],,[mm]';
         
                              %DATA
                              CSVout= [MJD week_R doy [phi_KAL,lam_KAL, h_KAL],[X_KAL, Y_KAL, Z_KAL], round([NORTH_UTM, EAST_UTM, h_ortho'],4),round(ZHD,4), round(ZWD,4), round(ZTD,4),round(ZHDvmf,4), round(ZWDvmf,4), round(ZTDvmf,4),round(Tm,3),round(PI,5),round(PWV,4)];  
                              
                               %WHEN VMF ZENITH DELAYS & PWV ARE NOT ICLUDED BUT PPP DELAYS ARE INCLUDED IN THE OUTPUT
                       elseif all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==1 ]) %GET PWV STATE(1 TO RETRIEVE & 0 NOT TO RETRIEVE)
                 
                              %CREATE COLUMN HEADERS FOR ALL EXTRACTED DATA
                              header1='MJD,GPS week,Doy,Latitude,Longitude,h (ellips.),ECEF X,ECEF Y,ECEF Z,UTM North,UTM East,h (orthom.),ZHD,ZWD,ZTD,ZWDppp,ZTDppp '; 
                              header2=',,,[deg],[deg],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m]';
         
                              %DATA
                              CSVout= [MJD week_R doy [phi_KAL,lam_KAL, h_KAL],[X_KAL, Y_KAL, Z_KAL], round([NORTH_UTM, EAST_UTM, h_ortho'],4),round(ZHD,4), round(ZWD,4), round(ZTD,4), round(ZWDppp,4), round(ZTDppp,4)];     
             
                               %WHEN PWV IS NOT ICLUDED BUT PPP & VMF ZENITH DELAYS &  ARE INCLUDED IN THE OUTPUT
                       elseif  all([getappdata(0,'pwv_state') == 0,getappdata(0,'getvmfZTDs_state')== 1,state.isTropoOn()==1])
                           
                               %CREATE COLUMN HEADERS FOR ALL EXTRACTED DATA
                               header1='MJD,GPS week,Doy,Latitude,Longitude,h (ellips.),ECEF X,ECEF Y,ECEF Z,UTM North,UTM East,h (orthom.),ZHD,ZWD,ZTD,ZWDppp,ZTDppp,ZHDvmf,ZWDvmf,ZTDvmf '; 
                               header2=',,,[deg],[deg],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m]';
         
                               %DATA
                               CSVout= [MJD week_R doy [phi_KAL,lam_KAL, h_KAL],[X_KAL, Y_KAL, Z_KAL], round([NORTH_UTM, EAST_UTM, h_ortho'],4),round(ZHD,4), round(ZWD,4), round(ZTD,4), round(ZWDppp,4), round(ZTDppp,4),round(ZHDvmf,4), round(ZWDvmf,4), round(ZTDvmf,4)];     
             
                              %WHEN PWV & PPP ZENITH DELAYS ARE INCLUDED BUT VMF ZTDs ARE NOT ICLUDED  IN THE OUTPUT      
                       elseif all([getappdata(0,'pwv_state') == 1,getappdata(0,'getvmfZTDs_state')== 0,state.isTropoOn()==1])
                      
                              %CREATE COLUMN HEADERS FOR ALL EXTRACTED DATA
                              header1='MJD,GPS week,Doy,Latitude,Longitude,h (ellips.),ECEF X,ECEF Y,ECEF Z,UTM North,UTM East,h (orthom.),ZHD,ZWD,ZTD,ZWDppp,ZTDppp,Tm,PI,PWV '; 
                              header2=',,,[deg],[deg],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[K],,[mm]';
         
                              %DATA
                              CSVout= [MJD week_R doy [phi_KAL,lam_KAL, h_KAL],[X_KAL, Y_KAL, Z_KAL], round([NORTH_UTM, EAST_UTM, h_ortho'],4),round(ZHD,4), round(ZWD,4), round(ZTD,4), round(ZWDppp,4), round(ZTDppp,4),round(Tm,3),round(PI,5),round(PWV,4)];     
                  
                       else %WHEN ALL(PWV & PPP ZENITH DELAYS,PWV) ARE INCLUDED IN THE OUTPUT
                              
                            %CREATE COLUMN HEADERS FOR ALL EXTRACTED DATA
                            header1='MJD,GPS week,Doy,Latitude,Longitude,h (ellips.),ECEF X,ECEF Y,ECEF Z,UTM North,UTM East,h (orthom.),ZHD,ZWD,ZTD,ZWDppp,ZTDppp,ZHDvmf,ZWDvmf,ZTDvmf,Tm,PI,PWV '; 
                            header2=',,,[deg],[deg],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[m],[K],,[mm]';
         
                            %DATA
                            CSVout= [MJD week_R doy [phi_KAL,lam_KAL, h_KAL],[X_KAL, Y_KAL, Z_KAL], round([NORTH_UTM, EAST_UTM, h_ortho'],4),round(ZHD,4), round(ZWD,4), round(ZTD,4), round(ZWDppp,4), round(ZTDppp,4),round(ZHDvmf,4), round(ZWDvmf,4), round(ZTDvmf,4),round(Tm,3),round(PI,5),round(PWV,4)];     
             
                       end
                       
                       %CREATE DIRECTORY 
                       fileoutput=fullfile([filerootOUT '_pos_ZTDs.csv']);
   
                       %OPEN FILE TO WRITE
                       fid_csv=fopen(fileoutput,'wt');

                       fprintf(fid_csv,'%s\n',header1);
                       fprintf(fid_csv,'%s',header2);
                       fclose(fid_csv);

                       dlmwrite (fileoutput,CSVout,'-append','delimiter',',','precision','%.8f','roffset',1)

                       %OPEN OUTPUT FILE IN WINDOWS
                       winopen([filerootOUT '_position.txt'])
                        

                       %IF USER WANTS TO OUTPUT OTHER REPORTS SUCH SAVING IN '.MAT' matlab format
                       if getappdata(0,'report_option') == 1
                           
                        if ~is_batch
                            % Save in matlab format all the outputs
                            if flag_tropo
                                
                                if getappdata(0,'pwv_state') == 1 %WHEN PWV CHECKBOX IS CHECKED
                                    
                                   save([filerootOUT '_position.mat'], 'date_R', 'week_R', 'tow', 'phi_KAL', 'lam_KAL', 'h_KAL', 'X_KAL', 'Y_KAL', 'Z_KAL', 'NORTH_UTM', 'EAST_UTM', 'utm_zone', 'nsat', 'HDOP', 'KHDOP', 'NORTH_KAL', 'EAST_KAL', 'UP_KAL', 'fixed_amb', 'succ_rate', 'estim_tropo', 'ZWD', 'PWV');
                                else
                                    save([filerootOUT '_position.mat'], 'date_R', 'week_R', 'tow', 'phi_KAL', 'lam_KAL', 'h_KAL', 'X_KAL', 'Y_KAL', 'Z_KAL', 'NORTH_UTM', 'EAST_UTM', 'utm_zone', 'nsat', 'HDOP', 'KHDOP', 'NORTH_KAL', 'EAST_KAL', 'UP_KAL', 'fixed_amb', 'succ_rate', 'estim_tropo', 'ZWD');
                                end
                            else
                                save([filerootOUT '_position.mat'], 'date_R', 'week_R', 'tow', 'phi_KAL', 'lam_KAL', 'h_KAL', 'X_KAL', 'Y_KAL', 'Z_KAL', 'NORTH_UTM', 'EAST_UTM', 'utm_zone', 'nsat', 'HDOP', 'KHDOP', 'NORTH_KAL', 'EAST_KAL', 'UP_KAL', 'fixed_amb', 'succ_rate');
                            
                            end %//if flag_tropo 
                            
                        end %//if ~is_batch
                        
                       end %//if getappdata(0,'report_option') == 1
                        
                    end %//if state.isOutPosition()
                    
                end %//if goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)
                
                %----------------------------------------------------------------------------------------------
                %% REPORT FILE (PDF)
                %----------------------------------------------------------------------------------------------
                
                %if any positioning was done (either post-processing or real-time)
                if (goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)) && (~isempty(EAST))
                    %display information
                    fprintf('Writing report file (PDF)...\n');
                    
                    if state.isOutDtPDF()
                        
                        if (exist('dtR','var'))
                            f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A4','Visible','off');
                            paperSize = get(f,'PaperSize');
                            set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
                            plot(dtR.*goGNSS.V_LIGHT,'.r');
                            grid on;
                            set(gca,'FontName','Verdana');
                            set(gca,'FontSize',10);
                            xlabel('Epoch','FontName','Verdana','FontSize',10,'FontWeight','Bold');
                            ylabel('Epoch-by-epoch code-estimated receiver clock (m)','FontName','Verdana','FontSize',10,'FontWeight','Bold');
                            %print PDF
                            print(f, '-dpdf', [filerootOUT '_dtR']);
                            close(f)
                        end
                    end
                    
                    if (mode == goGNSS.MODE_PP_KF_CP_SA)
                        
                        dtR_KAL = Xhat_t_t_OUT(end-nC+1,:)./goGNSS.V_LIGHT;
                        
                        if state.isOutDtPDF()
                            f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A4','Visible','off');
                            paperSize = get(f,'PaperSize');
                            set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
                            plot(dtR_KAL.*goGNSS.V_LIGHT,'.r');
                            grid on;
                            set(gca,'FontName','Verdana');
                            set(gca,'FontSize',10);
                            xlabel('Epoch','FontName','Verdana','FontSize',10,'FontWeight','Bold');
                            ylabel('Kalman-estimated receiver clock (m)','FontName','Verdana','FontSize',10,'FontWeight','Bold');
                            %print PDF
                            print(f, '-dpdf', [filerootOUT '_dtR_KAL']);
                            close(f)
                        end
                        
                        if (~isempty(ISBs) && state.isOutReportPDF)
                            sys = unique(constellations.systems);
                            for ISB = (nC-1) : -1 : 0
                                ISB_KAL = Xhat_t_t_OUT(end-ISB,:)./goGNSS.V_LIGHT;
                                switch sys(end-ISB)
                                    case 'R'
                                        sys_label = 'GLONASS';
                                    case 'E'
                                        sys_label = 'Galileo';
                                    case 'C'
                                        sys_label = 'BeiDou';
                                    case 'J'
                                        sys_label = 'QZSS';
                                end
                                
                                f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A4','Visible','off');
                                paperSize = get(f,'PaperSize');
                                set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
                                plot(ISB_KAL.*goGNSS.V_LIGHT,'.r');
                                grid on;
                                set(gca,'FontName','Verdana');
                                set(gca,'FontSize',10);
                                xlabel('Epoch','FontName','Verdana','FontSize',10,'FontWeight','Bold');
                                ylabel(['Kalman-estimated ' sys_label ' inter-system bias (m)'],'FontName','Verdana','FontSize',10,'FontWeight','Bold');
                                %print PDF
                                print(f, '-dpdf', [filerootOUT '_' sys_label '_ISB_KAL']);
                                close(f)
                            end
                        end
                    end
                    
                    if (any(estim_tropo) && state.isOutTropoPDF())
                        f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A4','Visible','off');
                        paperSize = get(f,'PaperSize');
                        set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
                        plot(estim_tropo,'.r');
                        grid on;
                        set(gca,'FontName','Verdana');
                        set(gca,'FontSize',10);
                        xlabel('Epoch','FontName','Verdana','FontSize',10,'FontWeight','Bold');
                        ylabel('Zenith Total Delay (m)','FontName','Verdana','FontSize',10,'FontWeight','Bold');
                        %print PDF
                        print(f, '-dpdf', [filerootOUT '_tropo']);
                        close(f)
                    end
                    
                    if goGNSS.isDD(mode && state.isOutDtPDF)
                        f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A4','Visible','off');
                        paperSize = get(f,'PaperSize');
                        set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
                        plot(dtM.*goGNSS.V_LIGHT,'.r');
                        set(gca,'FontName','Verdana');
                        set(gca,'FontSize',7);
                        ylabel('Receiver clock (m)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
                        %print PDF
                        print(f, '-dpdf', [filerootOUT '_dtM']);
                        close(f)
                    end
                    
                    if (~state.isModeBlock())
                        if (any(RES_PHASE1_FIXED(:)))
                            RES_PHASE1 = RES_PHASE1_FIXED;
                            RES_CODE1  = RES_CODE1_FIXED;
                            RES_PHASE2 = RES_PHASE2_FIXED;
                            RES_CODE2  = RES_CODE2_FIXED;
                        else
                            RES_PHASE1 = RES_PHASE1_FLOAT;
                            RES_CODE1  = RES_CODE1_FLOAT;
                            RES_PHASE2 = RES_PHASE2_FLOAT;
                            RES_CODE2  = RES_CODE2_FLOAT;
                        end
                    else
                        RES_PHASE1_FLOAT_MELSA = go_block.getPhaseResiduals()';
                        RES_PHASE1_FLOAT_MELSA(:, go_block.getEmptyEpochs()) = [];
                        RES_CODE1_FLOAT_MELSA = zeros(size(RES_PHASE1_FLOAT_MELSA));
                        RES_PHASE2_FLOAT_MELSA = zeros(size(RES_PHASE1_FLOAT_MELSA));
                        RES_CODE2_FLOAT_MELSA = zeros(size(RES_PHASE1_FLOAT_MELSA));
                        
                        RES_PHASE1 = RES_PHASE1_FLOAT_MELSA;
                        RES_CODE1  = RES_CODE1_FLOAT_MELSA;
                        RES_PHASE2 = RES_PHASE2_FLOAT_MELSA;
                        RES_CODE2  = RES_CODE2_FLOAT_MELSA;
                        RES_PHASE1(RES_PHASE1==0) = NaN;
                        RES_CODE1(RES_CODE1==0) = NaN;
                        RES_PHASE2(RES_PHASE2==0) = NaN;
                        RES_CODE2(RES_CODE2==0) = NaN;
                        outliers_PHASE1 = zeros(size(RES_PHASE1));
                        outliers_CODE1 = zeros(size(RES_CODE1));
                        outliers_PHASE2 = zeros(size(RES_PHASE2));
                        outliers_CODE2 = zeros(size(RES_CODE2));
                    end
                    
                    if state.isOutPhaseResPDF() || state.isOutCodeResPDF()
                        if (length(frequencies) > 1)
                            if (~strcmp(obs_comb, 'IONO_FREE'))
                                plot_residuals(constellations, RES_PHASE1, RES_CODE1, outliers_PHASE1, outliers_CODE1, [filerootOUT '_L' num2str(frequencies(1))]);
                                plot_residuals(constellations, RES_PHASE2, RES_CODE2, outliers_PHASE2, outliers_CODE2, [filerootOUT '_L' num2str(frequencies(2))]);
                            else
                                plot_residuals(constellations, RES_PHASE1, RES_CODE1, outliers_PHASE1, outliers_CODE1, [filerootOUT '_IF']);
                            end
                        else
                            plot_residuals(constellations, RES_PHASE1, RES_CODE1, outliers_PHASE1, outliers_CODE1, [filerootOUT '_L' num2str(frequencies(1))]);
                        end
                    end
                    
                    if state.isOutReportPDF()
                        f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','portrait','PaperUnits','centimeters','PaperType','A4','Visible','off');
                        paperSize = get(f,'PaperSize');
                        set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
                        
                        %settings
                        f1 = subplot(7,3,1);
                        set(f1,'Visible','off');
                        switch mode
                            case 1
                                text(0,1.00,sprintf('Mode: code\n        undifferenced'));
                                text(0,0.75,sprintf('Kalman filter: no'));
                            case 2
                                text(0,1.00,sprintf('Mode: code\n        undifferenced'));
                                text(0,0.75,sprintf('Kalman filter: yes'));
                            case 3
                                text(0,1.00,sprintf('Mode: code and phase\n        undifferenced'));
                                text(0,0.75,sprintf('Kalman filter: no'));
                            case 4
                                text(0,1.00,sprintf('Mode: code and phase\n        undifferenced'));
                                text(0,0.75,sprintf('Kalman filter: yes'));
                            case 11
                                text(0,1.00,sprintf('Mode: code\n        double difference'));
                                text(0,0.75,sprintf('Kalman filter: no'));
                            case 12
                                text(0,1.00,sprintf('Mode: code\n        double difference'));
                                text(0,0.75,sprintf('Kalman filter: yes'));
                            case 13
                                text(0,1.00,sprintf('Mode: code and phase\n        double difference'));
                                text(0,0.75,sprintf('Kalman filter: no'));
                            case 14
                                text(0,1.00,sprintf('Mode: code and phase\n        double difference'));
                                text(0,0.75,sprintf('Kalman filter: yes'));
                            case 24
                                text(0,1.00,sprintf('Mode: code and phase\n        double difference'));
                                text(0,0.75,sprintf('Kalman filter: yes'));
                        end
                        if (mode == goGNSS.MODE_PP_KF_C_SA) || (mode == goGNSS.MODE_PP_KF_CP_SA) || (mode == goGNSS.MODE_PP_KF_C_DD) || (mode == goGNSS.MODE_PP_KF_CP_DD) || (mode == goGNSS.MODE_RT_NAV)
                            switch order
                                case 1
                                    text(0,0.50,sprintf('Dynamics: static'));
                                case 2
                                    text(0,0.50,sprintf('Dynamics: constant\n            velocity'));
                                case 3
                                    text(0,0.50,sprintf('Dynamics: constant\n            acceleration'));
                            end
                        end
                        text(0,0.25,sprintf('Cutoff: %d deg', cutoff));
                        switch weights
                            case 0
                                text(0,0,sprintf('Weights: no weights'));
                            case 1
                                text(0,0,sprintf('Weights: elevation\n           (1/sin(el)^2)'));
                            case 2
                                text(0,0,sprintf('Weights: elevation\n           (exp)'));
                            case 3
                                text(0,0,sprintf('Weights: SNR'));
                            case 4
                                text(0,0,sprintf('Weights: elevation\n           and SNR'));
                        end
                        
                        %trajectory plotting
                        f3 = subplot(7,3,[2 3 5 6 8 9 11 12]);
                        %if relative positioning (i.e. with master station)
                        if (goGNSS.isDD(mode) || mode == goGNSS.MODE_RT_NAV)
                            EAST_O = EAST_KAL(end); NORTH_O = NORTH_KAL(end);
                        else
                            EAST_O = EAST_UTM(1); NORTH_O = NORTH_UTM(1);
                        end
                        plot(EAST-EAST_O, NORTH-NORTH_O, '.r');
                        axis equal
                        xlabel('EAST [m]'); ylabel('NORTH [m]'); grid on;
                        hold on
                        if (o1 == 1) && ((mode == goGNSS.MODE_PP_KF_C_SA) || (mode == goGNSS.MODE_PP_KF_CP_SA) || (mode == goGNSS.MODE_PP_KF_C_DD) || (mode == goGNSS.MODE_PP_KF_CP_DD) || (mode == goGNSS.MODE_RT_NAV))
                            %static positioning solution plotting
                            plot(EAST(end)-EAST_O, NORTH(end)-NORTH_O, '*b');
                            %        %covariance propagation
                            %        Cee_ENU = global2localCov(Cee_OUT([1 o1+1 o2+1],[1 o1+1 o2+1],end), Xhat_t_t_OUT([1 o1+1 o2+1],end));
                            
                            legend('Positioning (KF)','Final position','Location','SouthOutside');
                        else
                            legend('Positioning','Location','SouthOutside');
                        end
                        
                        if (mode == goGNSS.MODE_PP_LS_C_SA || mode == goGNSS.MODE_PP_LS_CP_SA || mode == goGNSS.MODE_PP_LS_CP_DD_L || mode == goGNSS.MODE_PP_LS_C_DD)
                            EAST_R = mean(EAST);
                            NORTH_R = mean(NORTH);
                            h_R = mean(h_KAL);
                            plot(EAST_R-EAST_O, NORTH_R-NORTH_O, '*b');
                        end
                        
                        %coordinate transformation (UTM)
                        [EAST_REF, NORTH_REF, h_REF, utm_zone] = cart2plan(pos_REF(1,1), pos_REF(2,1), pos_REF(3,1));
                        
                        %plot(EAST_M-EAST_O, NORTH_M-NORTH_O, 'xc', 'LineWidth', 2);
                        
                        %statistics
                        f2 = subplot(7,3,[4 7 10]);
                        set(f2,'Visible','off');
                        if (o1 == 1) && ((mode == goGNSS.MODE_PP_KF_C_SA) || (mode == goGNSS.MODE_PP_KF_CP_SA) || (mode == goGNSS.MODE_PP_KF_C_DD) || (mode == goGNSS.MODE_PP_KF_CP_DD) || (mode == goGNSS.MODE_RT_NAV))
                            text(0,0.95,'----------------');
                            text(0,0.90,'Final position (UTM)');
                            text(0,0.83,sprintf('E: %.3f m', EAST_UTM(end)));
                            text(0,0.78,sprintf('N: %.3f m', NORTH_UTM(end)));
                            text(0,0.73,sprintf('h(ell.): %.3f m', h_KAL(end)));
                            %         text(0,0.62,'----------------');
                            %         text(0,0.57,'Position estimation error');
                            %         text(0,0.50,sprintf('E: %.4f m', Cee_ENU(1,1)));
                            %         text(0,0.45,sprintf('N: %.4f m', Cee_ENU(2,2)));
                            %         text(0,0.40,sprintf('U: %.4f m', Cee_ENU(3,3)));
                        end
                        
                        if (mode == goGNSS.MODE_PP_LS_C_DD || mode == goGNSS.MODE_PP_LS_CP_DD_L)
                            text(0,0.95,'----------------');
                            text(0,0.90,'Baseline (mean)');
                            text(0,0.83,sprintf('E: %.3f m', mean(EAST_KAL)));
                            text(0,0.78,sprintf('N: %.3f m', mean(NORTH_KAL)));
                            text(0,0.73,sprintf('h: %.3f m', mean(UP_KAL)));
                        else
                            %if relative positioning (i.e. with master station)
                            if (goGNSS.isDD(mode) || mode == goGNSS.MODE_RT_NAV)
                                text(0,0.62,'----------------');
                                text(0,0.57,'Plot false origin');
                                text(0,0.52,'(last estimated baseline)');
                            else
                                text(0,0.57,'----------------');
                                text(0,0.52,'Plot false origin (UTM)');
                            end
                            text(0,0.45,sprintf('E: %.3f m', EAST_O));
                            text(0,0.40,sprintf('N: %.3f m', NORTH_O));
                        end
                        
                        %satellite number
                        f4 = subplot(7,3,[13 14 15]);
                        plot(nsat); grid on;
                        title('Number of satellites');
                        
                        if (~state.isModeBlock())
                            %EAST plot
                            f5 = subplot(7,3,[16 17 18]);
                            plot(EAST-EAST_O); grid on
                            hold on
                            pos = find(pivot_OUT == 0);
                            if (~isempty(pos))
                                plot(pos, EAST(pivot_OUT == 0)-EAST_O,'.y');
                            end
                            if (o1 == 1) && ((mode == goGNSS.MODE_PP_KF_C_SA) || (mode == goGNSS.MODE_PP_KF_CP_SA) || (mode == goGNSS.MODE_PP_KF_C_DD) || (mode == goGNSS.MODE_PP_KF_CP_DD) || (mode == goGNSS.MODE_RT_NAV))
                                plot([1, nSol], [EAST(end)-EAST_O EAST(end)-EAST_O],'r');
                                title('East coordinates (blue); Not processed / dynamics only (yellow); Final positioning (red)');
                            else
                                title('East coordinates (blue); Not processed / dynamics only (yellow)');
                            end
                            
                            %NORTH plot
                            f6 = subplot(7,3,[19 20 21]);
                            plot(NORTH-NORTH_O); grid on
                            hold on
                            pos = find(pivot_OUT == 0);
                            if (~isempty(pos))
                                plot(pos, NORTH(pivot_OUT == 0)-NORTH_O,'.y');
                            end
                            if (o1 == 1) && ((mode == goGNSS.MODE_PP_KF_C_SA) || (mode == goGNSS.MODE_PP_KF_CP_SA) || (mode == goGNSS.MODE_PP_KF_C_DD) || (mode == goGNSS.MODE_PP_KF_CP_DD) || (mode == goGNSS.MODE_RT_NAV))
                                plot([1, nSol], [NORTH(end)-NORTH_O NORTH(end)-NORTH_O],'r');
                                title('North coordinates (blue); Not processed / dynamics only (yellow); Final positioning (red)');
                            else
                                title('North coordinates (blue); Not processed / dynamics only (yellow)');
                            end
                        end
                        
                        %print PDF
                        print(f, '-dpdf', [filerootOUT '_report']);
                        
                        %remove figure
                        close(f)
                    end
                end
                
                %----------------------------------------------------------------------------------------------
                %% TROPOSPHERE FILE SAVING
                %----------------------------------------------------------------------------------------------
                az_el=reshape(serialize([az_r(:) el_r(:)]'),size(az_r,1)*2,size(az_r,2));
                if (goGNSS.isPP(mode) && state.isModeSA() && state.isModePh() && state.isModeKM() && flag_tropo && (~isempty(EAST)))
                    
                  %IF USER WANTS TO OUTPUT OTHER REPORTS SUCH AS ESTIMATED TROPO DELAYS 
                  if getappdata(0,'report_option') == 1 
                      
                    %display information
                    fprintf('Writing troposphere file...\n');
                    %file saving
                    if (strcmp(fsep_char,'default'))
                        head_str = '    Date          GPS time         GPS week          GPS tow           ZHD[m]           ZTD[m]          TGN[mm]          TGE[mm]           ZWD[m]          PWV[mm]';
                        row_str = '%02d/%02d/%02d    %02d:%02d:%06.3f %16d %16.3f %16.5f %16.5f %16.5f %16.5f %16.5f %16.5f';
                        for s = 1 : nSatTot
                            head_str = [head_str '     ' constellations.systems(s) num2str(constellations.PRN(s),'%02d') '-az[deg]' ]; %#ok<AGROW>
                            head_str = [head_str '      ' constellations.systems(s) num2str(constellations.PRN(s),'%02d') '-el[deg]' ]; %#ok<AGROW>
                            row_str  = [row_str  '%16.5f %16.5f']; %#ok<AGROW>
                        end
                        for s = 1 : nSatTot
                            head_str = [head_str '      STD[m] ' constellations.systems(s) num2str(constellations.PRN(s),'%02d')]; %#ok<AGROW>
                            row_str  = [row_str  '%16.5f']; %#ok<AGROW>
                        end
                        head_str = [head_str '\n']; %#ok<AGROW>
                        row_str  = [row_str  '\n']; %#ok<AGROW>
                    else
                        head_str = strcat('Date',fsep_char,'GPS time',fsep_char,'GPS week',fsep_char,'GPS tow',fsep_char,'ZHD[m]',fsep_char,'ZTD[m]',fsep_char,'TGN[mm]',fsep_char,'TGE[mm]',fsep_char,'ZWD[m]',fsep_char,'PWV[mm]');
                        row_str = strcat('%02d/%02d/%02d',fsep_char,'%02d:%02d:%f',fsep_char,'%d',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f',fsep_char,'%f');
                        for s = 1 : nSatTot
                            head_str = [head_str,fsep_char,'az[d]',fsep_char,'el[d]' constellations.systems(s) num2str(constellations.PRN(s),'%02d')]; %#ok<AGROW>
                            row_str  = [row_str  '%16.5f %16.5f']; %#ok<AGROW>
                        end
                        for s = 1 : nSatTot
                            head_str = [head_str,fsep_char,'STD[m] ',constellations.systems(s),num2str(constellations.PRN(s),'%02d')]; %#ok<AGROW>
                            row_str  = [row_str, fsep_char,'%16.5f']; %#ok<AGROW>
                        end
                        head_str = [head_str '\n']; %#ok<AGROW>
                        row_str  = [row_str  '\n']; %#ok<AGROW>
                    end
                    fid_tropo = fopen([filerootOUT '_tropo.txt'], 'wt');
                    fprintf(fid_tropo, head_str);
                    for i = 1 : nSol
                        
                        %***********FILE WRITING
                        if getappdata(0,'pwv_state') == 1 %WHEN PWV CHECKBOX IS CHECKED
                           fprintf(fid_tropo, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), tow(i), ZHD(i), estim_tropo(i), estim_gradN(i)*1e3, estim_gradE(i)*1e3, ZWD(i), PWV(i), az_el(:,i)', STDs(:,i)');
                           
                        else
                            fprintf(fid_tropo, row_str, date_R(i,1), date_R(i,2), date_R(i,3), date_R(i,4), date_R(i,5), date_R(i,6), week_R(i), tow(i), ZHD(i), estim_tropo(i), estim_gradN(i)*1e3, estim_gradE(i)*1e3, ZWD(i), az_el(:,i)', STDs(:,i)');
                           
                        end
                    end
                    
                    %file closing
                    fclose(fid_tropo);
                    
                  end %//if getappdata(0,'report_option') == 1 
                  
                end %//if (goGNSS.isPP(mode) && state.isModeSA() && state.isModePh() && state.isModeKM() && flag_tropo && (~isempty(EAST)))
                    
                
            end %//if goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)
            
            %----------------------------------------------------------------------------------------------
            %% NMEA FILE SAVING
            %----------------------------------------------------------------------------------------------
            
            %if any positioning was done (either post-processing or real-time)
            if (goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)) && (exist('is_batch','var') && is_batch == 0) && state.isOutNMEA()
                %display information
                fprintf('Writing NMEA file...\n');
                %file saving
                fid_nmea = fopen([filerootOUT '_NMEA.txt'], 'wt');
                
                for i = 1 : nSol / (1 + state.isForwardBackwardKF())
                    id = (state.getForwardBackwardKF() < 0) * (nSol/2) + (state.getForwardBackwardKF() > 0) * (nSol + 1 - 2 * i) + i;
                    
                    %active satellites
                    sat = find(abs(conf_sat_OUT(:,i)));
                    %number of active satellites
                    nsat = length(sat);
                    %visible satellites
                    vsat = find(elR(:,i) > 0);
                    
                    %NMEA string generation
                    GGAstring = NMEA_GGA_gen(pos_KAL(:,id), nsat, date_R(i,:), HDOP(id), mode);
                    if (pivot_OUT(id) ~= 0)
                        RMCstring = NMEA_RMC_gen(pos_KAL(:,id), date_R(i,:));
                        GSVstring = NMEA_GSV_gen(vsat, elR(vsat,id), azR(vsat,id), snr_R(vsat,i), constellations);
                        GSAstring = NMEA_GSA_gen(sat, PDOP(id), HDOP(id), VDOP(id), 'M', '3');
                        if (mode_vinc == 0) && ((mode == goGNSS.MODE_PP_KF_C_SA) || (mode == goGNSS.MODE_PP_KF_CP_SA) || (mode == goGNSS.MODE_PP_KF_C_DD) || (mode == goGNSS.MODE_PP_KF_CP_DD) || (mode == goGNSS.MODE_RT_NAV))
                            PGGPKstring = NMEA_PGGPK_gen(sat, KPDOP(id), KHDOP(id), KVDOP(id), 'S');
                        end
                    else
                        GSAstring = NMEA_GSA_gen(sat, PDOP(id), HDOP(id), VDOP(id), 'M', '1');
                        if (mode_vinc == 0) && ((mode == goGNSS.MODE_PP_KF_C_SA) || (mode == goGNSS.MODE_PP_KF_CP_SA) || (mode == goGNSS.MODE_PP_KF_C_DD) || (mode == goGNSS.MODE_PP_KF_CP_DD) || (mode == goGNSS.MODE_RT_NAV))
                            PGGPKstring = NMEA_PGGPK_gen(sat, KPDOP(id), KHDOP(id), KVDOP(id), 'D');
                        end
                    end
                    
                    %NMEA file write
                    if (pivot_OUT(id) ~= 0)
                        fprintf(fid_nmea, [RMCstring '\n']);
                        fprintf(fid_nmea, [GSVstring '\n']);
                    end
                    fprintf(fid_nmea, [GGAstring '\n']);
                    fprintf(fid_nmea, [GSAstring '\n']);
                    if (mode_vinc == 0) && ((mode == goGNSS.MODE_PP_KF_C_SA) || (mode == goGNSS.MODE_PP_KF_CP_SA) || (mode == goGNSS.MODE_PP_KF_C_DD) || (mode == goGNSS.MODE_PP_KF_CP_DD) || (mode == goGNSS.MODE_RT_NAV))
                        fprintf(fid_nmea, [PGGPKstring '\n']);
                    end
                end
                fclose(fid_nmea);
            end
            
            %----------------------------------------------------------------------------------------------
            %% GOOGLE EARTH FILE SAVING (KML FILE)
            %----------------------------------------------------------------------------------------------
            
            %if any positioning was done (either post-processing or real-time)
            if (goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)) && state.isOutKML()
                %display information
                fprintf('Writing KML file...\n');
                %"clampToGround" plots the points attached to the ground
                %"absolute" uses the height defined in the tag <coordinates>;
                %N.B. Google Earth uses orthometric heights
                z_pos = 'clampToGround';
                %z_pos = 'absolute';
                %URL to load the icon for the points
                iconR = 'http://maps.google.com/mapfiles/kml/pal2/icon26.png';
                iconM = 'http://maps.google.com/mapfiles/kml/shapes/square.png';
                iconP = 'http://maps.google.com/mapfiles/kml/shapes/square.png';
                good_point_colorR = 'fff5005a';
                bad_point_colorR = 'ff0000ff';
                dyn_point_colorR = 'ff00ffff';
                point_colorM = 'ff00ffff';
                point_colorP = 'ff32cd32';
                %point size
                scaleR = 0.2;
                scaleM = 0.8;
                scaleP = 0.8;
                %line color and thickness (rover track)
                line_colorR = 'fff5005a';
                line_widthR = 1;
                %line color and thickness (stop-go-stop direction)
                line_colorG = 'ff0000ff';
                line_widthG = 4;
                %label color
                label_colorM = point_colorM;
                label_colorP = point_colorP;
                %label size
                label_scaleM = 0.7;
                label_scaleP = 0.7;
                %initialization
                
                phiM = zeros(1, nSol / (1 + state.isForwardBackwardKF()));
                lamM = zeros(1, nSol / (1 + state.isForwardBackwardKF()));
                hM = zeros(1, nSol / (1 + state.isForwardBackwardKF()));
                
                %threshold on KHDOP
                if (o1 == 1)
                    KHDOP_thres = median(KHDOP);
                else
                    KHDOP_thres = 2;
                end
                
                %if relative positioning (i.e. with master station)
                if (goGNSS.isDD(mode) || mode == goGNSS.MODE_RT_NAV)
                    %master station coordinates
                    for i = 1 : nSol / (1 + state.isForwardBackwardKF())
                        if (sum(abs(pos_M(:,i))) ~= 0)
                            XM = pos_M(1,i);
                            YM = pos_M(2,i);
                            ZM = pos_M(3,i);
                            
                            %conversion from cartesian to geodetic coordinates
                            [phiM(i), lamM(i), hM(i)] = cart2geod(XM, YM, ZM);
                            
                            %conversion from radians to degrees
                            lamM(i) = lamM(i)*180/pi;
                            phiM(i) = phiM(i)*180/pi;
                        else
                            lamM(i) = 0;
                            phiM(i) = 0;
                            hM(i) = 0;
                        end
                    end
                end
                
                pos = find(filerootOUT == '/');
                if (isempty(pos))
                    pos = find(filerootOUT == '\');
                end
                kml_name = checkPath(filerootOUT(pos(end)+1:end));
                
                %file saving (Google Earth KML)
                fid_kml = fopen([filerootOUT '.kml'], 'wt');
                fprintf(fid_kml, '<?xml version="1.0" encoding="UTF-8"?>\n');
                fprintf(fid_kml, '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n');
                fprintf(fid_kml, '<Document>\n');
                fprintf(fid_kml, '\t<name>%s</name>\n', [kml_name '.kml']);
                fprintf(fid_kml, '\t<snippet>created by goGPS</snippet>\n');
                fprintf(fid_kml, '\t\t<Style id="go1">\n');
                fprintf(fid_kml, '\t\t\t<IconStyle>\n');
                fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',good_point_colorR);
                fprintf(fid_kml, '\t\t\t\t<scale>%.2f</scale>\n',scaleR);
                fprintf(fid_kml, '\t\t\t\t<Icon>\n');
                fprintf(fid_kml, '\t\t\t\t\t<href>%s</href>\n',iconR);
                fprintf(fid_kml, '\t\t\t\t</Icon>\n');
                fprintf(fid_kml, '\t\t\t</IconStyle>\n');
                fprintf(fid_kml, '\t\t</Style>\n');
                fprintf(fid_kml, '\t\t<Style id="go2">\n');
                fprintf(fid_kml, '\t\t\t<IconStyle>\n');
                fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',bad_point_colorR);
                fprintf(fid_kml, '\t\t\t\t<scale>%.2f</scale>\n',scaleR);
                fprintf(fid_kml, '\t\t\t\t<Icon>\n');
                fprintf(fid_kml, '\t\t\t\t\t<href>%s</href>\n',iconR);
                fprintf(fid_kml, '\t\t\t\t</Icon>\n');
                fprintf(fid_kml, '\t\t\t</IconStyle>\n');
                fprintf(fid_kml, '\t\t</Style>\n');
                fprintf(fid_kml, '\t\t<Style id="go3">\n');
                fprintf(fid_kml, '\t\t\t<IconStyle>\n');
                fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',dyn_point_colorR);
                fprintf(fid_kml, '\t\t\t\t<scale>%.2f</scale>\n',scaleR);
                fprintf(fid_kml, '\t\t\t\t<Icon>\n');
                fprintf(fid_kml, '\t\t\t\t\t<href>%s</href>\n',iconR);
                fprintf(fid_kml, '\t\t\t\t</Icon>\n');
                fprintf(fid_kml, '\t\t\t</IconStyle>\n');
                fprintf(fid_kml, '\t\t</Style>\n');
                fprintf(fid_kml, '\t\t<Style id="master">\n');
                fprintf(fid_kml, '\t\t\t<IconStyle>\n');
                fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',point_colorM);
                fprintf(fid_kml, '\t\t\t\t<colorMode>normal</colorMode>\n');
                fprintf(fid_kml, '\t\t\t\t<scale>%.2f</scale>\n',scaleM);
                fprintf(fid_kml, '\t\t\t\t<Icon>\n');
                fprintf(fid_kml, '\t\t\t\t\t<href>%s</href>\n',iconM);
                fprintf(fid_kml, '\t\t\t\t</Icon>\n');
                fprintf(fid_kml, '\t\t\t</IconStyle>\n');
                fprintf(fid_kml, '\t\t\t<LabelStyle>\n');
                fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',label_colorM);
                fprintf(fid_kml, '\t\t\t\t<scale>%s</scale>\n',label_scaleM);
                fprintf(fid_kml, '\t\t\t</LabelStyle>\n');
                fprintf(fid_kml, '\t\t</Style>\n');
                fprintf(fid_kml, '\t\t<Style id="ppos">\n');
                fprintf(fid_kml, '\t\t\t<IconStyle>\n');
                fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',point_colorP);
                fprintf(fid_kml, '\t\t\t\t<colorMode>normal</colorMode>\n');
                fprintf(fid_kml, '\t\t\t\t<scale>%.2f</scale>\n',scaleP);
                fprintf(fid_kml, '\t\t\t\t<Icon>\n');
                fprintf(fid_kml, '\t\t\t\t\t<href>%s</href>\n',iconP);
                fprintf(fid_kml, '\t\t\t\t</Icon>\n');
                fprintf(fid_kml, '\t\t\t</IconStyle>\n');
                fprintf(fid_kml, '\t\t\t<LabelStyle>\n');
                fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',label_colorP);
                fprintf(fid_kml, '\t\t\t\t<scale>%s</scale>\n',label_scaleP);
                fprintf(fid_kml, '\t\t\t</LabelStyle>\n');
                fprintf(fid_kml, '\t\t</Style>\n');
                fprintf(fid_kml, '\t\t<Style id="goLine1">\n');
                fprintf(fid_kml, '\t\t\t<LineStyle>\n');
                fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',line_colorR);
                fprintf(fid_kml, '\t\t\t\t<width>%d</width>\n',line_widthR);
                fprintf(fid_kml, '\t\t\t</LineStyle>\n');
                fprintf(fid_kml, '\t\t</Style>\n');
                if (flag_stopGOstop && goGNSS.isPP(mode)) %stop-go-stop and post-processing
                    fprintf(fid_kml, '\t\t<Style id="goLine2">\n');
                    fprintf(fid_kml, '\t\t\t<LineStyle>\n');
                    fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',line_colorG);
                    fprintf(fid_kml, '\t\t\t\t<width>%d</width>\n',line_widthG);
                    fprintf(fid_kml, '\t\t\t</LineStyle>\n');
                    fprintf(fid_kml, '\t\t</Style>\n');
                end
                if (goGNSS.isDD(mode) || mode == goGNSS.MODE_RT_NAV) %relative positioning
                    for i = 1 : length(phiM)
                        if (lamM(i) ~= 0 || phiM(i) ~= 0 || hM(i) ~= 0)
                            if (i == 1) || (lamM(i)~=lamM(i-1) || phiM(i)~=phiM(i-1) || hM(i)~=hM(i-1))
                                fprintf(fid_kml, '\t\t<Placemark>\n');
                                fprintf(fid_kml, '\t\t\t<name>Master station</name>\n');
                                fprintf(fid_kml, '\t\t\t<description><![CDATA[ <i>Latitude:</i> %.8f &#176;<br/> <i>Longitude:</i> %.8f &#176;<br/> <i>Elevation (ellips.):</i> %.1f m<br/>]]></description>\n',phiM(i),lamM(i),hM(i));
                                fprintf(fid_kml, '\t\t\t<styleUrl>#master</styleUrl>\n');
                                fprintf(fid_kml, '\t\t\t<Point>\n');
                                fprintf(fid_kml, '\t\t\t\t<altitudeMode>%s</altitudeMode>\n',z_pos);
                                fprintf(fid_kml, '\t\t\t\t<coordinates>%.8f,%.8f,%.3f</coordinates>\n',lamM(i),phiM(i),h_ortho(i));
                                fprintf(fid_kml, '\t\t\t</Point>\n');
                                fprintf(fid_kml, '\t\t</Placemark>\n');
                            end
                        end
                    end
                end
                fprintf(fid_kml, '\t\t<Placemark>\n');
                fprintf(fid_kml, '\t\t<name>Rover track</name>\n');
                fprintf(fid_kml, '\t\t\t<styleUrl>#goLine1</styleUrl>\n');
                fprintf(fid_kml, '\t\t\t<LineString>\n');
                fprintf(fid_kml, '\t\t\t\t<coordinates>\n\t\t\t\t\t');
                for i = 1 : nSol / (1 + state.isForwardBackwardKF())
                    id = (state.getForwardBackwardKF() < 0) * (nSol/2) + (state.getForwardBackwardKF() > 0) * (nSol + 1 - 2 * i) + i;
                    fprintf(fid_kml, '%.8f,%.8f,0 ',lam_KAL(id),phi_KAL(id));
                end
                fprintf(fid_kml, '\n\t\t\t\t</coordinates>\n');
                fprintf(fid_kml, '\t\t\t</LineString>\n');
                fprintf(fid_kml, '\t\t</Placemark>\n');
                if (flag_stopGOstop && flag_var_dyn_model && mode == goGNSS.MODE_PP_KF_CP_DD)
                    
                    [P1Lat, P1Lon] = cart2geod(P1_GLB(1), P1_GLB(2), P1_GLB(3));
                    [P2Lat, P2Lon] = cart2geod(P2_GLB(1), P2_GLB(2), P2_GLB(3));
                    
                    fprintf(fid_kml, '\t\t<Placemark>\n');
                    fprintf(fid_kml, '\t\t<name>Estimated direction</name>\n');
                    fprintf(fid_kml, '\t\t\t<styleUrl>#goLine2</styleUrl>\n');
                    fprintf(fid_kml, '\t\t\t<LineString>\n');
                    fprintf(fid_kml, '\t\t\t\t<coordinates>\n\t\t\t\t\t');
                    fprintf(fid_kml, '%.8f,%.8f,0 ',P1Lon*180/pi,P1Lat*180/pi);
                    fprintf(fid_kml, '%.8f,%.8f,0 ',P2Lon*180/pi,P2Lat*180/pi);
                    fprintf(fid_kml, '\n\t\t\t\t</coordinates>\n');
                    fprintf(fid_kml, '\t\t\t</LineString>\n');
                    fprintf(fid_kml, '\t\t</Placemark>\n');
                end
                fprintf(fid_kml, '\t\t<Folder>\n');
                fprintf(fid_kml, '\t\t<name>Rover positioning</name>\n');
                for i = 1 : nSol / (1 + state.isForwardBackwardKF())
                    id = (state.getForwardBackwardKF() < 0) * (nSol/2) + (state.getForwardBackwardKF() > 0) * (nSol + 1 - 2 * i) + i;
                    
                    fprintf(fid_kml, '\t\t<Placemark>\n');
                    if (pivot_OUT(id) == 0)
                        fprintf(fid_kml, '\t\t\t<styleUrl>#go3</styleUrl>\n');
                    elseif ~(mode == goGNSS.MODE_PP_BLK_CP_DD_STATIC) && (KHDOP(id)>KHDOP_thres)
                        fprintf(fid_kml, '\t\t\t<styleUrl>#go2</styleUrl>\n');
                    else
                        fprintf(fid_kml, '\t\t\t<styleUrl>#go1</styleUrl>\n');
                    end
                    fprintf(fid_kml, '\t\t\t<Point>\n');
                    fprintf(fid_kml, '\t\t\t\t<altitudeMode>%s</altitudeMode>\n',z_pos);
                    fprintf(fid_kml, '\t\t\t\t<coordinates>%.8f,%.8f,%.3f</coordinates>\n',lam_KAL(i),phi_KAL(i),h_KAL(i));
                    fprintf(fid_kml, '\t\t\t</Point>\n');
                    fprintf(fid_kml, '\t\t</Placemark>\n');
                end
                fprintf(fid_kml, '\t\t</Folder>\n');
                
                if (mode_vinc == 0) && ((mode == goGNSS.MODE_PP_KF_C_SA) || (mode == goGNSS.MODE_PP_KF_CP_SA) || (mode == goGNSS.MODE_PP_KF_C_DD) || (mode == goGNSS.MODE_PP_KF_CP_DD) || (mode == goGNSS.MODE_RT_NAV) || (mode == goGNSS.MODE_PP_BLK_CP_DD_STATIC))
                    if (o1 == 1) && (nSol ~= 0)
                        %static positioning coordinates
                        phiP = phi_KAL(end);
                        lamP = lam_KAL(end);
                        hP   = h_KAL(end);
                        
                        fprintf(fid_kml, '\t\t<Placemark>\n');
                        fprintf(fid_kml, '\t\t\t<name>Static positioning</name>\n');
                        fprintf(fid_kml, '\t\t\t<description><![CDATA[ <i>Latitude:</i> %.8f &#176;<br/> <i>Longitude:</i> %.8f &#176;<br/> <i>Elevation (ellips.):</i> %.1f m<br/>]]></description>\n',phiP,lamP,hP);
                        fprintf(fid_kml, '\t\t\t<styleUrl>#ppos</styleUrl>\n');
                        fprintf(fid_kml, '\t\t\t<Point>\n');
                        fprintf(fid_kml, '\t\t\t\t<altitudeMode>%s</altitudeMode>\n',z_pos);
                        fprintf(fid_kml, '\t\t\t\t<coordinates>%.8f,%.8f,%.3f</coordinates>\n',lamP,phiP,hP);
                        fprintf(fid_kml, '\t\t\t</Point>\n');
                        fprintf(fid_kml, '\t\t</Placemark>\n');
                    end
                end
                fprintf(fid_kml, '</Document>\n</kml>');
                fclose(fid_kml);
            end
            
            
           %-----------------------------------------------------
           %% ***********************OUTPUT OTHER REPORTS & GRAPHS
           %-----------------------------------------------------
           
           %GET STATE OF OPTION FOR OTHER REPORTS IN goGPS OUTPUTS
           report_option = getappdata(0,'report_option');
% %            any([any([state.isOutTropoPDF(),state.isOutDtPDF(),state.isOutPhaseResPDF(),...
% %                      state.isOutCodeResPDF(),state.isOutReportPDF(),state.isOutKML(),state.isOutNMEA()]),...
% %                 all([state.isOutTropoPDF(),state.isOutDtPDF(),state.isOutPhaseResPDF(), state.isOutCodeResPDF(),...
% %                      state.isOutReportPDF(),state.isOutKML(),state.isOutNMEA()])])
               
          if report_option == 1 %IF USER WANTS TO OUTPUT OTHER REPORTS 
            
            %----------------------------------------------------------------------------------------------
            %% REPRESENTATION OF THE ESTIMATED ERROR COVARIANCE (AND TEXT FILE SAVING)
            %----------------------------------------------------------------------------------------------
            
            %if any positioning was done (either post-processing or real-time, not constrained)
            if ((goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)) && (~isempty(EAST)) && (mode_vinc == 0))
                
                %display information
                fprintf('Writing estimated error covariance files...\n');
                %covariance propagation
                Cee_XYZ = Cee_OUT([1 o1+1 o2+1],[1 o1+1 o2+1],:);
                Cee_ENU = global2localCov(Cee_OUT([1 o1+1 o2+1],[1 o1+1 o2+1],:), Xhat_t_t_OUT([1 o1+1 o2+1],:));
                
                if (flag_cov == 1 && mode_user == 1)
                    %trajectory plotting
                    figure
                    plot(EAST, NORTH, '.r'); axis equal
                    xlabel('EAST [m]'); ylabel('NORTH [m]'); grid on;
                    
                    hold on
                    for i = 1 : size(Cee_ENU,3)         % ellipse definition
                        T = chol(Cee_ENU(1:2,1:2,i));   % Cholesky decomposition
                        n = size(x_circle,1);
                        x_ellipse = zeros(n,2);         % pre-allocation
                        for j = 1 : n                   % ellipse definition
                            x_ellipse(j,:) = x_circle(j,:) * T + [EAST(i), NORTH(i)];
                        end
                        he = plot(x_ellipse(:,1),x_ellipse(:,2),'b');
                        if (any(fixed_amb) && fixed_amb(i) == 0)
                            set(he,'Color',[0.8 0.5 0])
                        end
                    end
                    hold off
                end
                
                %file saving (XYZ covariance)
                fid_cov = fopen([filerootOUT '_cov_XYZ.txt'], 'wt');
                fprintf(fid_cov, '       XX              XY              XZ              YY              YZ              ZZ\n');
                for i = 1 : length(phi_KAL)
                    fprintf(fid_cov, '%.3e% 16.3e% 16.3e% 16.3e% 16.3e% 16.3e\n', Cee_XYZ(1,1,i), Cee_XYZ(1,2,i), ...
                        Cee_XYZ(1,3,i), Cee_XYZ(2,2,i), Cee_XYZ(2,3,i), Cee_XYZ(3,3,i));
                end
                fclose(fid_cov);
                
                %file saving (ENU covariance)
                fid_cov = fopen([filerootOUT '_cov_ENU.txt'], 'wt');
                fprintf(fid_cov, ' EastEast       EastNorth          EastUp      NorthNorth         NorthUp            UpUp\n');
                for i = 1 : length(phi_KAL)
                    fprintf(fid_cov, '%.3e% 16.3e% 16.3e% 16.3e% 16.3e% 16.3e\n', Cee_ENU(1,1,i), Cee_ENU(1,2,i), ...
                        Cee_ENU(1,3,i), Cee_ENU(2,2,i), Cee_ENU(2,3,i), Cee_ENU(3,3,i));
                end
                fclose(fid_cov);
            end
            
            %----------------------------------------------------------------------------------------------
            %% REPRESENTATION OF THE ESTIMATED COORDINATES TIME SERIES
            %----------------------------------------------------------------------------------------------
            
            if (mode ~= goGNSS.MODE_PP_LS_CP_VEL) && (goGNSS.isPP(mode) && (~isempty(time_GPS)))  && (mode_user == 1) && ~state.isModeBlock()
                epochs = (time_GPS-time_GPS(1))/interval;
                
                time = GPS_Time(GPS_Time.GPS_ZERO, time_GPS);
                id = 1 : numel(epochs);
                switch state.getForwardBackwardKF()
                    case -1
                        epochs = [flipud(epochs); epochs]; %#ok<AGROW>
                        id = (numel(epochs)/2 +1) : numel(epochs);
                    case  1
                        epochs = [epochs; flipud(epochs)]; %#ok<AGROW>
                        time = GPS_Time(GPS_Time.GPS_ZERO, flipud(time_GPS));
                        id = (numel(epochs)/2 +1) : numel(epochs);
                end
                figure;
                ax = zeros(3,1);
                ax(1) = subplot(3,1,1); hold on; grid on
                ax(2) = subplot(3,1,2); hold on; grid on
                ax(3) = subplot(3,1,3); hold on; grid on
                xlabel('epoch')
                
                %if relative positioning (i.e. with master station)
                if (goGNSS.isDD(mode) || mode == goGNSS.MODE_RT_NAV)
                    subplot(ax(1))
                    plot(epochs, EAST,'.'); title('EAST'); ylabel('[m]')
                    subplot(ax(2))
                    plot(epochs, NORTH,'.'); title('NORTH'); ylabel('[m]')
                    subplot(ax(3))
                    plot(epochs, UP_KAL,'.'); title('UP'); ylabel('[m]')
                    
                    pos = find(fixed_amb == 1);
                    plot(ax(1), epochs(pos), EAST(pos),'xr')
                    plot(ax(2), epochs(pos), NORTH(pos),'xr')
                    plot(ax(3), epochs(pos), UP_KAL(pos),'xr')
                    
                    pos = find(pivot_OUT == 0);
                    plot(ax(1), epochs(pos), EAST(pos),'.y')
                    plot(ax(2), epochs(pos), NORTH(pos),'.y')
                    plot(ax(3), epochs(pos), UP_KAL(pos),'.y')
                    
                    Fig_Lab.plotENU(time,[EAST(id) NORTH(id) UP_KAL(id)]);
                else
                    subplot(ax(1))
                    plot(epochs, EAST_UTM,'.'); title('EAST UTM'); ylabel('[m]')
                    subplot(ax(2))
                    plot(epochs, NORTH_UTM,'.'); title('NORTH UTM'); ylabel('[m]')
                    subplot(ax(3))
                    plot(epochs, h_KAL,'.'); title('ellipsoidal height'); ylabel('[m]')
                    
                    pos = find(fixed_amb == 1);
                    plot(ax(1), epochs(pos), EAST_UTM(pos),'xr')
                    plot(ax(2), epochs(pos), NORTH_UTM(pos),'xr')
                    plot(ax(3), epochs(pos), h_KAL(pos),'xr')
                    
                    pos = find(pivot_OUT == 0);
                    plot(ax(1), epochs(pos), EAST_UTM(pos),'.y')
                    plot(ax(2), epochs(pos), NORTH_UTM(pos),'.y')
                    plot(ax(3), epochs(pos), h_KAL(pos),'.y')
                    
                    Fig_Lab.plotENU(time,[EAST(id) NORTH(id) h_KAL(id)]);
                end
                linkaxes(ax,'x');
            end
            
            %----------------------------------------------------------------------------------------------
            %% REPRESENTATION OF THE REFERENCE TRAJECTORY
            %----------------------------------------------------------------------------------------------
            
            % %if any positioning was done (either post-processing or real-time, but constrained)
            % if ((goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)) && mode_ref == 1)
            %     [EAST_ref, NORTH_ref, h_ref] = cart2plan(ref_path(:,1), ref_path(:,2),ref_path(:,3));
            %
            %     %reference data plot
            %     plot(EAST_ref, NORTH_ref, 'm', 'LineWidth', 2)
            %
            %     hold off
            % end
            
            %----------------------------------------------------------------------------------------------
            %% REPRESENTATION OF THE 2D TRAJECTORY
            %----------------------------------------------------------------------------------------------
            
            %if any positioning was done (either post-processing or real-time, not constrained)
            if ((goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)) && (~isempty(EAST)) && mode_user == 1) && ~state.isModeBlock()
                %2D plot
                figure
                plot(EAST, NORTH, '.r');
                xlabel('EAST [m]'); ylabel('NORTH [m]'); grid on
            end
            
            %----------------------------------------------------------------------------------------------
            %% REPRESENTATION OF THE 3D TRAJECTORY
            %----------------------------------------------------------------------------------------------
            
            %if any positioning was done (either post-processing or real-time, not constrained)
            if ((goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)) && (~isempty(EAST)) && mode_user == 1) && ~state.isModeBlock()
                
                %     %3D plot (XYZ)
                %     figure
                %     plot3(X_KAL, Y_KAL, Z_KAL, '.r');
                %     xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]'); grid on
                
                %3D plot (ENU)
                figure
                plot3(EAST, NORTH, h_KAL, '.r');
                xlabel('EAST [m]'); ylabel('NORTH [m]'); zlabel('h [m]'); grid on
            end
            
            %----------------------------------------------------------------------------------------------
            %% REPRESENTATION OF RESIDUALS MEAN AND STANDARD DEVIATION
            %----------------------------------------------------------------------------------------------
            
            if (goGNSS.isPP(mode) || (mode == goGNSS.MODE_RT_NAV)) && (~isempty(EAST))
                %code
                
                switch state.getForwardBackwardKF()
                    case -1, id = (size(RES_CODE1,2)/2 +1) : size(RES_CODE1,2);
                    case  1, id = size(RES_CODE1,2) : -1 : (size(RES_CODE1,2)/2 +1);
                    case  0, id = 1 : size(RES_CODE1,2);
                end
                
                if (any(RES_CODE1(:)))
                    RES_CODE1_mean = nan(size(RES_CODE1,1),1);
                    RES_CODE1_stdv = nan(size(RES_CODE1,1),1);
                    for s = 1 : nSatTot
                        row = RES_CODE1(s,id);
                        row(outliers_CODE1(s,id) == 1) = [];
                        row(isnan(row(1,:))) = [];
                        RES_CODE1_mean(s,1) = mean(row);
                        RES_CODE1_stdv(s,1) = std(row);
                    end
                    figure;
                    errorbar(1:nSatTot, RES_CODE1_mean, RES_CODE1_stdv,'k*');
                    grid on;
                    title('RESIDUALS MEAN AND ST.DEV: CODE 1st FREQ.');
                    xlabel('Satellite ID');
                    ylabel('[m]');
                end
                if (any(RES_CODE2(:)))
                    RES_CODE2_mean = nan(size(RES_CODE2,1),1);
                    RES_CODE2_stdv = nan(size(RES_CODE2,1),1);
                    for s = 1 : nSatTot
                        row = RES_CODE2(s,id);
                        row(outliers_CODE2(s,id) == 1) = [];
                        row(isnan(row(1,:))) = [];
                        RES_CODE2_mean(s,1) = mean(row);
                        RES_CODE2_stdv(s,1) = std(row);
                    end
                    figure;
                    errorbar(1:nSatTot, RES_CODE2_mean, RES_CODE2_stdv,'k*');
                    grid on;
                    title('RESIDUALS MEAN AND ST.DEV: CODE 2nd FREQ.');
                    xlabel('Satellite ID');
                    ylabel('[m]');
                end
                %phase
                if (any(RES_PHASE1(:)))
                    RES_PHASE1_mean = nan(size(RES_PHASE1,1),1);
                    RES_PHASE1_stdv = nan(size(RES_PHASE1,1),1);
                    for s = 1 : nSatTot
                        row = RES_PHASE1(s,id);
                        row(outliers_PHASE1(s,id) == 1) = [];
                        row(isnan(row(1,:))) = [];
                        RES_PHASE1_mean(s,1) = mean(row);
                        RES_PHASE1_stdv(s,1) = std(row);
                    end
                    figure;
                    errorbar(1:nSatTot, RES_PHASE1_mean, RES_PHASE1_stdv,'k*');
                    grid on;
                    title('RESIDUALS MEAN AND ST.DEV: PHASE 1st FREQ.');
                    xlabel('Satellite ID');
                    ylabel('[m]');
                end
                if (any(RES_PHASE2(:)))
                    RES_PHASE2_mean = nan(size(RES_PHASE2,1),1);
                    RES_PHASE2_stdv = nan(size(RES_PHASE2,1),1);
                    for s = 1 : nSatTot
                        row = RES_PHASE2(s,id);
                        row(outliers_PHASE2(s,id) == 1) = [];
                        row(isnan(row(1,:))) = [];
                        RES_PHASE2_mean(s,1) = mean(row);
                        RES_PHASE2_stdv(s,1) = std(row);
                    end
                    figure;
                    errorbar(1:nSatTot, RES_PHASE2_mean, RES_PHASE2_stdv,'k*');
                    grid on;
                    title('RESIDUALS MEAN AND ST.DEV: PHASE 2nd FREQ.');
                    xlabel('Satellite ID');
                    ylabel('[m]');
                end
            end
            
            %----------------------------------------------------------------------------------------------
            %% BUILD FILE_NAME_BASE
            %----------------------------------------------------------------------------------------------
            if (session == 1)
                % Composing the name of the batch output
                [dir_name, file_name, ext] = fileparts(f_trg_rec{1}{1});
                marker_trg = file_name(1:4);
                if state.isModeDD() && ~state.isModeSEID()
                    [dir_name, file_name, ext] = fileparts(f_mst_rec{1}{1});
                    marker_mst = [file_name(1:4) '_'];
                else
                    marker_mst = '';
                end
                
                file_name_base = fnp.dateKeyRep(fnp.checkPath(fullfile(state.getOutDir(), sprintf('%s_%s${YYYY}${DOY}', marker_trg, marker_mst))), state.getSessionStart());
                file_name_base = [fnp.dateKeyRep(sprintf('%s_${YYYY}${DOY}',file_name_base), state.getSessionStop()) filesep];
                if ~(exist(file_name_base, 'dir') == 7)
                    mkdir(file_name_base);
                end
            end
            
            %----------------------------------------------------------------------------------------------
            %% STATISTICS COMPUTATION AND VISUALIZATION
            %----------------------------------------------------------------------------------------------
            
            if goGNSS.isPP(mode) && (mode_vinc == 0) && (~isempty(gs.getReferencePath().path)) && (~isempty(EAST_KAL))
                
                switch state.getForwardBackwardKF()
                    case -1, id = (size(RES_CODE1,2)/2 +1) : size(RES_CODE1,2);
                    case  1, id = size(RES_CODE1,2) : -1 : (size(RES_CODE1,2)/2 +1);
                    case  0, id = 1 : size(RES_CODE1,2);
                end
                
                %coordinate transformation
                ref_path = gs.getReferencePath().path;
                [EAST_REF, NORTH_REF, h_REF] = cart2plan(ref_path(:,1), ref_path(:,2), ref_path(:,3));
                
                ref = [EAST_REF NORTH_REF h_REF];
                
                [dist2D, ~] = ref_2d_projection(ref,EAST_UTM(id),NORTH_UTM(id));
                
                fprintf('\n');
                fprintf('-------- STATISTICS ------------');
                fprintf('\n');
                fprintf('Mean2D: %7.4f m\n',mean(dist2D));
                fprintf('Std2D:  %7.4f m\n',std(dist2D,1));
                fprintf('RMS2D:  %7.4f m\n\n',sqrt(std(dist2D)^2+mean(dist2D)^2));
                
                [dist3D,proj] = ref_3d_projection(ref,EAST_UTM,NORTH_UTM,h_KAL);
                
                fprintf('Mean3D: %7.4f m\n',mean(dist3D));
                fprintf('Std3D:  %7.4f m\n',std(dist3D,1));
                fprintf('RMS3D:  %7.4f m\n',sqrt(std(dist3D)^2+mean(dist3D)^2));
                fprintf('--------------------------------\n\n');
            end
            
            %----------------------------------------------------------------------------------------------
            
            if (exist('time_GPS', 'var') && isempty(time_GPS))
                fprintf('\n');
                fprintf('... WARNING: no epochs were available/usable for processing. The result files will be empty.\n');
                fprintf('\n');
            end
            
            %----------------------------------------------------------------------------------------------
            % write report
            report_generator(report);
            %----------------------------------------------------------------------------------------------
            
            if (exist('fout_report','var')), fclose(fout_report); end
           
            if (mode_user && ~is_batch)
                % close all the opened files
                fclose('all');
            end
            
          end %//if report_option == 1
          
            %re-enable MATLAB warnings
            warning on
            
            %evaluate computation time
            
            if is_batch
                nSol = size(date_R,1);
                if (~state.isModeBlock())
                    switch state.getForwardBackwardKF()
                        case -1
                            id_time = nSol;
                            id_data = nSol*2;
                        case  1
                            id_time = 1;
                            id_data = nSol*2;
                        case  0
                            id_time = nSol;
                            id_data = nSol;
                    end
                else
                    id_time = nSol;
                    id_data = 1;
                end
                
                tropo_vec_ZTD = nan(1,86400/interval);
                tropo_vec_ZWD = nan(1,86400/interval);
                if ~flag_init_out
                    flag_init_out = true;
                    fid_extract = fopen(sprintf('%sextraction.txt', file_name_base),'w');
                    
                    if (state.isModeBlock())
                        s_rate = state.getSolutionRate();
                        for s = 1 : numel(s_rate)
                            fid_extract_hr(s) = fopen(sprintf('%sextraction_hr_at%d.txt', file_name_base, s_rate(s)),'w');
                        end
                        fid_extract_float = fopen(sprintf('%sextraction_float.txt', file_name_base),'w');
                        fid_extract_fix = fopen(sprintf('%sextraction_fix.txt', file_name_base),'w');
                    end
                    
                    fid_extract_ZTD = fopen(sprintf('%sZTD.txt', file_name_base),'w');
                    fid_extract_ZWD = fopen(sprintf('%sZWD.txt', file_name_base),'w');
                    
                    fid_extract_POS = fopen(sprintf('%sposition.txt', file_name_base),'w');
                    fprintf(fid_extract_POS,' yyyy-ddd   date          time           UTM east         UTM north      ellips. height        ZTD\n');
                    fprintf(fid_extract_POS,'+--------+----------+----------------+----------------+----------------+----------------+----------------+\n');
                    
                    fid_extract_OBS = fopen(sprintf('%squalityOBS.txt', file_name_base),'w');
                    fprintf(fid_extract_OBS,' yy-ddd  Rover observation file            Rate  #Sat   #Epoch    #Frq   #C1/P1  #C2/P2     #L1     #L2   #DOP1   #DOP2  %%Epoch %%L2/L1    Master observation file            Rate  #Sat   #Epoch    #Frq   #C1/P1  #C2/P2     #L1     #L2   #DOP1   #DOP2  %%Epoch %%L2/L1\n');
                    fprintf(fid_extract_OBS,'+------+------------------------------+--------+-----+--------+-------+--------+-------+-------+-------+-------+-------+-------+------+---------------------------------+--------+-----+--------+-------+--------+-------+-------+-------+-------+-------+-------+------+\n');
                end
                
                % Take the central value of the RINEX as reference time for the solution
                if (state.isModeBlock())
                    date_R = datevec(GPS_Time.GPS_ZERO + (zero_time + go_block.getTimeHR(1e100))/86400);
                    id_time = 1;
                end
                if exist('estim_tropo','var') && (state.isTropoOn())
                    
                    tropo_vec_ZTD(1,1:length(estim_tropo)) = estim_tropo;
                    fprintf(fid_extract_ZTD,'%.6f ', tropo_vec_ZTD);
                    fprintf(fid_extract_ZTD,'\n');
                    if (~isempty(ZHD))
                        tropo_vec_ZWD(1,1:length(estim_tropo)) = estim_tropo-ZHD';
                        fprintf(fid_extract_ZWD,'%.6f ', tropo_vec_ZWD);
                        fprintf(fid_extract_ZWD,'\n');
                    end
                end
                
                if exist('X_KAL','var')
                    fprintf(fid_extract,'%s  %02d/%02d/%02d    %02d:%02d:%06.3f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f\n', fnp.dateKeyRep('${YYYY}-${DOY}',cur_date_start), date_R(id_time,1), date_R(id_time,2), date_R(id_time,3), date_R(id_time,4), date_R(id_time,5), date_R(id_time,6), X_KAL(id_data), Y_KAL(id_data), Z_KAL(id_data), EAST_UTM(id_data), NORTH_UTM(id_data), h_KAL(id_data));
                    nSol = size(Xhat_t_t_OUT,2);
                    for e = 1 : nSol / (1 + state.isForwardBackwardKF())
                        id = (state.getForwardBackwardKF() > 0) * (nSol + 1 - 2 * e) + (state.getForwardBackwardKF() < 0) * (nSol/2) + e;
                        fprintf(fid_extract_POS,' %s  %02d/%02d/%02d    %02d:%02d:%06.3f %16.6f %16.6f %16.6f %15.6f\n', fnp.dateKeyRep('${YYYY}-${DOY}',cur_date_start), date_R(e,1), date_R(e,2), date_R(e,3), date_R(e,4), date_R(e,5), date_R(e,6), EAST_UTM(id), NORTH_UTM(id), h_KAL(id), Xhat_t_t_OUT(end-nC-2,id));
                    end
                    delete([filerootOUT '_*.bin']);
                end
                
                if (state.isModeBlock())
                    s_rate = state.getSolutionRate();
                    for s = 1 : numel(s_rate)
                        if (s > 1)
                            go_block.solveHighRate(s_rate(s));
                        end
                        pos = go_block.getPosHR()';
                        ko_pos = isnan(sum(pos));
                        pos(:,ko_pos) = [];
                        tmp = (GPS_Time.GPS_ZERO + (zero_time + go_block.getTimeHR())/86400);
                        date_s = datevec(tmp(~ko_pos));
                        
                        if ~isempty(pos)
                            [~, ~, h_BLK] = cart2geod(pos(1, :), pos(2, :), pos(3, :));
                            
                            %coordinate transformation (UTM)
                            [EAST_UTM, NORTH_UTM] = cart2plan(pos(1, :)', pos(2, :)', pos(3, :)');
                            fprintf(fid_extract_hr(s), '%c%c%c%c%c%c%c%c  %02d/%02d/%02d    %02d:%02d:%06.3f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f\n', ...
                                [double(repmat(fnp.dateKeyRep('${YYYY}-${DOY}',cur_date_start), size(pos,2), 1)), ...
                                date_s(:,1), date_s(:,2), date_s(:,3), date_s(:,4), date_s(:,5), date_s(:,6), ...
                                pos(1, :)', pos(2, :)', pos(3, :)', EAST_UTM(:), NORTH_UTM(:), h_BLK(:)]');
                        end
                    end
                    
                    pos = go_block.getFloatPos()';
                    [~, ~, h_BLK] = cart2geod(pos(1, :), pos(2, :), pos(3, :));
                    
                    %coordinate transformation (UTM)
                    [EAST_UTM, NORTH_UTM] = cart2plan(pos(1, :), pos(2, :), pos(3, :));
                    fprintf(fid_extract_float,'%s  %02d/%02d/%02d    %02d:%02d:%06.3f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f\n', fnp.dateKeyRep('${YYYY}-${DOY}',cur_date_start), date_R(id_time,1), date_R(id_time,2), date_R(id_time,3), date_R(id_time,4), date_R(id_time,5), date_R(id_time,6), pos(1, id_data), pos(2, id_data), pos(3, id_data), EAST_UTM(id_data), NORTH_UTM(id_data), h_BLK(id_data));
                    
                    pos = go_block.getFixPos()';
                    
                    [~, ~, h_BLK] = cart2geod(pos(1, :), pos(2, :), pos(3, :));
                    
                    %coordinate transformation (UTM)
                    [EAST_UTM, NORTH_UTM] = cart2plan(pos(1, :), pos(2, :), pos(3, :));
                    fprintf(fid_extract_fix,'%s  %02d/%02d/%02d    %02d:%02d:%06.3f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f\n', fnp.dateKeyRep('${YYYY}-${DOY}',cur_date_start), date_R(id_time,1), date_R(id_time,2), date_R(id_time,3), date_R(id_time,4), date_R(id_time,5), date_R(id_time,6), pos(1, id_data), pos(2, id_data), pos(3, id_data), EAST_UTM(id_data), NORTH_UTM(id_data), h_BLK(id_data));
                end
                
                % append report information in the database
                fid_rep_i=fopen([filerootOUT,'_report.txt'],'rt');
                if fid_rep_i~=-1
                    line=fgetl(fid_rep_i);
                    while (isempty(strfind(line,'Observations (RAW)              Start time'))) %#ok<STREMP>
                        line = fgetl(fid_rep_i);
                    end
                    line=fgetl(fid_rep_i);
                    line1=fgetl(fid_rep_i);
                    line1(33:80)='';
                    line2=fgetl(fid_rep_i);
                    if ~isempty(line2)
                        line2(33:80)='';
                    end
                    
                    fprintf(fid_extract_OBS,' %s  %s     %s\n',fnp.dateKeyRep('${YY}-${DOY}',cur_date_start), line1, line2);
                    fclose(fid_rep_i);
                end
                close all
            end
        end
    end
    log.newLine();
    toc
    
    state.setMode(initial_mode);
    mode = initial_mode;
end

if flag_init_out && is_batch && ~state.isModeSEID()
    fclose(fid_extract);
    if (state.isModeBlock())
        for i = 1 : numel(fid_extract_hr)
            fclose(fid_extract_hr(i));
        end
        fclose(fid_extract_float);
        fclose(fid_extract_fix);
    end
    
    fclose(fid_extract_ZTD);
    fclose(fid_extract_ZWD);
    fclose(fid_extract_POS);
    fclose(fid_extract_OBS);
end
end
