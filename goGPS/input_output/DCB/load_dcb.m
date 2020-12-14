function [DCB] = load_dcb(data_dir_dcb, gps_week, time_R, codeC1_R, constellations)

% SYNTAX:
%   [DCB] = load_dcb(data_dir_dcb, gps_week, time_R, codeC1_R, constellations);
%
% INPUT:
%   data_dir_dcb = path to the directory containing DCB files [string]
%   gps_week = reference vector of GPS week numbers
%   time_R = reference vector of GPS time
%   codeC1_R = flag to indicate the need of P1C1 DCBs as well
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   DCB = matrix containing DCB data
%
% DESCRIPTION:
%   Tool for loading .DCB files and providing P1P2 (and if needed P1C1) DCB data in output.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
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

if (isempty(constellations)) %then use only GPS as default
    [constellations] = goGNSS.initConstellation(1, 0, 0, 0, 0, 0);
end

%starting index in the total array for the various constellations
idGPS = constellations.GPS.indexes(1);
idGLONASS = constellations.GLONASS.indexes(1);
idGalileo = constellations.Galileo.indexes(1);
idBeiDou = constellations.BeiDou.indexes(1);
idQZSS = constellations.QZSS.indexes(1);

%output initialization
DCB = [];
DCB.P1C1 = [];
DCB.P1P2 = [];
DCB.P1C1.time  = 0;
DCB.P1C1.value = zeros(constellations.nEnabledSat, 1);
DCB.P1C1.rms   = zeros(constellations.nEnabledSat, 1);
DCB.P1C1.prn   = zeros(constellations.nEnabledSat, 1);
DCB.P1C1.sys   = zeros(constellations.nEnabledSat, 1);
DCB.P1P2.value = zeros(constellations.nEnabledSat, 1);
DCB.P1P2.rms   = zeros(constellations.nEnabledSat, 1);
DCB.P1P2.prn   = zeros(constellations.nEnabledSat, 1);
DCB.P1P2.sys   = zeros(constellations.nEnabledSat, 1);

%convert GPS time to time-of-week
gps_tow = weektime2tow(gps_week, time_R);

%detect starting and ending year/month
date = gps2date(gps_week, gps_tow);
year_start  = two_digit_year(date(1,1));
year_end    = two_digit_year(date(end,1));
month_start = date(1,2);
month_end   = date(end,2);

%directory containing DCB files
data_dir = dir(data_dir_dcb);

%check the number of files contained in the directory
nmax = size(data_dir,1);

%file counter
n = 0;
n_P1C1 = 0;
n_P1P2 = 0;

%find files with ".DCB" extension
for j = 1 : nmax

    %read the name of the j-th file
    dcb_file_name = getfield(data_dir,{j,1},'name');

    %get the number of characters in the filename
    dcb_fn_length = size(dcb_file_name,2);

    if (dcb_fn_length < 12)
        continue
    end

    year = str2num(dcb_file_name(5:6));
    month = str2num(dcb_file_name(7:8));

    %check if the filename corresponds to that expected from a standard DCB file required by goGPS (e.g. "P1C1xxyy.DCB",
    % with 'xx' = two-digit year and 'yy' = two-digit month)
    if ((strcmpi(dcb_file_name(1:4), 'P1P2') || strcmpi(dcb_file_name(1:4), 'P1C1')) && ...
       ((year >  year_start && year  <  year_end)    || ...
        (year == year_start && month >= month_start) && ...
        (year == year_end   && month <= month_end))  && ...
        ((dcb_fn_length == 12 && strcmpi(dcb_file_name(dcb_fn_length - 3 : dcb_fn_length), '.DCB')) || ...
         (dcb_fn_length == 16 && strcmpi(dcb_file_name(dcb_fn_length - 7 : dcb_fn_length), '.DCB_TMP')))) %#ok<*ST2NM>

        n = n + 1;

        switch dcb_file_name(3:4)
            case 'C1'
                n_P1C1 = n_P1C1 + 1;
            case 'P2'
                n_P1P2 = n_P1P2 + 1;
        end

        %full path to the target file
        dcb_file_target  = strcat(data_dir_dcb, '/', dcb_file_name);

        %open .dcb file
        fid_fd = fopen(dcb_file_target,'r');

        %warnings
        if (fid_fd ~= -1)
            %fprintf(['Reading DCB file ', dcb_file_name, '\n']);
            if (n == 1)
                fprintf('Reading DCB files...\n');
            end
        else
            fprintf(['WARNING: impossible to open DCB file ', dcb_file_name, '\n']);
            break
        end

        line = '';
        while(~feof(fid_fd) && ~strcmp(line, '***   ****************    *****.***   *****.***'))
            line = fgetl(fid_fd);
        end

        while(~feof(fid_fd))
            line = fgetl(fid_fd);

            if (isempty(line))
                continue
            end

            sys_id = line(1);
            if (strcmp(sys_id,'G') && constellations.GPS.enabled || ...
                strcmp(sys_id,'R') && constellations.GLONASS.enabled || ...
                strcmp(sys_id,'E') && constellations.Galileo.enabled || ...
                strcmp(sys_id,'C') && constellations.BeiDou.enabled || ...
                strcmp(sys_id,'J') && constellations.QZSS.enabled)

                PRN   = sscanf(line(2:3),'%f');
                value = sscanf(line(30:35),'%f');
                rms   = sscanf(line(43:47),'%f');

                switch (sys_id)
                    case 'G'
                        index = idGPS;
                        system = 'GPS';
                    case 'R'
                        index = idGLONASS;
                        system = 'GLONASS';
                    case 'E'
                        index = idGalileo;
                        system = 'Galileo';
                    case 'C'
                        index = idBeiDou;
                        system = 'BeiDou';
                    case 'J'
                        index = idQZSS;
                        system = 'QZSS';
                end

                if(~ismember(PRN,constellations.(system).PRN))
                    continue
                end

                index = index + PRN - 1;

                [w, s] = date2gps([four_digit_year(year) month 15 0 0 0]);

                switch dcb_file_name(3:4)
                    case 'C1'
                        DCB.P1C1.time(n_P1C1, 1)      = weektow2time(w, s, 'G');
                        DCB.P1C1.value(index, n_P1C1) = value;
                        DCB.P1C1.rms(index, n_P1C1)   = rms;
                        DCB.P1C1.prn(index)           = PRN;
                        DCB.P1C1.sys(index)           = sys_id;
                    case 'P2'
                        DCB.P1P2.time(n_P1P2, 1)      = weektow2time(w, s, 'G');
                        DCB.P1P2.value(index, n_P1P2) = value;
                        DCB.P1P2.rms(index, n_P1P2)   = rms;
                        DCB.P1P2.prn(index)           = PRN;
                        DCB.P1P2.sys(index)           = sys_id;
                end
            end
        end

        fclose(fid_fd);
    end
end

%if P1C1 data are needed but not available, return empty
if (any(codeC1_R(:)) && isempty(DCB.P1C1))
    DCB = [];
    fprintf(['The required P1C1 DCB file(s) were not found in ' data_dir_dcb ' directory.\n'])
    return
end

%if P1P2 data are needed but not available, return empty
if (isempty(DCB.P1P2))
    DCB = [];
    fprintf(['The required P1P2 DCB file(s) were not found in ' data_dir_dcb ' directory.\n'])
    return
end
