function [download_successful, compressed] = download_nav(filename, nav_path)

% SYNTAX:
%   [download_successful, compressed] = download_nav(filename, nav_path);
%
% INPUT:
%   filename = name of the RINEX navigation file to be downloaded
%   (brdcddd0.yyn, brdmddd0.yyp and CGIMddd0.yyN supported)
%   nav_path = download path
%
% OUTPUT:
%   download_successful = flag to identify unsuccessful downloads
%   compressed = flag to let the calling function know whether the
%                downloaded files are still compressed
%
% DESCRIPTION:
%   Download of RINEX navigation files from the IGS or AIUB FTP servers.

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

% Pointer to the global settings:
state = Go_State.getCurrentSettings();

download_successful = 0;
compressed = 0;

%IGS FTP server URL
igs_url = 'cddis.gsfc.nasa.gov';

%AIUB FTP server URL
aiub_url = 'ftp.aiub.unibe.ch';

%download directory
if (nargin < 2)
    down_dir = state.eph_dir;
else
    down_dir = nav_path;
end

% Check / create output folder
if not(exist(down_dir, 'dir'))
    mkdir(down_dir);
end

%identify requested file type
if (strcmp(filename(1:4),'brdc'))
    url = igs_url;
    name = 'IGS';
    path = '/pub/gps/data/daily/';
    subdir = '/brdc/';
elseif (strcmp(filename(1:4),'brdm'))
    url = igs_url;
    name = 'IGS';
    path = '/pub/gps/data/campaign/mgex/daily/rinex3/';
    subdir = '/brdm/';
elseif (strcmp(filename(1:4),'CGIM'))
    url = aiub_url;
    name = 'AIUB';
    path = '/aiub/CODE/';
    subdir = '';
else
    error('Only "brdc", "brdm" (IGS) and "CGIM" (AIUB) files are supported.');
end

fprintf(['FTP connection to the ' name ' server (ftp://' url '). Please wait...'])

%connect to the server
try
    ftp_server = ftp(url);
catch
    fprintf(' connection failed.\n');
    return
end

fprintf('\n');

%target directory
s = [path num2str(four_digit_year(str2num(filename(end-2:end-1)))) subdir];

cd(ftp_server, '/');
cd(ftp_server, s);

filename = [filename '.Z'];

try
    mget(ftp_server,filename,down_dir);
    if (isunix())
        system(['uncompress -f ' down_dir filesep filename]);
    else
        try
            [status, result] = system(['".\utility\thirdParty\7z1602-extra\7za.exe" -y x ' '"' down_dir filename '"' ' -o' '"' down_dir '"']); %#ok<ASGLU>
            delete([down_dir filename]);
            filename = filename(1:end-2);
        catch
            fprintf(['Please decompress the ' filename ' file before trying to use it in goGPS.\n']);
            compressed = 1;
        end
    end
    fprintf(['Downloaded ' filename(1:4) ' file: ' filename '\n']);
    download_successful = 1;
catch
end

close(ftp_server);

fprintf('Download complete.\n')
