function [sys, prn] = find_sat_system(sat, constellations)

% SYNTAX:
%   [sys, prn] = find_sat_system(sat, constellations);
%
% INPUT:
%   sat = satellite array indexes
%   constellations = struct with multi-constellation settings
%                   (see goGNSS.initConstellation - empty if not available)
%
% OUTPUT:
%   sys = system id (e.g. 'GGGRRREEJ') [string]
%   prn = system-specific PRN code
%
% DESCRIPTION:
%   Given a list of satellite array indexes, find to which system they
%   refer to, and assign the correct PRN code.

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

sys = zeros(1,length(sat));
prn = zeros(1,length(sat));

for j = 1 : length(sat)
    if (ismember(sat(j),constellations.GPS.indexes))
        sys(1,j) = 'G';
        prn(1,j) = find(constellations.GPS.indexes == sat(j));
    elseif (ismember(sat(j),constellations.GLONASS.indexes))
        sys(1,j) = 'R';
        prn(1,j) = find(constellations.GLONASS.indexes == sat(j));
    elseif (ismember(sat(j),constellations.Galileo.indexes))
        sys(1,j) = 'E';
        prn(1,j) = find(constellations.Galileo.indexes == sat(j));
    elseif (ismember(sat(j),constellations.BeiDou.indexes))
        sys(1,j) = 'C';
        prn(1,j) = find(constellations.BeiDou.indexes == sat(j));
    elseif (ismember(sat(j),constellations.QZSS.indexes))
        sys(1,j) = 'J';
        prn(1,j) = find(constellations.QZSS.indexes == sat(j));
    elseif (ismember(sat(j),constellations.SBAS.indexes))
        sys(1,j) = 'S';
        prn(1,j) = find(constellations.SBAS.indexes == sat(j));
    end
end

sys = char(sys);
