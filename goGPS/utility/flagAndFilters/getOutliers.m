function [flagIntervals] = getOutliers(flags)
% SYNTAX:
%    [flagIntervals] = getOutliers(flags)
%
% DESCRIPTION:
%    returns start and end of flagged intervals
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti
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


    % convert flags into array
    if isstruct(flags)
        flagArray = int8(struct2flagVec(flags, max(flags.pos+1)));
        % add padding to avoid problem with flags on the borders
    else
        flagArray = flags;
    end
    flagArray = flagArray(:);
    flagArray = [0; flagArray; 0];
    flagArray(flagArray ~= 0) = 1;
    diff = flagArray(1:end-1) - flagArray(2:end);
    clear flagArray;
    flagIntervals = [find(diff<0), find(diff>0)-1];
end

function [flagArray] = struct2flagVec(flags, maxSize)
flagArray = false(maxSize,1);
flagArray(flags.pos) = flags.val;
end
