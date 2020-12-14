function [ String ] = RemoveUnwantedTrailingSpaces( String )

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Damiano Triglione
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

%REMOVEUNWANTEDTRAILINGSPACES Summary of this function goes here
%   Detailed explanation goes here

    LengthOfString = length(String);
    IndexOfLatestCharInString = LengthOfString;

    while (IndexOfLatestCharInString >= 1) && isspace(String(IndexOfLatestCharInString))
        IndexOfLatestCharInString = IndexOfLatestCharInString-1;
    end %while

    if IndexOfLatestCharInString == 0
        String = [];
    else
        String = String(1:IndexOfLatestCharInString);
    end %if

end

