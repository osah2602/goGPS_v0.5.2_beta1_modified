function dockAllFigures(fig_handle)
% Change the WindowStyle of all the figures to docked
% SYNTAX:
%
% SINTAX:
%   dockAllFigures(<fig_handle>);
%   dockAllFigures()
%
% EXAMPLE:
%   dockAllFigures(gcf);
%
% INPUT:
%   fig_handle = handler to the axes to modify           <optional argument>
%
% DEFAULT VALUES:
%   fig_handle all the figures
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
%  Contributors:     Andrea Gatti ...
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
    if (nargin == 1)
        set(fig_handle,'WindowStyle','docked');
    else
        figHandles = findall(0,'Type','figure');
        for fig_handle = 1:length(figHandles)
            set(figHandles(fig_handle),'WindowStyle','docked')
        end
    end
end
