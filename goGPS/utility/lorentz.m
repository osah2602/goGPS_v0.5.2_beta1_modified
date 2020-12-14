function [p] = lorentz(x,y)

% SYNTAX:
%   [p] = lorentz(x,y);
%
% INPUT:
%   x = vector of size 4x1
%   y = vector of size 4x1
%
% OUTPUT:
%   p = Lorentz inner product of the two 4x1 vectors
%
% DESCRIPTION:
%   Computation of the Lorentz inner product.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) Kai Borre
%  Written by:       Kai Borre 04-22-95
%  Contributors:     Mirko Reguzzoni, Eugenio Realini, 2009
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

p = x(1)*y(1) + x(2)*y(2) + x(3)*y(3) - x(4)*y(4);
