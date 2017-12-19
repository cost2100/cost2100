function matrix = rotate_matrix(phi,theta) % azi ele

%ROTATE_MATRIX Function to generate a 3x3 rotation matrix 
%matrix = rotate_matrix(phi,theta) return a 3x3 rotation matrix by given
%rotation in azimuth phi and elevation theta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C)2008 LIU Ling-Feng, ICTEAM, UCL, Belgium 
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% matrix = [cos(phi)*cos(theta) -sin(phi) cos(phi)*sin(theta);...
%           sin(phi)*cos(theta)  cos(phi) sin(phi)*cos(theta);...
%           -sin(theta)         0         cos(theta)];

delta = 0;
Tx=[1    0           0       ;...
    0  cos(delta) -sin(delta);...
    0  sin(delta)  cos(delta)];
        
Ty=[cos(theta) 0 sin(theta);...
       0     1   0     ;...
   -sin(theta) 0 cos(theta)];

Tz=[cos(phi)  sin(phi) 0 ;...
    -sin(phi)  cos(phi) 0 ;...
        0          0       1];

matrix=Tx*Ty*Tz;

