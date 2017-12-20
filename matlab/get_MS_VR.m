function MS_VR = get_MS_VR(VRtable, paraEx)
%GET_MS_VR Function to generate the VR
%
%Default call: 
%MS_VR = get_MS_VR(VRtable, paraEx)
%------
%Input:
%------
%paraEx: External parameters
%VRtable: VR assignment table
%------
%Output:
%------
%VR: the VR distribution (numVR, [x y])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This file is a part of the COST2100 channel model.
%
%This program, the COST2100 channel model, is free software: you can 
%redistribute it and/or modify it under the terms of the GNU General Public 
%License as published by the Free Software Foundation, either version 3 of 
%the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful, but 
%WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
%or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
%for more details.
%
%If you use it for scientific purposes, please consider citing it in line 
%with the description in the Readme-file, where you also can find the 
%contributors.
%
%You should have received a copy of the GNU General Public License along 
%with this program. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nVR = length(VRtable(1,:,1)); % VR number
posBS = paraEx.pos_BS(1:2); % BS position (x,y) [m]
net_radii = paraEx.net_radii; % Radius of cell
rangeX = [posBS(1)-net_radii posBS(1)+net_radii]; % VR distribution range, x
rangeY = [posBS(2)-net_radii posBS(2)+net_radii]; % VR distribution range, y

angleVR = 2*pi*rand(nVR, 1);
distVR  = sqrt(rand(nVR, 1))*net_radii;
x       = distVR.*cos(angleVR)+posBS(1);
y       = distVR.*sin(angleVR)+posBS(2);            
MS_VR(:, 1) = x;
MS_VR(:, 2) = y;
