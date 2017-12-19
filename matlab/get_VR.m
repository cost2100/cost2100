function VR = get_VR( VRtable, paraEx, paraSt)
%GET_VR function to generate the VR
%Default call: VR = get_VR( VRtable, paraEx, paraSt)
%------
%Input:
%------
%paraEx,paraSt: external parameters and stochastic parameters
%VRtable: VR assignment table
%------
%Output:
%------
%VR: the VR distribution (numVR, [x y])
%
%See also: cost2100, get_para, get_VRtable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C)2008 LIU Ling-Feng, ICTEAM, UCL, Belgium 
%This file is part of cost2100.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nVR = length(VRtable(1,:,1)); %VR number
posBS = paraEx.pos_BS(1:2); %BS position (x,y) [m]
net_radii = paraEx.net_radii; %Radius of cell
rangeX = [posBS(1)-net_radii posBS(1)+net_radii]; %VR distribution range
rangeY = [posBS(2)-net_radii posBS(2)+net_radii]; %VR distribution range

angleVR = 2*pi* rand(nVR,1);
distVR  = sqrt(rand(nVR,1))*net_radii;
x       = distVR.*cos(angleVR)+posBS(1);
y       = distVR.*sin(angleVR)+posBS(2);            
VR(:,1) = x;
VR(:,2) = y;
    
% x = (rangeX(2)-rangeX(1))*rand(nVR,1)+rangeX(1); %Uniform distribution
% y = (rangeY(2)-rangeY(1))*rand(nVR,1)+rangeY(1); %Uniform distribution
% VR(:,1) = x;
% VR(:,2) = y;