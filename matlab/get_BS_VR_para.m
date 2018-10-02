function [BS_VR, BS_VR_len, BS_VR_pow_slope] = get_BS_VR_para(VRtable, paraEx, paraSt)
%GET_BS_VR_PARA Get parameters for BS visibility regions
%
%Default call:
%[BS_VR, BS_VR_len, BS_VR_pow_slope] = get_BS_VR_para(VRtable, paraEx, paraSt)
%
%------
%Input:
%------
%VRtable: Visibility region table
%paraEx: External parameters
%paraSt: Stochastic parameters
%------
%Output:
%------
%BS_VR: BS visibility regions
%BS_VR_len: BS visibility regions length
%BS_VR_pow_slope: BS visibility regions power slope

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

if (size(paraEx.pos_BS,1) > 1)
    error('Obs! So far, only support for single BS.');
end

Ncluster = max(VRtable(1,:,2)); % Number of clusters

% (1) Generate BS-VR length
BS_VR_len = mexprnd_own(paraSt.bs_vr_mu, paraEx.BS_range, Ncluster, 1); % [m]
%BS_VR_len = lognrnd_own(.7, 2, 1, Ncluster); % [m]

% (2) Determine BS-VR location
BS_VR = zeros(Ncluster, 2);
BS_VR(:,1) = paraEx.pos_BS(1, 1); % BS-VR location, x
BS_VR(:,2) = paraEx.pos_BS(1, 2); % BS-VR location, y
if (paraSt.bs_vr_mu < Inf)
    % Offset center of BS-VR along the x-axis.
    BS_VR_pos_span = paraEx.BS_range + BS_VR_len;
    BS_VR_pos_offset = BS_VR_pos_span.*(rand(Ncluster,1)-.5);
    BS_VR(:,1) = BS_VR(:,1) +  BS_VR_pos_offset; % BS-VR location, x
end

% (3) BS-VR power slope
BS_VR_pow_slope = randn(Ncluster,1).*paraSt.bs_vr_slope_sigma+paraSt.bs_vr_slope_mu; % Normal distribution [dB/m]