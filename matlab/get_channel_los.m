function h_los = get_channel_los(channel, BS, BS_pos, MS, paraEx, paraSt)
%GET_CHANNEL_LOS Compute the LOS of the channel in BS-MS link
%
%Default call: 
%h_los = get_channel_los(channel, BS, BS_pos, MS, paraEx, paraSt)
%-------
%Input:
%-------
%channel: channel recording of link BS-MS
%BS: BS information
%BS_pos: Position at BS
%MS: MS information
%paraEx External parameters
%paraSt: Stochastic parameters
%------
%Output:
%------
%h_los: LOS channel reponse

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

% Consider just one position
BS.pos = BS_pos;

lambda = paraEx.c0/paraEx.center_freq; % Wavelength

if (calc_dist(BS.pos,MS.pos)>paraSt.d_co)
    power_los = 0;  
    power_factor = 0;
else
    % Get the transition function
    d_MS_VRLOS = calc_dist(MS.pos(1:2),BS.VRLOS);
    if (d_MS_VRLOS >paraSt.r_l)
        power_los = 0;
        power_factor = 0;
    else
        y = d_MS_VRLOS-(paraSt.r_l-paraSt.l_l);
        A_los = 1/2-atan(2*sqrt(2)*y/sqrt(lambda*paraSt.l_l))/pi;

        power_other = sum(abs(channel.h(:,6)).^2);
        
        power_factor = paraSt.power_factor;        
        power_los = abs(power_factor)*power_other;
    end
end

[phi_BS_MS,theta_BS_MS,~]=cart2sph(MS.pos(1)-BS.pos(1), MS.pos(2)-BS.pos(2), MS.pos(3)-BS.pos(3));
[phi_MS_BS,theta_MS_BS,~]=cart2sph(BS.pos(1)-MS.pos(1), BS.pos(2)-MS.pos(2), BS.pos(3)-MS.pos(3));
d_BS_MS = calc_dist(BS.pos(1:2),MS.pos(1:2)); % Distance BS to MS, obs in X-Y plane!
tau_0=d_BS_MS/paraEx.c0; % Delay of the LOS

% Magnitude of the MPC amplitude
amp_LOS = sqrt(power_los);

% Polarization matrix
amp_vv = 1/sqrt(2);
amp_vh = 0;
amp_hh = 1/sqrt(2);
amp_hv = 0;

% The unfiltered LOS impulse response, angle in radians
h_los = [phi_BS_MS theta_BS_MS phi_MS_BS theta_MS_BS tau_0 amp_LOS amp_vv amp_vh amp_hh amp_hv]; 
