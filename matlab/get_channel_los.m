function h_los = get_channel_los( channel, BS, MS, paraEx, paraSt)
%GET_LOS Compute the LOS of the channel in BS-MS link
%Default call: h_los = get_channel_los( channel, BS, MS, paraEx, paraSt)
%-------
%Input:
%-------
%channel: channel recording of link BS-MS
%BS: BS information
%MS: MS information
%paraEx, paraSt: external and stochastic parameters
%------
%Output:
%------
%h_los: LOS channel reponse
%
%See also: cost2100, get_para, get_VRLOS, calc_pathloss, get_channel

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


lambda = paraEx.c0/paraEx.freq;%Wavelength

if calc_dist(BS.pos,MS.pos)>paraSt.d_co
    power_los = 0;  
    power_factor = 0;
else
    d_MS_VRLOS = calc_dist(MS.pos(1:2),BS.VRLOS);%Get the transition function
    if d_MS_VRLOS >paraSt.r_l
        power_los = 0;
        power_factor = 0;
    else
        y = paraSt.l_l+d_MS_VRLOS-paraSt.r_l;
        A_los = 1/2-atan(2*sqrt(2)*y/sqrt(lambda*paraSt.l_l))/pi;

        power_other = abs(sum(channel.h(:,6)))^2;
        %OR
        %power_other = sqrt(sum(abs(channel.h(:,6)).^2));
        
        power_factor = paraSt.power_factor;        
        power_los = abs(power_factor)*power_other;
    end
end

[phi_BS_MS,theta_BS_MS,tmp]=cart2sph(MS.pos(1)-BS.pos(1), MS.pos(2)-BS.pos(2), MS.pos(3)-BS.pos(3));
[phi_MS_BS,theta_MS_BS,tmp]=cart2sph(BS.pos(1)-MS.pos(1), BS.pos(2)-MS.pos(2), BS.pos(3)-MS.pos(3));
d_BS_MS = calc_dist(BS.pos(1:2),MS.pos(1:2)); %distance BS to MS, in X-Y plane
tau_0=d_BS_MS/paraEx.c0; % delay of the LOS
h_los = [phi_BS_MS theta_BS_MS phi_MS_BS theta_MS_BS tau_0 sqrt(power_los)*exp(-j*2*pi*d_BS_MS/paraEx.c0*paraEx.freq)]; %The unfiltered LOS impulse response, angle in radii;