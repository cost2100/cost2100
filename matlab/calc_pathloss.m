function pathloss = calc_pathloss(pos_BS, pos_MS, paraEx)
%CALC_PATHLOSS Function to calculate the pathloss between BS and MS

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dist = calc_dist(pos_BS,pos_MS);

switch paraEx.network
    case 'macrocell' % The Walfisch-Ikegami-Model
        delta_base = paraEx.h_BS-paraEx.h_rooftop;
        pl_0 = 32.4+20*log10( dist/1000 )+20*log10(paraEx.freq/1e6);
        phi_road = paraEx.phi_road*180/pi;
        if phi_road>=0&&phi_road<35
            pl_ori = -10+0.354*phi_road;
        elseif phi_road>=35&&phi_road<55
            pl_ori = 2.5+0.075*(phi_road-35);
        elseif phi_road>=55&&phi_road<90
            pl_ori = 4.0-0.114*(phi_road-55);
        end
        pl_lts = -16.9-10*log10(paraEx.w_r)+10*log10(paraEx.freq/1e6)+20*log10(paraEx.h_rooftop-paraEx.h_MS)+pl_ori;
        if paraEx.h_BS>paraEx.h_rooftop
            pl_bsh = -18*log10(1+delta_base);
        else
            pl_bsh = 0;
        end
        if paraEx.h_BS>paraEx.h_rooftop
            k_a = 54;
        elseif paraEx.h_BS<=paraEx.h_rooftop && d>=500
            k_a = 54-0.8*delta_base;
        elseif paraEx.h_BS<=paraEx.h_rooftop && d<500
            k_a = 54-0.8*delta_base*dist/500;
        end
        if paraEx.h_BS>paraEx.h_rooftop
            k_d = 18;
        else
            k_d = 18-15*delta_base/(paraEx.h_rooftop);
        end
        k_f = -4+0.7*(paraEx.freq/1e6/925-1);
        %         k_f = 1.5*(paraEx.freq/1e6/925-1);
        pl_msd = pl_bsh+k_a+k_d*log10(dist/1000)+k_f*log10(paraEx.freq/1e6)-9*log10(paraEx.w_b);% building seperation?        
        if pl_msd+pl_lts >= 0

            pathloss = pl_0+pl_msd+pl_lts;
        else
            pathloss = pl_0;
        end
        pathloss = 10^(-pathloss/20); %Don't forget The complex phase
    case 'microcell' %COST259(GSN)        
        pl = 10*2.6*log10(dist)+20*log10(4*pi*1/(paraEx.c0/paraEx.freq)); 
        pathloss = 10^(-pl/20);
    case 'picocell' %COST2100 TD02-055
        pl_0 = 20*log10(dist)+20*log10(paraEx.freq)-147.56; %Free space pathloss
        pl_d = paraEx.n_floor*30;
        pathloss = 10^(-(pl_0+pl_d)/20);
    case 'IndoorHall_5GHz' %IndoorHall_5GHz
        pl_0 = 20*log10(dist)+20*log10(paraEx.freq)-147.56; %Free space pathloss
        pl_d = paraEx.n_floor*30; %Average number of floors * 30 dB extra pathloss
        pathloss = 10^(-(pl_0+pl_d)/20);
    case 'SemiUrban_300MHz' %SemiUrban_300MHz
        pl_0 = 20*log10(dist)+20*log10(paraEx.freq)-147.56; %Free space pathloss
        pl_d = paraEx.n_floor*30; %Average number of floors * 30 dB extra pathloss
        pathloss = 10^(-(pl_0+pl_d)/20);        
end        