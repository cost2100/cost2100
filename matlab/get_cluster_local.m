function cluster = get_cluster_local(type, ind, paraEx,paraSt)
%GET_CLUSTER_LOCAL Generate the local cluster at BS/MS side
%default call: cluster = get_cluster_local(type, ind, paraEx,paraSt)
%-------
%Input:
%-------
%type: local cluster type: 'BS', 'MS'
%ind: index of BS/MS for the local cluster
%paraEx, paraSt: external and stochastic parameters
%------
%Output:
%------
%cluster: local cluster information, the same structure of far clusters
%
%See also: get_cluster, cost2100, get_mpc, get_para

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


corr_randn = randn(1,6)*paraSt.corr_mat;

shadow_f = 10.^(0.1*corr_randn(1)*paraSt.sigma_sf); %Shadow fading
tau_c = paraSt.mu_tau.*10.^(0.1*paraSt.sigma_tau*corr_randn(2)); %Delay spread
theta_c_BS = paraSt.mu_theta_BS*10.^(0.1*paraSt.sigma_theta_BS*corr_randn(3)); %elevation spread BS
phi_c_BS = paraSt.mu_phi_BS*10.^(0.1*paraSt.sigma_phi_BS*corr_randn(4)); %azimtuh spread BS
theta_c_MS = paraSt.mu_theta_MS*10.^(0.1*paraSt.sigma_theta_MS*corr_randn(5)); %elevation spread MS
phi_c_MS = paraSt.mu_phi_MS*10.^(0.1*paraSt.sigma_phi_MS*corr_randn(6)); %azimuth spread MS
d_tau = tau_c*paraEx.c0/2; % spatial delay spread

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The spread, shadowing, power, delay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cluster.a_c_BS = d_tau; %Delay spread
cluster.b_c_BS = d_tau; %Azimuth spread
cluster.h_c_BS = d_tau/20; % Elevation spread
cluster.a_c_MS = cluster.a_c_BS;
cluster.b_c_MS = cluster.b_c_BS;
cluster.h_c_MS = d_tau/20;
cluster.shadow_f = shadow_f;
cluster.tau_c_link = 0.;%Single bounce for the tau_c_link
        
switch type
    case 'BS'                
        pos_BS = paraEx.pos_BS(ind,:);
        cluster.pos_c_BS = pos_BS; %Determine cluster at BS side
        cluster.pos_c_MS = pos_BS; %Determin cluster at MS side
        cluster.type = 0;
        
        cluster.active = paraSt.BSLocal;
    case 'MS'
        pos_MS = paraEx.pos_MS(ind,:);
        cluster.pos_c_BS = pos_MS; %Determine cluster at BS side
        cluster.pos_c_MS = pos_MS; %Determin cluster at MS side
        cluster.type = 0;
        
        cluster.active = paraSt.MSLocal;
        
    otherwise 
        error('Local cluster type error.\n');
end
