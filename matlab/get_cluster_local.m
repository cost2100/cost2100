function cluster = get_cluster_local(type, ind, paraEx, paraSt)
%GET_CLUSTER_LOCAL Generate the local cluster at BS/MS side
%Default call: 
%cluster = get_cluster_local(type, ind, paraEx, paraSt)
%-------
%Input:
%-------
%type: Local cluster type: 'BS', 'MS'
%ind: Index of BS/MS for the local cluster
%paraEx: External parameters
%paraSt: Stochastic parameters
%------
%Output:
%------
%cluster: local cluster information, the same structure as of far clusters

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

corr_randn = randn(1, 6)*paraSt.corr_mat;

shadow_f = 10.^(0.1*corr_randn(1)*paraSt.sigma_sf); % Shadow fading
tau_c = paraSt.mu_tau.*10.^(0.1*paraSt.sigma_tau*corr_randn(2)); % Delay spread
theta_c_BS = paraSt.mu_theta_BS*10.^(0.1*paraSt.sigma_theta_BS*corr_randn(3)); % Elevation spread BS
phi_c_BS = paraSt.mu_phi_BS*10.^(0.1*paraSt.sigma_phi_BS*corr_randn(4)); % Azimtuh spread BS
theta_c_MS = paraSt.mu_theta_MS*10.^(0.1*paraSt.sigma_theta_MS*corr_randn(5)); % Elevation spread MS
phi_c_MS = paraSt.mu_phi_MS*10.^(0.1*paraSt.sigma_phi_MS*corr_randn(6)); % Azimuth spread MS
d_tau = tau_c*paraEx.c0/2; % Spatial delay spread

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The spread, shadowing, power, delay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cluster.a_c_BS = d_tau; % Delay spread
cluster.b_c_BS = d_tau; % Width of cluster [m]
cluster.h_c_BS = d_tau/20; % Height of cluster
cluster.a_c_MS = cluster.a_c_BS; 
cluster.b_c_MS = cluster.b_c_BS;
cluster.h_c_MS = d_tau/20;
cluster.shadow_f = shadow_f;
cluster.tau_c_link = 0.; % Single bounce for the tau_c_link
        
switch type
    case 'BS' % Local cluster at BS side                
        pos_BS = paraEx.pos_BS(ind,:);
        cluster.pos_c_BS = pos_BS; % Determine cluster at BS side
        cluster.pos_c_MS = pos_BS; % Determine cluster at MS side
        cluster.type = 0;
        
        cluster.active = paraSt.BSLocal;
    case 'MS' % Local cluster at MS side
        pos_MS = paraEx.pos_MS(ind,:);
        cluster.pos_c_BS = pos_MS; % Determine cluster at BS side
        cluster.pos_c_MS = pos_MS; % Determine cluster at MS side
        cluster.type = 0;
        
        cluster.active = paraSt.MSLocal;
        
    otherwise 
        error('Local cluster type error.\n');
end
