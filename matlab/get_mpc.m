function mpc = get_mpc(cluster, paraSt)
%GET_MPC Generate MPCs for clusters
%
%Default call: 
%mpc = get_mpc(cluster, paraSt)
%-------
%Input:
%-------
%cluster: Cluster information
%paraSt: Stochastic parameters
%-------
%Output:
%-------
%mpc: mpc information
% .pos_BS(numMPC,[x y z]): Position of MPC at BS side
% .pos_MS(numMPC,[x y z]): Position of MPC at MS side
% .a_mpc(numMPC): Complex fading of each MPC, normalized

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

n_c = length(cluster); % Number of far clusters
n_mpc = paraSt.n_mpc; % Number of MPCs per cluster

for m = 1:n_c % Loop for number of clusters
    r = randn(n_mpc, 3); % MPC distribution
    mpc_tmp = r*diag([cluster(m).a_c_BS/3 cluster(m).b_c_BS/3 cluster(m).h_c_BS/3]); % [spatial delay spread, width of cluster, height of cluster]
    
    switch cluster(m).type
        case 0 % Local cluster
            mpc_BS_tmp = mpc_tmp;
            mpc_BS_tmp = mpc_BS_tmp + repmat(cluster(m).pos_c_BS, n_mpc,1);
            mpc_MS_tmp = mpc_BS_tmp; 
        case 1 % Single cluster
            mpc_BS_tmp = mpc_tmp*rotate_matrix(cluster(m).Phi_c_BS, cluster(m).Theta_c_BS); % According to cluster spread
            mpc_BS_tmp = mpc_BS_tmp + repmat(cluster(m).pos_c_BS, n_mpc, 1); % Refer to cluster position
            mpc_MS_tmp = mpc_BS_tmp;
        case 2 % Twin cluster
            mpc_BS_tmp = mpc_tmp*rotate_matrix(cluster(m).Phi_c_BS, cluster(m).Theta_c_BS);
            mpc_BS_tmp = mpc_BS_tmp + repmat(cluster(m).pos_c_BS, n_mpc, 1);
            mpc_tmp = r*diag([cluster(m).a_c_MS/3 cluster(m).b_c_MS/3 cluster(m).h_c_MS/3]); % [spatial delay spread, width of cluster, height of cluster]
            mpc_MS_tmp = mpc_tmp*rotate_matrix(cluster(m).Phi_c_MS, cluster(m).Theta_c_MS);
            mpc_MS_tmp = mpc_MS_tmp + repmat(cluster(m).pos_c_MS, n_mpc, 1);                  
    end
    mpc(m).pos_BS = mpc_BS_tmp;    
    mpc(m).pos_MS = mpc_MS_tmp;
    
    % Amplitudes shall be either Rayleigh-faded or modulated by the gain
    % function depending on the network parameters
    if (paraSt.mpc_gain_mu == Inf)
        P_S = randn(1, n_mpc)+1j*randn(1, n_mpc); % Complex Rayleigh distribution
        a_mpc = P_S/sqrt(sum(abs(P_S).^2)); % Normalized attenuation of each MPCs, complex amplitude
    else
        a_mpc = sqrt(1/n_mpc)*exp(1j*2*pi*rand(1, n_mpc));
    end
    mpc(m).a_mpc = a_mpc; % Complex amplitude of MPC (flat fading)    
    mpc(m).idx_c = m; % Belongs to which cluster
    
    % Mapping to the center of MPC gain function in cluster MS-VR
    % Origin at the center of cluster MS-VR
    angle_gain_center = angle(exp(1j*2*pi*rand(n_mpc, 1)));
    % MPC gain peak location - relative to the cluster VR center
    if cluster(m).type==0 % Local cluster
        mpc(m).gain_center = mpc(m).pos_BS(:,1:2); % [x, y], consistent with update_chan.m
        mpc(m).gain_radius = lognrnd_own(paraSt.mpc_gain_mu*(log(10)/10), paraSt.mpc_gain_sigma*(log(10)/10), n_mpc, 1);
    else % Far cluster
        rr = rand(n_mpc, 1);
        mpc(m).gain_center = [paraSt.r_c*sqrt(rr).*cos(angle_gain_center), paraSt.r_c*sqrt(rr).*sin(angle_gain_center)];
        mpc(m).gain_radius = lognrnd_own(paraSt.mpc_gain_mu*(log(10)/10), paraSt.mpc_gain_sigma*(log(10)/10), n_mpc, 1);
        mpc(m).gain_radius(mpc(m).gain_radius<paraSt.l_c) = paraSt.l_c;
        mpc(m).gain_radius(mpc(m).gain_radius>paraSt.r_c & mpc(m).gain_radius<Inf) = paraSt.r_c;
    end
    
    % MPC polarization matrix; see (3.45) in COST 2100 book.
    XPDV = 10.^((paraSt.mu_xpdv+paraSt.sigma_xpdv.*randn(n_mpc, 1))/10);
    XPDH = 10.^((paraSt.mu_xpdh+paraSt.sigma_xpdh.*randn(n_mpc, 1))/10);
    CPR = 10.^((paraSt.mu_cpr+paraSt.sigma_cpr.*randn(n_mpc, 1))/10);
    mpc_vv = exp(1j*2*pi*rand(n_mpc, 1));
    mpc_vh = 1./sqrt(XPDV).*exp(1j*2*pi*rand(n_mpc, 1));
    mpc_hv = 1./sqrt(XPDH.*CPR).*exp(1j*2*pi*rand(n_mpc, 1));
    mpc_hh = 1./sqrt(CPR).*exp(1j*2*pi*rand(n_mpc, 1));
    
    % Polarization power normalization
    for idx_mpc = 1:n_mpc
         norm_factor = norm([mpc_vv(idx_mpc), mpc_vh(idx_mpc); mpc_hv(idx_mpc), mpc_hh(idx_mpc)], 'fro');
         
         mpc(m).vv(idx_mpc) = mpc_vv(idx_mpc)/norm_factor;
         mpc(m).vh(idx_mpc) = mpc_vh(idx_mpc)/norm_factor;
         mpc(m).hh(idx_mpc) = mpc_hh(idx_mpc)/norm_factor;
         mpc(m).hv(idx_mpc) = mpc_hv(idx_mpc)/norm_factor;
    end
end