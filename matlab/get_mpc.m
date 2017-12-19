function mpc = get_mpc(cluster, paraSt)
%GET_MPC Generate MPCs for clusters
%Default call: mpc = get_mpc(cluster, paraSt)
%-------
%Input:
%-------
%cluster: cluster information
%paraSt: stochastic parameters
%-------
%Output:
%-------
%mpc: mpc information
% .pos_BS(numMPC,[x y z]): position of MPC at BS side
% .pos_MS(numMPC,[x y z]): position of MPC at MS side
% .a_mpc(numMPC): complex fading of each mpc, normalized
%
%See also: get_cluster, get_cluster_local, get_para, cost2100

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

n_c = length(cluster); %Number of far clusters
n_mpc = paraSt.n_mpc; %Number of MPCs per cluster

for m = 1:n_c
    r = randn(n_mpc,3); %MPC distribution
    mpc_tmp = r*diag([cluster(m).a_c_BS/3 cluster(m).b_c_BS/3 cluster(m).h_c_BS/3]);
    
    switch cluster(m).type
        case 0 %Local cluster
            mpc_BS_tmp = mpc_tmp;
            mpc_BS_tmp = mpc_BS_tmp + repmat(cluster(m).pos_c_BS,n_mpc,1);
            mpc_MS_tmp = mpc_BS_tmp;
        case 1 %Single cluster
            mpc_BS_tmp = mpc_tmp*rotate_matrix(cluster(m).Phi_c_BS,cluster(m).Theta_c_BS);
            mpc_BS_tmp = mpc_BS_tmp + repmat(cluster(m).pos_c_BS,n_mpc,1);
            mpc_MS_tmp = mpc_BS_tmp;
        case 2 %Twin cluster
            mpc_BS_tmp = mpc_tmp*rotate_matrix(cluster(m).Phi_c_BS,cluster(m).Theta_c_BS);
            mpc_BS_tmp = mpc_BS_tmp + repmat(cluster(m).pos_c_BS,n_mpc,1);
            mpc_MS_tmp = mpc_tmp*rotate_matrix(cluster(m).Phi_c_MS,cluster(m).Theta_c_MS);
            mpc_MS_tmp = mpc_MS_tmp + repmat(cluster(m).pos_c_MS,n_mpc,1);                   
    end
    mpc(m).pos_BS = mpc_BS_tmp;    
    mpc(m).pos_MS = mpc_MS_tmp;
    
    P_S = randn(1,n_mpc)+j*randn(1,n_mpc); %Complex Rayleigh distribution
    a_mpc = P_S/sqrt(sum(abs(P_S).^2)); %Normalized attenuation of each MPC's complex amplitude
    mpc(m).a_mpc = a_mpc; %Complex amplitude of MPC (flat fading)    
    mpc(m).idx_c=m;    
end