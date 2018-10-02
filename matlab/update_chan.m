function MS_new = update_chan(MS, paraEx, paraSt)
%UPDATE_CHAN Update the channel according to the MS movement
%
%Default call: 
%MS_new = update_chan(MS, paraEx,paraSt)
%
%-------
%Input:
%-------
%MS: MS information
%paraEx: External parameters 
%paraSt: Stochastic parameters
%
%------
%Output:
%------
%MS_new: Updated MS information

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

MS_new = MS;
pos_MS_old = MS.pos;
pos_MS_new = pos_MS_old + MS.velo/paraEx.snap_rate;

% Only update the local cluster at MS side 
d_local_mpc_MS = sqrt(sum((repmat(pos_MS_new(1:2),paraSt.n_mpc,1) - MS.mpc_local.pos_MS(:,1:2)).^2,2));

replace_mpc = find(d_local_mpc_MS > MS.cluster_local.a_c_MS); % Find MPCs out range of the local cluster

if length(replace_mpc)==0 % Then do nothing
    
else % Update MPCs in local cluster
    for m = 1:length(replace_mpc)
        new_mpc(m, :) = randn(1, 3)*diag([MS.cluster_local.a_c_MS/3 MS.cluster_local.b_c_MS/3,MS.cluster_local.h_c_MS/3])+pos_MS_new;
        while calc_dist(new_mpc(m, :),pos_MS_new) > MS.cluster_local.a_c_MS
            new_mpc(m, :) = randn(1, 3)*diag([MS.cluster_local.a_c_MS/3 MS.cluster_local.b_c_MS/3,MS.cluster_local.h_c_MS/3])+pos_MS_new;
        end
    end
    MS_new.mpc_local.pos_BS(replace_mpc,:) = new_mpc; % Copy out-range MPCs
    MS_new.mpc_local.pos_MS(replace_mpc,:) = new_mpc;
    
    % MPC gain function peak position
    MS_new.mpc_local.gain_center(replace_mpc,:) = MS_new.mpc_local.pos_BS(replace_mpc,1:2); % [x, y]
    MS_new.mpc_local.gain_radius(replace_mpc,:) = lognrnd_own(paraSt.mpc_gain_mu*(log(10)/10), paraSt.mpc_gain_sigma*(log(10)/10), length(replace_mpc), 1);
end

MS_new.cluster_local.p_c_MS = pos_MS_new;
MS_new.cluster_local.p_c_BS = pos_MS_new;
MS_new.pos = pos_MS_new;