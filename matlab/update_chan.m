function MS_new = update_chan( MS, paraEx,paraSt )
%UPDATE_CHAN Update the channel according to the MS movement.
%
%Default call: MS_new = update_chan( MS, paraEx,paraSt )
%
%-------
%Input:
%-------
%
%MS: MS information
%
%paraEx, paraSt: external and stochastic parameters
%
%------
%Output:
%------
%MS_new: updated MS information
%
%See also: cost2100, get_para, get_cluster, get_MPC

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

MS_new = MS;
pos_MS_old = MS.pos;
pos_MS_new = pos_MS_old + MS.velo*paraEx.snap_rate.^(-1);

for m = 1:paraSt.n_mpc %Only update the local cluster at MS side 
    d_local_mpc_MS(m) = calc_dist(MS.mpc_local.pos_MS(m,1:2),pos_MS_new(1:2)'); %Distance of old mpcs to new MS    
end

replace_mpc = find (d_local_mpc_MS > MS.cluster_local.a_c_MS); %Find MPCs out range of the local cluster

if length(replace_mpc)==0 %Then do nothing
else %Update
    for m = 1:length(replace_mpc)
        new_mpc(m,:) = randn(1,3)*diag([MS.cluster_local.a_c_MS/3 MS.cluster_local.b_c_MS/3,MS.cluster_local.h_c_MS/3])+pos_MS_new;
        while calc_dist(new_mpc(m,:),pos_MS_new)<MS.cluster_local.a_c_MS
            new_mpc(m,:) = randn(1,3)*diag([MS.cluster_local.a_c_MS/3 MS.cluster_local.b_c_MS/3,MS.cluster_local.h_c_MS/3])+pos_MS_new;
        end
    end
    MS_new.mpc_local.pos_BS(replace_mpc,:) = new_mpc;%Copy out-range MPCs
    MS_new.mpc_local.pos_MS(replace_mpc,:) = new_mpc;
end

MS_new.cluster_local.p_c_MS = pos_MS_new;
MS_new.cluster_local.p_c_BS = pos_MS_new;
MS_new.pos = pos_MS_new;