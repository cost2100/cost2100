function [paraEx, paraSt, link, env] = cost2100(network, scenario, freq, snapRate, snapNum, BSPosCenter, BSPosSpacing, BSPosNum, MSPos, MSVelo)
%COST2100 The main body of the COST 2100 channel model implementation
%
%Default call:
%[paraEx, paraSt, link, env] = cost2100(network, scenario, freq, snapRate, 
%snapNum, BSPosCenter, BSPosSpacing, BSPosNum, MSPos, MSVelo)
%
%------
%Input:
%------
%network: 'IndoorHall_5GHz','SemiUrban_300MHz','SemiUrban_CloselySpacedUser_2_6GHz', 'Indoor_CloselySpacedUser_2_6GHz', or 'SemiUrban_VLA_2_6GHz'
%scenario: 'LOS' or 'NLOS'   
%freq: Frequency band [Hz]
%snapRate: Number of channel snapshot per s
%snapNum: Number of simulated snapshots
%BSPosCenter: Center position of BS array [x, y, z] [m]
%BSPosSpacing: Inter-position spacing, for large arrays [m]
%BSPosNum: Number of positions at each BS site, for large arrays
%MSPos: Initial position matrix of MSs, size = (number of MS, [x y z]) [m]
%MSVelo: Velocity vector of MSs, size = (number of MS, [x y z]) [m/s]
%
%------
%Output:
%------
%paraEx: External parameters
%paraSt: Stochastic parameters
%link(numBS,numMS): Link information matrix
% .channel: Channel information
% .MS: MS information
% .BS: BS information
%env: Radio propagation environment record
% .VRtable: VR assignment table 
% .MS_VR: MS visibility region
% .BS_VR: BS visibility region
% .BS_VR_len: BS visibility region length
% .BS_VR_slope: BS visibility region slope
% .cluster: Cluster information
% .mpc: MPC information
% .dmc: DMC information
%
%For more information about the cost2100 model, references and version history, 
%see 'Readme.txt'.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(freq)~=2
    error('Please specify freq as [freq_start freq_stop].\n');
end

chkVal = size(BSPosCenter);

if chkVal(2)~=3
    error('Please specify each BS position as [x y z].\n');
end

chkVal = size(MSPos);

if chkVal(2)~=3
    error('Please specify each MS position as [x y z].\n');
end

if ~all(size(MSPos)==size(MSVelo))
    error('Please check MSPos and MSVelo.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numBS = length(BSPosCenter(:,1)); % Number of BSs
numMS = length(MSPos(:,1)); % Number of MSs

[paraEx, paraSt] = get_para(network, scenario, freq, snapRate, snapNum, BSPosCenter, BSPosSpacing, BSPosNum, MSPos, MSVelo);
VRtable = get_VRtable(paraEx,paraSt); % Get VR table
MS_VR = get_MS_VR(VRtable, paraEx); % Get MS-VR locations
PRUNE_VR = 0; % Activate for speed if there are lots of useless MS-VRs.
if PRUNE_VR; [VRtable, MS_VR] = prune_MS_VR(MS_VR, VRtable, paraEx, paraSt); end
cluster = get_cluster(MS_VR, VRtable, paraEx, paraSt);
mpc = get_mpc(cluster, paraSt); % Get MPCs - scattering points in each cluster
dmc = get_dmc(cluster, paraSt); % Get DMCs - diffuse multipath components
[BS_VR, BS_VR_len, BS_VR_slope] = get_BS_VR_para(VRtable, paraEx, paraSt); % Get BS_VR parameters (large array extension)

% For every BS
for m = 1:numBS 
    BS(m).idx = m; % Label of BS(m)
    BS(m).pos = BSPosCenter(m,:); % Position of BS
    BS(m).VRLOS = get_VRLOS(BS(m), paraEx, paraSt); % LOS VR at BS
    BS(m).cluster_local = get_cluster_local('BS', m, paraEx, paraSt); % BS local cluster
    BS(m).mpc_local = get_mpc( BS(m).cluster_local, paraSt); % MPC in BS local cluster
    BS(m).dmc_local = get_dmc( BS(m).cluster_local, paraSt); % DMC in BS local cluster
       
    % Position of the antennas within one BS, for very large arrays
    pos_vec = (-1*BSPosSpacing(m, 1)/2-BSPosSpacing(m, 1)*(BSPosNum(m)/2-1)):BSPosSpacing(m, 1):(BSPosSpacing(m, 1)/2+BSPosSpacing(m,1)*(BSPosNum(m)/2-1));
    
    if ~isempty(pos_vec)
        BS(m).BS_pos(:, 1) = BS(m).pos(1)+pos_vec; % [x]
    else
        BS(m).BS_pos(:, 1) = BS(m).pos(1);
    end
    
    BS(m).BS_pos(:, 2) = BS(m).pos(2); % [y]
    BS(m).BS_pos(:, 3) = BS(m).pos(3); % [z]
end

% For every MS
for m = 1:numMS
    MS(m).idx = m; % Label of MS(m)       
    MS(m).pos = MSPos(m,:); % Position of MS
    MS(m).velo = MSVelo(m,:); % Velocity of MS
    
    MS(m).cluster_local = get_cluster_local('MS', m, paraEx, paraSt); % MS local cluster
    MS(m).mpc_local = get_mpc( MS(m).cluster_local, paraSt); % MPC in MS local cluster
    MS(m).dmc_local = get_dmc( MS(m).cluster_local, paraSt); % DMC in MS local cluster          
end 


% Get the channel of each link
for nSS = 1:snapNum  
    for nB = 1:numBS % Loop for every BS
        for nB_pos = 1:BSPosNum(nB) % Loop for every position at the same BS
            for nM = 1:numMS % Loop for every MS  
                link(nB, nM).channel{nB_pos, nSS} = get_channel(BS(nB), BS.BS_pos(nB_pos, :), MS(nM), VRtable, MS_VR, BS_VR, BS_VR_len, BS_VR_slope,...
                                                                cluster, mpc, dmc, paraEx, paraSt);
                link(nB, nM).MS(nB_pos, nSS) = MS(nM); % Record of MS at snapshot m
                link(nB, nM).BS(nB_pos, nSS) = BS(nB); % Record of BS at snapshot m
            end
        end    
    end
   
    % Update MS information according to movement
    for nM = 1:numMS
        MS(nM) = update_chan(MS(nM), paraEx, paraSt);
    end
end

% Saving the environment information
env.VRtable = VRtable;
env.MS_VR = MS_VR;
env.BS_VR = BS_VR; 
env.BS_VR_len = BS_VR_len; 
env.BS_VR_slope = BS_VR_slope; 
env.cluster = cluster; % Far clusters
env.mpc = mpc;
env.dmc = dmc;

function [VRtable, MS_VR] = prune_MS_VR(MS_VR, VRtable, paraEx, paraSt)
prune_list = ones(size(MS_VR,1),1);
for MS_VR_idx = 1:size(MS_VR,1)
    for MS_idx = 1:paraEx.num_MS
        velo_MS = paraEx.velo_MS(MS_idx,1:2);
        pos_MS = paraEx.pos_MS(MS_idx,1:2);
        if (norm(velo_MS)==0)
            t_star = 0;
        else
            t_star = - 1/norm(velo_MS)^2 * velo_MS*(pos_MS - MS_VR(MS_VR_idx,1:2)).';
        end
        if (t_star<0)
            t_min = 0;
        elseif (t_star>(paraEx.snap_num-1)/paraEx.snap_rate)
            t_min = (paraEx.snap_num-1)/paraEx.snap_rate;
        else
            t_min = t_star;
        end
        d_min = norm(pos_MS - MS_VR(MS_VR_idx,1:2) + velo_MS*t_min);
        if (d_min<paraSt.r_c)
            prune_list(MS_VR_idx) = 0;
        end
    end
end
VRtable(:,prune_list==1,:) = [];
MS_VR(prune_list==1,:) = [];
% MS-VR indices should be contiguous.
MS_VR_indices = VRtable(:,:,2);
[C,~,IC] = unique(MS_VR_indices(:));
MS_VR_indices_min = 1:length(C);
VRtable(:,:,2) = MS_VR_indices_min(IC);
