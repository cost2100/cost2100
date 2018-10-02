function channel = get_channel(BS, BS_pos, MS, VRtable, MS_VR, BS_VR,...
    BS_VR_len, BS_VR_slope, cluster, mpc, dmc, paraEx, paraSt)
%GET_CHANNEL Get the channel of a single BS-MS link
%
%Default call:
%function channel = get_channel(BS ,BS_pos, MS, VRtable, MS_VR, BS_VR,...
%    BS_VR_len, BS_VR_slope, cluster, mpc, dmc, paraEx, paraSt)
%
%-------
%Input:
%-------
%BS: BS information
%BS_pos: BS position within a VLA setting, to be processed
%MS: MS information
%VRtable: VR allocation table
%MS_VR: MS-VRs locations
%BS_VR: BS-VRs locations
%BS_VR_len: Length of the BS-VRs.
%BS_VR_slope: Power slopes of the BS-VRs
%cluster: Cluster information
%mpc: All MPCs of the far clusters
%dmc: Dense multipath components
%paraEx: External parameters
%paraSt: Stochastic parameters
%------
%Output:
%------
%channel: Channel recording
% .a_c: Cluster power attenuation
% .h(numAllMPC,[AoD EoD AoA EoA delay amp. vv vh hh hv]): MPCs
% .pathloss: Channel pathloss
% .active_c: Indices of active far clusters
% .a_VR: VR attenuation
% .h_los: LOS channel response

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = paraEx.c0/paraEx.center_freq; % Wavelength 
pathloss = calc_pathloss(BS.pos,MS.pos,paraEx); % Compute the pathloss 

d_MS_BS = calc_dist(BS.pos, MS.pos); % Distance between MS and BS
tau_0 = d_MS_BS/paraEx.c0; % Delay of LOS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the cluster activity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[active_MS_VR, active_BS_VR, active_c] = get_active_MS_VR(BS, MS, VRtable, MS_VR, BS_VR, BS_VR_len, paraSt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-allocated data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_VR = [];

if ~isempty(active_MS_VR) % If there are active far cluster

    % MS-VR gain
    gain_MS_VR = zeros(1, length(active_MS_VR));
    for m = 1:length(active_MS_VR)
        d_MS_VR = calc_dist(MS.pos(1:2), MS_VR(active_MS_VR(m),:));
        y = d_MS_VR - (paraSt.r_c-paraSt.l_c);
        gain_MS_VR(m) = 1/2-atan(2*sqrt(2)*y/sqrt(lambda*paraSt.l_c))/pi;
    end

    % A_VR for cluster.    
    % Now, it can happen that several MS VRs are associated with the same
    % cluster, if `MS VR groups' of size larger than one are used. (This
    % situation is not desirable since, for a given MS, it would mean the
    % VRs are overlapping, and contributions from the cluster are counted
    % twice!) We deal with it by defining a `effective VR MS', which is the
    % one with the maximum gain
    a_MS_VR = zeros(1, length(active_c));
    a_MS_VR_idx = zeros(1, length(active_c));
    for m = 1:length(active_c)
        active_MS_VR_cluster_m_bitmap = VRtable(BS.idx,active_MS_VR,2)==active_c(m);
        [a_MS_VR(m),a_MS_VR_idx(m)] = max(gain_MS_VR.*active_MS_VR_cluster_m_bitmap);
    end
    
    % To all effects, in what follows we can restrict ourselves to
    % active_MS_VR_eff, and forget about active_MS_VR
    active_MS_VR_eff = active_MS_VR(a_MS_VR_idx);
    active_MS_VR = active_MS_VR_eff;
    
    % BS-VR gain
    a_BS_VR = zeros(1, length(active_c));
    for m = 1:length(active_c)
        d_BS_VR = calc_dist(BS.pos(1:2),BS_VR(active_c(m),:));
        a_BS_VR(m) = sqrt(10^(d_BS_VR*BS_VR_slope(active_c(m))/10)); % Convert from dB to linear, and then to amplitude
    end
    
    % A_VR including BS-BR gain
    a_VR = a_MS_VR.*a_BS_VR;

    h_far = zeros(paraSt.n_mpc, length(active_c), 10); % Store all MPC information
    
    % Loop over active clusters
    for m = 1:length(active_c) 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine the cluster power attenuation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d_BS_c_MS = calc_dist(BS.pos, cluster(active_c(m)).pos_c_BS)+calc_dist(MS.pos, cluster(active_c(m)).pos_c_MS); % Propagation distance (BS-cluster + cluster-MS)
        tau_m = d_BS_c_MS/paraEx.c0+cluster(active_c(m)).tau_c_link; % Cluster delay
        a_c = max(exp((-paraSt.k_tau*(tau_m-tau_0)*1e6)), exp((-paraSt.k_tau*(paraSt.tau_b-tau_0)*1e6))); % Attenuation of clusters

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create the channels
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        direct_mpc_BS = mpc(active_c(m)).pos_BS - repmat(BS.pos,paraSt.n_mpc,1); % Vector from MPC to BS
        [phi_mpc_BS, theta_mpc_BS, r_mpc_BS] = cart2sph(direct_mpc_BS(:,1), direct_mpc_BS(:,2), direct_mpc_BS(:,3));
        direct_mpc_MS = squeeze(mpc(active_c(m)).pos_MS) - repmat(MS.pos,paraSt.n_mpc,1); % Vector from MPC to MS
        [phi_mpc_MS, theta_mpc_MS, r_mpc_MS] = cart2sph(direct_mpc_MS(:,1), direct_mpc_MS(:,2), direct_mpc_MS(:,3));
        tau_mpc = (r_mpc_BS+r_mpc_MS)/paraEx.c0 + cluster(active_c(m)).tau_c_link;
            
        % MPC gain function.
        mpc_gain_pos = mpc(active_c(m)).gain_center + repmat(MS_VR(active_MS_VR(m),:),paraSt.n_mpc,1); % Center location of MPC gain function [x, y]
                        
        % Two dimensional Gauss function
        sigma_x = mpc(active_c(m)).gain_radius;
        sigma_y = mpc(active_c(m)).gain_radius;
        a_mpc_gain = exp(-1*((mpc_gain_pos(:,1)-MS.pos(1)).^2./(2*sigma_x.^2)+(mpc_gain_pos(:,2)-MS.pos(2)).^2./(2*sigma_y.^2)));

        % Magnitude of the MPC amplitude.
        % Note that pathloss has amplitude `dimension'
        amp_h = a_VR(m)*sqrt( a_c*cluster(active_c(m)).shadow_f )*a_mpc_gain.*(mpc(active_c(m)).a_mpc).'*pathloss;
        %           |  distance from peak
        %      X    |  ------------------
        %           |      gain radius
        % ----------------------------------
        %  -3.00 dB |         0.83
        %  -4.34 dB |         1.00
        % -10.00 dB |         1.52
        % -20.00 dB |         2.15
        % -40.00 dB |         3.03
        % -50.00 dB |         3.39
        % -60.00 dB |         3.71
        % -80.00 dB |         4.29
        weak_mpcs = pow2db(abs(a_mpc_gain).^2) < -10;
        amp_h(weak_mpcs) = 0;
        
        % Porization matrix.
        amp_vv = mpc(active_c(m)).vv.';
        amp_vh = mpc(active_c(m)).vh.';
        amp_hh = mpc(active_c(m)).hh.';
        amp_hv = mpc(active_c(m)).hv.';

        % The unfiltered impulse response, angle in radian.
        h_far(:, m, :) = [phi_mpc_BS, theta_mpc_BS, phi_mpc_MS, theta_mpc_MS, tau_mpc, amp_h, amp_vv, amp_vh, amp_hh, amp_hv];
        
        % Record the cluster power attenuation.
        channel.a_c(m) = a_c; 
    end
    
    h = reshape(h_far, [], 10);
else
    h = [];
end

% Local cluster at BS side
if BS.cluster_local.active 
    
    h_local_BS = zeros(paraSt.n_mpc, 10);
       
    for n = 1:paraSt.n_mpc
        if (norm(MS.mpc_local.pos_BS(n,:) - MS.mpc_local.pos_MS(n,:))>0)
            error('Local clusters are single clusters.');
        end
        direct_mpc_BS = squeeze(BS.mpc_local.pos_BS(n,:)) - BS.pos; % Vector from MPC to BS
        [phi_mpc_BS, theta_mpc_BS, r_mpc_BS] = cart2sph(direct_mpc_BS(1), direct_mpc_BS(2), direct_mpc_BS(3));
        direct_mpc_MS = squeeze(BS.mpc_local.pos_MS(n,:)) - MS.pos; % Vector from MPC to MS
        [phi_mpc_MS, theta_mpc_MS, r_mpc_MS] = cart2sph(direct_mpc_MS(1), direct_mpc_MS(2), direct_mpc_MS(3));
        tau_mpc = (r_mpc_BS+r_mpc_MS)/paraEx.c0 + BS.cluster_local.tau_c_link;

        % Magnitude of the MPC amplitude
        amp_h = sqrt( BS.cluster_local.shadow_f )*BS.mpc_local.a_mpc(n)*pathloss;
        
        % Polarization matrix
        amp_vv = mpc(active_c(m)).vv(n);
        amp_vh = mpc(active_c(m)).vh(n);
        amp_hh = mpc(active_c(m)).hh(n);
        amp_hv = mpc(active_c(m)).hv(n);
        
         % The unfiltered impulse response, angle in radian
        h_local_BS(n, :) = [phi_mpc_BS theta_mpc_BS phi_mpc_MS theta_mpc_MS tau_mpc amp_h amp_vv amp_vh amp_hh amp_hv];
    end
else
    h_local_BS = [];
end

% Local cluster at MS side
if MS.cluster_local.active 
    
    h_local_MS = zeros(paraSt.n_mpc, 10);
      
    for n = 1:paraSt.n_mpc
        if (norm(MS.mpc_local.pos_BS(n,:) - MS.mpc_local.pos_MS(n,:))>0)
            error('Local clusters are single clusters.');
        end
        direct_mpc_BS = squeeze(MS.mpc_local.pos_BS(n,:)) - BS.pos; % Vector from MPC to BS
        [phi_mpc_BS, theta_mpc_BS, r_mpc_BS] = cart2sph(direct_mpc_BS(1),direct_mpc_BS(2),direct_mpc_BS(3));
        direct_mpc_MS = squeeze(MS.mpc_local.pos_MS(n,:)) - MS.pos; % Vector from MPC to MS
        [phi_mpc_MS, theta_mpc_MS, r_mpc_MS] = cart2sph(direct_mpc_MS(1),direct_mpc_MS(2),direct_mpc_MS(3));
        tau_mpc = (r_mpc_BS+r_mpc_MS)/paraEx.c0 + MS.cluster_local.tau_c_link;

        % MPC gain function
        if (norm(MS.mpc_local.pos_BS(n,1:2) - MS.mpc_local.gain_center(n,:))>0)
            error('With local clusters, the cluster center is the gain center.');
        end
        mpc_gain_pos = squeeze(MS.mpc_local.gain_center(n, :, :)); % Center location of MPC gain function [x, y]

        % Two dimensional Gauss function
        sigma_x = MS.mpc_local.gain_radius(n);
        sigma_y = MS.mpc_local.gain_radius(n);
        a_mpc_gain = exp(-1*((mpc_gain_pos(1)-MS.pos(1))^2/(2*sigma_x^2)+(mpc_gain_pos(2)-MS.pos(2))^2/(2*sigma_y^2)));
        
        if pow2db(abs(a_mpc_gain)^2) < -10
        %           |  distance from peak
        %      X    |  ------------------
        %           |      gain radius
        % ----------------------------------
        %  -3.00 dB |         0.83
        %  -4.34 dB |         1.00
        % -10.00 dB |         1.52
        % -20.00 dB |         2.15
        % -40.00 dB |         3.03
        % -50.00 dB |         3.39
        % -60.00 dB |         3.71
        % -80.00 dB |         4.29
            amp_h = 0;
        else
            % Magnitude of the MPC amplitude
            amp_h = sqrt( MS.cluster_local.shadow_f )*MS.mpc_local.a_mpc(n)*a_mpc_gain*pathloss;
        end
        
        % Polarization matrix
        amp_vv = mpc(active_c(m)).vv(n);
        amp_vh = mpc(active_c(m)).vh(n);
        amp_hh = mpc(active_c(m)).hh(n);
        amp_hv = mpc(active_c(m)).hv(n);

        % The unfiltered impulse response, angle in radians
        h_local_MS(n, :) = [phi_mpc_BS theta_mpc_BS phi_mpc_MS theta_mpc_MS tau_mpc amp_h amp_vv amp_vh amp_hh amp_hv]; 
    end
else
    h_local_MS = [];
end

% Contains information of all MPCs, except for the LOS
h = [h; h_local_BS; h_local_MS]; 

if isempty(h)
    error('Empty h, cluster');
end

% NOTE: The following collides with code in visual_channel.m
% Only keep MPCs with non-zero amplitude
h = h(abs(h(:, 6))>0, :);

if isempty(h)
    error('Empty h, mpc');
end

channel.h = h; % H record
channel.pathloss = pathloss; % Pathloss record
channel.active_c = active_c; % Index of active clusters record
channel.a_VR = a_VR; % VR attenuation record
channel.h_los = get_channel_los(channel, BS, BS_pos, MS, paraEx,paraSt); % LOS

function [active_MS_VR_Idx,active_BS_VR_Idx,active_c_Idx] = get_active_MS_VR(BS,MS,VRtable,MS_VR,BS_VR,BS_VR_len,paraSt)
% We return 3 sets, active_MS_VR_Idx, active_BS_VR_Idx, and active_c_Idx,
% which are the set of active MS-VRs, the set of active BS-VRS, and the set
% of active (far) clusters. For the current version of the COST 2100 model,
% we have that (active_BS_VR_Idx == active_c_Idx) holds. In general,
% however, active_MS_VR_Idx may be larger than active_c_Idx; this happens
% when `MS VR groups' of size larger than one are used.

% Find which MS-VRs are, via interaction objects (IOs), reachable from BS
MS_VR_Idx = find(VRtable(BS.idx, :, 1)==1).';

% Only MS-VRs sufficiently close to the MS current position are active
d_MS_VR = sqrt(sum((repmat(MS.pos(1:2),length(MS_VR_Idx),1) - MS_VR).^2,2));
active_MS_VR_Idx = MS_VR_Idx(d_MS_VR < paraSt.r_c);

% NOTE: Because currently BS-VRs only support one BS, all BS VRs belong to
% it. This might need to be changed when support for multiple BSs is added. 

% Only BS-VRs sufficiently close to the BS current array section are
% active (for compact arrays we have infinite BR_VR_len) 
d_BS_VR = sqrt(sum((repmat(BS.pos(1:2),size(BS_VR, 1),1) - BS_VR).^2,2));
active_BS_VR_Idx = find(d_BS_VR < BS_VR_len/2);

% A given MS-VR is active if a BS-VR can be found which illuminates the
% same cluster, and vice versa
active_MS_c_Idx = VRtable(BS.idx, active_MS_VR_Idx, 2);
% Note the simple mapping BS-VR to illuminated IO. This is because no
% 'BS VR groups' have been introduced
active_BS_c_Idx = active_BS_VR_Idx;
active_c_Idx = intersect(active_MS_c_Idx,active_BS_c_Idx);
if isempty(active_c_Idx)
    error('No active cluster! Please restart simulation.');
end

% Last, find the inverse mappings to map active clusters to active MS VRs, and
% active clusters to active BS VRs
active_BS_VR_Idx = active_c_Idx;
% The latter is a bit involved, since the forward mapping MS-VR to cluster
% is not injective (one-to-one). In particular, active_MS_VR_Idx might be
% longer than active_c_Idx if 'MS VR groups' of size larger than one are
% used. This will be handled by the parent routine by defining a
% 'reference' MS-VR.
active_MS_VR_Idx = active_MS_VR_Idx(ismember(active_MS_c_Idx,active_c_Idx));
