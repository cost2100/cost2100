function channel = get_channel( BS, MS, VRtable, VR, cluster, mpc, dmc, paraEx, paraSt)
%GET_CHANNEL get the channel in link BS-MS
%Default call: channel = get_channel( BS, MS, VRtable, VR, cluster, mpc,
%dmc, paraEx, paraSt)
%-------
%Input:
%-------
%BS: BS information
%MS: MS information
%paraEx, paraSt: external and stochastic parameters
%------
%Output:
%------
%channel: channel recording, struct.
% .h(numAllMPC,[AoD EoD AoA EoA delay amp.]): unfiltered channel response
% .pathloss: channel pathloss
% .active_c: index of active far clusters
% .a_VR: VR attenuation 
% .h_los: LOS channel reponse
%
%See also: cost2100, get_para, get_VR, get_cluster, get_mpc, calc_pathloss,
%get_channel_los

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_MS_BS = calc_dist(BS.pos, MS.pos); %Distance between MS and BS
lambda = paraEx.c0/paraEx.freq; %Wavelength 
pathloss = calc_pathloss(BS.pos,MS.pos,paraEx);%Compute the pathloss 
tau_0=d_MS_BS/paraEx.c0; %Delay of LOS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine the cluster activity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VRBS = VR(find(VRtable(BS.idx,:,1)==1),:); %VRs belong to BS
VRIdx = find(VRtable(BS.idx,:,1)==1); %VRs index
numVRBS = length(VRBS); %Number of VRs belong to BS
% save VRBS VRBS 
d_MS_VR = zeros(1,numVRBS);%Distance MS-VRs
for m = 1:numVRBS
    d_MS_VR(m) = calc_dist(MS.pos(1:2),VRBS(m,:));
end
% save d_MS_VR d_MS_VR 
active_VR = VRIdx(find(d_MS_VR < paraSt.r_c)); %Index of active VR
% save active_VR active_VR
% active_VR
[tmpVal active_VR_ref] = find(d_MS_VR< paraSt.r_c); %Reference index of the active VR
active_c = VRtable(BS.idx,active_VR,2); %Index of active clusters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-allocated data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_actc = 0;
a_VR = [];
h = [];

if length(active_VR)>0 %If there are active far cluster
    active_c = active_c([true diff(sort(active_c))>0]); %Index of active clusters
 
    for m = 1:length(active_VR)
        y = paraSt.l_c+d_MS_VR(active_VR_ref(m))-paraSt.r_c;
        gainVR(m) = 1/2-atan(2*sqrt(2)*y/sqrt(lambda*paraSt.l_c))/pi; %VR gain
    end

    for m = 1:length(active_c)
        a_VR(m) = max(gainVR( (find(VRtable(BS.idx,active_VR,2)==active_c(m)) ))); %A_VR for cluster
    end

    for m = 1:length(active_c)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Determine the cluster power attenuation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d_BS_c_MS = calc_dist(BS.pos,cluster(active_c(m)).pos_c_BS)+calc_dist(MS.pos,cluster(active_c(m)).pos_c_MS); %Propagation distance
        tau_m = d_BS_c_MS/paraEx.c0+cluster(active_c(m)).tau_c_link; %Cluster delay
        a_c=max(exp((-paraSt.k_tau*(tau_m-tau_0)*1e6)),exp((-paraSt.k_tau*(paraSt.tau_b-tau_0)*1e6))); % Attenuation of clusters

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Create the channels
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for n = 1:paraSt.n_mpc
            direct_mpc_BS = squeeze(mpc(active_c(m)).pos_BS(n,:)) - BS.pos; %Vector from mpc to BS;
            [phi_mpc_BS theta_mpc_BS r_mpc_BS] = cart2sph(direct_mpc_BS(1),direct_mpc_BS(2),direct_mpc_BS(3));
            direct_mpc_MS = squeeze(mpc(active_c(m)).pos_MS(n,:)) - MS.pos; %Vector from mpc to MS;
            [phi_mpc_MS theta_mpc_MS r_mpc_MS] = cart2sph(direct_mpc_MS(1),direct_mpc_MS(2),direct_mpc_MS(3));
            tau_mpc = (r_mpc_BS+r_mpc_MS)/paraEx.c0 + cluster(active_c(m)).tau_c_link;
            amp_h = a_VR(m)*sqrt( cluster(active_c(m)).shadow_f*a_c )*mpc(active_c(m)).a_mpc(n)*pathloss*exp(-1j*2*pi*paraEx.freq*tau_mpc);
            %Original%Transition*sqrt(shadow fading*cluster attenuation*mpc attenuation)*pathloss*mpc phase

            h( (m-1)*paraSt.n_mpc+1+n,: ) = [phi_mpc_BS theta_mpc_BS phi_mpc_MS theta_mpc_MS tau_mpc amp_h]; %The unfiletered impulse response, angle in radii
        end
        channel.a_c(m) = a_c; %The cluster power attenuation
    end

    num_actc = length(active_c);
end

%Get the MPCs from local clusters
if BS.cluster_local.active
    for n = 1:paraSt.n_mpc
        direct_mpc_BS = squeeze(BS.mpc_local.pos_BS(n,:)) - BS.pos; %Vector from mpc to BS;
        [phi_mpc_BS theta_mpc_BS r_mpc_BS] = cart2sph(direct_mpc_BS(1),direct_mpc_BS(2),direct_mpc_BS(3));
        direct_mpc_MS = squeeze(BS.mpc_local.pos_MS(n,:)) - MS.pos; %Vector from mpc to MS;
        [phi_mpc_MS theta_mpc_MS r_mpc_MS] = cart2sph(direct_mpc_MS(1),direct_mpc_MS(2),direct_mpc_MS(3));
        tau_mpc = (r_mpc_BS+r_mpc_MS)/paraEx.c0 + BS.cluster_local.tau_c_link;
        
        amp_h = sqrt( BS.cluster_local.shadow_f )*BS.mpc_local.a_mpc(n)*pathloss*exp(-1j*2*pi*paraEx.freq*tau_mpc);

        h( num_actc*paraSt.n_mpc+1+n,: ) = [phi_mpc_BS theta_mpc_BS phi_mpc_MS theta_mpc_MS tau_mpc amp_h]; %The unfiletered impulse response, angle in radii
    end
    num_actc = num_actc+1;
end

if MS.cluster_local.active
    for n = 1:paraSt.n_mpc
        direct_mpc_BS = squeeze(MS.mpc_local.pos_BS(n,:)) - BS.pos; %Vector from mpc to BS;
        [phi_mpc_BS theta_mpc_BS r_mpc_BS] = cart2sph(direct_mpc_BS(1),direct_mpc_BS(2),direct_mpc_BS(3));
        direct_mpc_MS = squeeze(MS.mpc_local.pos_MS(n,:)) - MS.pos; %Vector from mpc to MS;
        [phi_mpc_MS theta_mpc_MS r_mpc_MS] = cart2sph(direct_mpc_MS(1),direct_mpc_MS(2),direct_mpc_MS(3));
        tau_mpc = (r_mpc_BS+r_mpc_MS)/paraEx.c0 + MS.cluster_local.tau_c_link;

        amp_h = sqrt( MS.cluster_local.shadow_f )*MS.mpc_local.a_mpc(n)*pathloss*exp(-1j*2*pi*mod(paraEx.freq*tau_mpc,1));

        h( num_actc*paraSt.n_mpc+1+n,: ) = [phi_mpc_BS theta_mpc_BS phi_mpc_MS theta_mpc_MS tau_mpc amp_h]; %The unfiletered impulse response, angle in radii
    end
end

channel.h = h; %H record
channel.pathloss = pathloss; %Pathloss record
channel.active_c = active_c; %Index of active clusters record
channel.a_VR = a_VR; %VR attenuation record
channel.h_los = get_channel_los(channel, BS, MS, paraEx,paraSt); %LOS 